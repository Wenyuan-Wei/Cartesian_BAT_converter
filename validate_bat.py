#!/usr/bin/env python3
"""validate_bat.py — Independent Python cross-validation of BAT coordinate values.

This script proves that every value stored in a .bat file is geometrically
correct with respect to the original PDB atom positions, using only the Python
standard library (math module).

Strategy
--------
1. Read the original PDB file, replicating the same parsing rules as the C++
   PDBReader::read() — ATOM records only, chain filter, altLoc/occupancy
   deduplication — to get atom positions in the same 0-based index order.
2. Read the .bat file to extract anchor indices, tree refs (p, g, h), and the
   stored (bond, angle, torsion) values for each entry.
3. Apply the virtual-reference resolution convention (same as compute_bat /
   reconstruct_cartesian):
     • missing g  (raw -1) → anchor.N_idx
     • missing h  (raw -1) → anchor.C_idx, unless C == p or C == g,
                             in which case anchor.N_idx
4. For each entry, independently compute bond / angle / torsion from PDB
   positions and compare to the stored values.
5. Report statistics and exit 0 on PASS, 1 on FAIL.

The dihedral formula matches geometry::dihedral(a, b, c, d) exactly:
  b1 = b-a, b2 = c-b, b3 = d-c
  n1 = cross(b1, b2), n2 = cross(b2, b3)
  x = dot(n1, n2)
  y = dot(cross(n1, n2), normalise(b2))
  return atan2(y, x)

Usage
-----
  python3 validate_bat.py <input.pdb> <input.bat> [chain]

Exit codes
----------
  0  All entries match to within tolerance.
  1  One or more mismatches, or a parsing / I/O error.
"""

import sys
import math
import collections

# ---------------------------------------------------------------------------
# Tolerance — should be well below floating-point rounding noise for the
# same arithmetic performed in Python vs C++ (both IEEE 754 double).
# Both programs parse identical ASCII PDB coordinates, so the only source
# of difference is multi-step floating-point rounding order (< 1e-12).
# 1e-8 is therefore comfortably tight while being immune to platform noise.
# ---------------------------------------------------------------------------
TOLERANCE = 1e-8

# ---------------------------------------------------------------------------
# Pure-Python geometry (mirrors geometry.cpp exactly)
# ---------------------------------------------------------------------------

def _dist(a, b):
    dx, dy, dz = a[0]-b[0], a[1]-b[1], a[2]-b[2]
    return math.sqrt(dx*dx + dy*dy + dz*dz)


def _dot(a, b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]


def _cross(a, b):
    return (a[1]*b[2] - a[2]*b[1],
            a[2]*b[0] - a[0]*b[2],
            a[0]*b[1] - a[1]*b[0])


def _norm(v):
    return math.sqrt(_dot(v, v))


def _angle(a, b, c):
    """geometry::angle(a, b, c): b is the vertex atom."""
    ba = (a[0]-b[0], a[1]-b[1], a[2]-b[2])
    bc = (c[0]-b[0], c[1]-b[1], c[2]-b[2])
    nba = _norm(ba)
    nbc = _norm(bc)
    if nba == 0.0 or nbc == 0.0:
        return 0.0
    cos_t = _dot(ba, bc) / (nba * nbc)
    return math.acos(max(-1.0, min(1.0, cos_t)))


def _dihedral(a, b, c, d):
    """geometry::dihedral(a, b, c, d): dihedral around b-c bond."""
    b1 = (b[0]-a[0], b[1]-a[1], b[2]-a[2])
    b2 = (c[0]-b[0], c[1]-b[1], c[2]-b[2])
    b3 = (d[0]-c[0], d[1]-c[1], d[2]-c[2])

    n1 = _cross(b1, b2)
    n2 = _cross(b2, b3)

    b2n = _norm(b2)
    if b2n == 0.0:
        return 0.0
    b2h = (b2[0]/b2n, b2[1]/b2n, b2[2]/b2n)

    x = _dot(n1, n2)
    y = _dot(_cross(n1, n2), b2h)
    return math.atan2(y, x)


# ---------------------------------------------------------------------------
# PDB reader — replicates PDBReader::read() parsing rules
# ---------------------------------------------------------------------------

_ALLOWED_ELEMENTS = {'H', 'C', 'N', 'O', 'S'}

_AtomKey = collections.namedtuple('_AtomKey', ['chain', 'resid', 'name'])


def _altloc_rank(alt):
    """Lower rank is preferred (mirrors C++ altloc_rank)."""
    if alt in (' ', '\x00', ''):
        return 0
    if alt == 'A':
        return 1
    return 2


def _infer_element(raw4):
    """Infer element from the 4-char PDB atom-name field (columns 13–16)."""
    if len(raw4) < 2:
        return None
    if raw4[0].isalpha():
        return raw4[0].upper()
    if raw4[1].isalpha():
        return raw4[1].upper()
    return None


def read_pdb_atoms(path, chain_filter=None):
    """
    Parse ATOM records from *path*, deduplicate by (chain, resid, name)
    keeping the best alt-loc / occupancy, then optionally filter by chain.

    Returns a list of {'name', 'resname', 'chain_id', 'resid', 'pos', 'element'}
    dicts in the same 0-based order as C++ Structure::atoms.

    Raises ValueError for unsupported elements (mirroring C++ behaviour).
    """
    seen   = {}   # _AtomKey → index in `atoms`
    atoms  = []

    with open(path) as fh:
        for line in fh:
            if len(line) < 54:
                continue
            if line[:6] != 'ATOM  ':
                continue

            try:
                raw_name  = line[12:16]               # 4-char, NOT stripped
                name      = raw_name.strip()
                alt_loc   = line[16] if len(line) > 16 else ' '
                resname   = line[17:20].strip()
                chain_id  = line[21] if len(line) > 21 else ' '
                resid     = int(line[22:26])
                x         = float(line[30:38])
                y         = float(line[38:46])
                z         = float(line[46:54])
                occ       = float(line[54:60]) if len(line) >= 60 else 1.0
            except (ValueError, IndexError):
                continue

            # Element column (columns 77–78, 0-based 76–78)
            element = ''
            if len(line) >= 78:
                element = line[76:78].strip().upper()
            if not element:
                element = _infer_element(raw_name) or ''
            else:
                element = element.upper()

            if element not in _ALLOWED_ELEMENTS:
                raise ValueError(
                    f"Unsupported element '{element}' in "
                    f"{resname} {chain_id}{resid} {name}"
                )

            atom = dict(name=name, resname=resname, chain_id=chain_id,
                        resid=resid, pos=(x, y, z), element=element,
                        alt_loc=alt_loc, occ=occ)

            key = _AtomKey(chain=chain_id, resid=resid, name=name)
            if key not in seen:
                seen[key] = len(atoms)
                atoms.append(atom)
            else:
                existing = atoms[seen[key]]
                if atom['occ'] > existing['occ']:
                    atoms[seen[key]] = atom
                elif atom['occ'] == existing['occ']:
                    if _altloc_rank(atom['alt_loc']) < _altloc_rank(existing['alt_loc']):
                        atoms[seen[key]] = atom
                # else: keep existing (lower occ → discard)

    # Chain filter (mirrors Structure::select_chain)
    if chain_filter:
        atoms = [a for a in atoms if a['chain_id'] == chain_filter]

    return atoms


# ---------------------------------------------------------------------------
# BAT file reader
# ---------------------------------------------------------------------------

def read_bat_file(path):
    """
    Parse a .bat file.

    Returns
    -------
    anchor : dict with keys N_idx, N_pos, CA_idx, CA_pos, C_idx, C_pos
    entries : list of dicts with keys atom_idx, p, g, h, bond, angle, torsion
              (p/g/h are raw values: -1 means BFSTree::npos)
    natoms  : int
    """
    anchor  = {}
    entries = []
    natoms  = None

    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            tag   = parts[0]

            if tag == 'BAT':
                if int(parts[1]) != 1:
                    raise ValueError(f"Unsupported BAT version: {parts[1]}")

            elif tag == 'NATOMS':
                natoms = int(parts[1])

            elif tag in ('ANCHOR_N', 'ANCHOR_CA', 'ANCHOR_C'):
                idx = int(parts[1])
                pos = (float(parts[2]), float(parts[3]), float(parts[4]))
                label = tag.split('_', 1)[1]   # 'N', 'CA', or 'C'
                anchor[f'{label}_idx'] = idx
                anchor[f'{label}_pos'] = pos

            elif tag == 'ENTRY':
                entries.append({
                    'atom_idx': int(parts[1]),
                    'p':        int(parts[2]),
                    'g':        int(parts[3]),
                    'h':        int(parts[4]),
                    'bond':     float(parts[5]),
                    'angle':    float(parts[6]),
                    'torsion':  float(parts[7]),
                })

    return anchor, entries, natoms


# ---------------------------------------------------------------------------
# Virtual-reference resolution (mirrors bat_coords.cpp / bat_reconstruct.cpp)
# ---------------------------------------------------------------------------

def resolve_refs(p_raw, g_raw, h_raw, aN, aC):
    """
    Apply the virtual-reference substitution used by compute_bat and
    reconstruct_cartesian.  -1 represents BFSTree::npos.

    Rules:
      g == -1  →  g = aN   (anchor N)
      h == -1  →  h = aC   if aC != p and aC != g
                  h = aN   otherwise
    """
    g = g_raw if g_raw != -1 else aN
    if h_raw == -1:
        h = aC if (aC != p_raw and aC != g) else aN
    else:
        h = h_raw
    return p_raw, g, h   # p is always valid for non-root atoms


# ---------------------------------------------------------------------------
# Torsion wrap-around helper
# ---------------------------------------------------------------------------

def _torsion_diff(a, b):
    """Absolute angular difference, accounting for ±π wrap-around."""
    d = abs(a - b)
    if d > math.pi:
        d = abs(d - 2.0 * math.pi)
    return d


# ---------------------------------------------------------------------------
# Main validation
# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <input.pdb> <input.bat> [chain]",
              file=sys.stderr)
        sys.exit(1)

    pdb_path   = sys.argv[1]
    bat_path   = sys.argv[2]
    chain_filt = sys.argv[3] if len(sys.argv) >= 4 else None

    print(f"=== validate_bat: {pdb_path}  ×  {bat_path} ===")
    if chain_filt:
        print(f"Chain filter      : {chain_filt}")

    # ── Read PDB ──────────────────────────────────────────────────────────
    try:
        atoms = read_pdb_atoms(pdb_path, chain_filt)
    except Exception as exc:
        print(f"FAIL  PDB read: {exc}", file=sys.stderr)
        sys.exit(1)

    print(f"PDB atoms loaded  : {len(atoms)}")
    if not atoms:
        print("FAIL  No atoms loaded (wrong chain filter?)", file=sys.stderr)
        sys.exit(1)

    pos = [a['pos'] for a in atoms]   # 0-based, matches C++ Structure::atoms

    # ── Read BAT ──────────────────────────────────────────────────────────
    try:
        anchor, entries, natoms = read_bat_file(bat_path)
    except Exception as exc:
        print(f"FAIL  BAT read: {exc}", file=sys.stderr)
        sys.exit(1)

    print(f"BAT entries       : {len(entries)}")
    aN  = anchor['N_idx']
    aCA = anchor['CA_idx']
    aC  = anchor['C_idx']
    print(f"Anchor indices    : N={aN}  CA={aCA}  C={aC}")

    # ── Sanity: atom count consistency ────────────────────────────────────
    if natoms is not None and natoms != len(atoms):
        print(f"FAIL  NATOMS in .bat ({natoms}) != PDB atom count ({len(atoms)})",
              file=sys.stderr)
        sys.exit(1)

    expected_entries = len(atoms) - 3
    if len(entries) != expected_entries:
        print(f"FAIL  entry count {len(entries)} != atoms-3 ({expected_entries})",
              file=sys.stderr)
        sys.exit(1)
    print(f"Count invariant   : {len(entries)} == {len(atoms)} - 3  PASS")

    # ── Sanity: anchor positions match PDB ────────────────────────────────
    for label, idx_key, pos_key in [('N',  'N_idx',  'N_pos'),
                                     ('CA', 'CA_idx', 'CA_pos'),
                                     ('C',  'C_idx',  'C_pos')]:
        idx     = anchor[idx_key]
        bat_pos = anchor[pos_key]
        pdb_pos = pos[idx]
        d = _dist(bat_pos, pdb_pos)
        if d > 1e-6:
            print(f"FAIL  anchor {label} pos mismatch: "
                  f"bat={bat_pos}  pdb={pdb_pos}  dist={d:.3e}",
                  file=sys.stderr)
            sys.exit(1)

    print("Anchor pos check  : PASS")

    # ── Per-entry geometry validation ─────────────────────────────────────
    bond_errs    = []
    angle_errs   = []
    torsion_errs = []
    failures     = []

    for e in entries:
        v       = e['atom_idx']
        p, g, h = resolve_refs(e['p'], e['g'], e['h'], aN, aC)

        pv = pos[v]
        pp = pos[p]
        pg = pos[g]
        ph = pos[h]

        # Independent computation from PDB coordinates
        bond_calc    = _dist(pv, pp)
        angle_calc   = _angle(pv, pp, pg)          # angle at pp, looking v–p–g
        torsion_calc = _dihedral(pv, pp, pg, ph)   # dihedral around p–g bond

        bond_err    = abs(bond_calc    - e['bond'])
        angle_err   = abs(angle_calc   - e['angle'])
        torsion_err = _torsion_diff(torsion_calc, e['torsion'])

        bond_errs.append(bond_err)
        angle_errs.append(angle_err)
        torsion_errs.append(torsion_err)

        if bond_err > TOLERANCE or angle_err > TOLERANCE or torsion_err > TOLERANCE:
            failures.append({
                'atom_idx':    v,
                'name':        atoms[v]['name'],
                'resid':       atoms[v]['resid'],
                'bond_calc':   bond_calc,
                'bond_bat':    e['bond'],
                'bond_err':    bond_err,
                'angle_calc':  angle_calc,
                'angle_bat':   e['angle'],
                'angle_err':   angle_err,
                'torsion_calc': torsion_calc,
                'torsion_bat': e['torsion'],
                'torsion_err': torsion_err,
            })

    # ── Report statistics ─────────────────────────────────────────────────
    def _stats(errs):
        return max(errs), sum(errs) / len(errs)

    bmax, bmean = _stats(bond_errs)
    amax, amean = _stats(angle_errs)
    tmax, tmean = _stats(torsion_errs)

    print(f"\nBond   error (Å)  : max={bmax:.3e}  mean={bmean:.3e}")
    print(f"Angle  error (rad): max={amax:.3e}  mean={amean:.3e}")
    print(f"Torsion error(rad): max={tmax:.3e}  mean={tmean:.3e}")

    if failures:
        print(f"\nFAIL  {len(failures)} / {len(entries)} entries exceed "
              f"tolerance {TOLERANCE:.0e}:", file=sys.stderr)
        for f in failures[:20]:
            print(f"  atom {f['atom_idx']:5d} {f['name']:<5s} res {f['resid']:4d} : "
                  f"bond_err={f['bond_err']:.2e}  "
                  f"angle_err={f['angle_err']:.2e}  "
                  f"torsion_err={f['torsion_err']:.2e}",
                  file=sys.stderr)
        if len(failures) > 20:
            print(f"  ... and {len(failures)-20} more.", file=sys.stderr)
        sys.exit(1)

    print(f"\nAll {len(entries)} entries match within tolerance {TOLERANCE:.0e} :  PASS")
    sys.exit(0)


if __name__ == '__main__':
    main()
