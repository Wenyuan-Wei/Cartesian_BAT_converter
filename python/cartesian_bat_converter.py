from typing import cast
from MDAnalysis.core.groups import AtomGroup
import numpy as np
from numba import njit


def ensure_radians(angles: np.ndarray, *, name: str = "angles") -> np.ndarray:
    """
    Ensure angle array is in radians.

    Detects whether angles are in degrees or radians based on magnitude
    and converts to radians if necessary.

    Parameters
    ----------
    angles : np.ndarray
        Angle values from an external source.
    name : str
        Name used in error messages.

    Returns
    -------
    np.ndarray
        Angles in radians.

    Raises
    ------
    ValueError
        If units cannot be determined unambiguously.
    """
    max_abs = np.nanmax(np.abs(angles))

    if max_abs <= np.pi + 1e-6:
        # Already radians
        return angles

    if max_abs <= 180.0 + 1e-6:
        # Degrees → radians
        return np.deg2rad(angles)

    raise ValueError(
        f"Cannot determine units for {name}: max |value| = {max_abs}. "
        "Expected radians (≤ π) or degrees (≤ 180)."
    )

def generate_torsion_indices(sel: AtomGroup) -> np.ndarray:
    """
    Generate torsion index array for a given AtomGroup.

    Only dihedrals fully contained within the selection are kept.

    Returns
    -------
    np.ndarray, shape (n_torsions, 4)
        Torsion indices relative to sel.positions.
    """
    # Universe atom indices present in selection
    sel_atom_indices = sel.atoms.indices
    sel_atom_set = set(sel_atom_indices)

    # Map Universe atom index -> selection-local index
    atom_to_sel = {
        atom_idx: i
        for i, atom_idx in enumerate(sel_atom_indices)
    }

    # Dihedral definitions in Universe indices
    dih_idx = cast(np.ndarray, sel.dihedrals.indices)

    torsion_idx = []

    for a, b, c, d in dih_idx:
        # keep only dihedrals fully inside the selection
        if (
            a in sel_atom_set
            and b in sel_atom_set
            and c in sel_atom_set
            and d in sel_atom_set
        ):
            torsion_idx.append([
                atom_to_sel[a],
                atom_to_sel[b],
                atom_to_sel[c],
                atom_to_sel[d],
            ])

    return np.asarray(torsion_idx, dtype=np.int32)

def compute_torsions(xyz: np.ndarray, torsion_idx: np.ndarray) -> np.ndarray:
    """
    Compute torsional (dihedral) angles for a set of atom quadruplets.

    Parameters
    ----------
    xyz : np.ndarray, shape (n_atoms, 3)
        Cartesian coordinates.
    torsion_idx : np.ndarray, shape (n_torsions, 4)
        Atom indices (i, j, k, l) defining each torsion.

    Returns
    -------
    np.ndarray, shape (n_torsions,)
        Torsional angles in radians, in the range (-pi, pi].
    """
    i = torsion_idx[:, 0]
    j = torsion_idx[:, 1]
    k = torsion_idx[:, 2]
    l = torsion_idx[:, 3]

    r_i = xyz[i]
    r_j = xyz[j]
    r_k = xyz[k]
    r_l = xyz[l]

    b1 = r_i - r_j
    b2 = r_k - r_j
    b3 = r_l - r_k

    # Normalize b2 for stability
    b2_norm = b2 / np.linalg.norm(b2, axis=1)[:, None]

    # Compute normals
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)

    n1 /= np.linalg.norm(n1, axis=1)[:, None]
    n2 /= np.linalg.norm(n2, axis=1)[:, None]

    # Orthogonal vector
    m1 = np.cross(n1, b2_norm)

    x = np.sum(n1 * n2, axis=1)
    y = np.sum(m1 * n2, axis=1)

    return np.arctan2(y, x)

@njit(fastmath=True, cache=True)
def compute_torsions_numba(xyz: np.ndarray, torsion_idx: np.ndarray) -> np.ndarray:
    n_torsions = torsion_idx.shape[0]
    out = np.empty(n_torsions, dtype=np.float64)

    for t in range(n_torsions):
        i = torsion_idx[t, 0]
        j = torsion_idx[t, 1]
        k = torsion_idx[t, 2]
        l = torsion_idx[t, 3]

        # coordinates
        ri = xyz[i]
        rj = xyz[j]
        rk = xyz[k]
        rl = xyz[l]

        # bond vectors
        b1 = ri - rj
        b2 = rk - rj
        b3 = rl - rk

        # normalize b2
        b2_norm = b2 / np.sqrt(b2[0]**2 + b2[1]**2 + b2[2]**2)

        # plane normals
        n1 = np.cross(b1, b2)
        n2 = np.cross(b2, b3)

        n1 /= np.sqrt(n1[0]**2 + n1[1]**2 + n1[2]**2)
        n2 /= np.sqrt(n2[0]**2 + n2[1]**2 + n2[2]**2)

        # orthogonal vector
        m1 = np.cross(n1, b2_norm)

        x = n1[0]*n2[0] + n1[1]*n2[1] + n1[2]*n2[2]
        y = m1[0]*n2[0] + m1[1]*n2[1] + m1[2]*n2[2]

        out[t] = np.arctan2(y, x)

    return out

def generate_angle_indices(sel: AtomGroup) -> np.ndarray:
    """
    Generate bond-angle index array for a given AtomGroup.

    Only angles fully contained within the selection are kept.

    Returns
    -------
    np.ndarray, shape (n_angles, 3)
        Angle indices relative to sel.positions.
    """
    sel_atom_indices = sel.atoms.indices
    sel_atom_set = set(sel_atom_indices)

    atom_to_sel = {
        atom_idx: i
        for i, atom_idx in enumerate(sel_atom_indices)
    }

    ang_idx = cast(np.ndarray, sel.angles.indices)

    angle_idx = []

    for i, j, k in ang_idx:
        if (
            i in sel_atom_set
            and j in sel_atom_set
            and k in sel_atom_set
        ):
            angle_idx.append([
                atom_to_sel[i],
                atom_to_sel[j],
                atom_to_sel[k],
            ])

    return np.asarray(angle_idx, dtype=np.int32)

def compute_angles(xyz: np.ndarray, angle_idx: np.ndarray) -> np.ndarray:
    """
    Compute bond angles for a set of atom triplets.

    Parameters
    ----------
    xyz : np.ndarray, shape (n_atoms, 3)
        Cartesian coordinates.
    angle_idx : np.ndarray, shape (n_angles, 3)
        Atom indices (i, j, k) defining each angle.

    Returns
    -------
    np.ndarray, shape (n_angles,)
        Bond angles in radians, in [0, pi].
    """
    i = angle_idx[:, 0]
    j = angle_idx[:, 1]
    k = angle_idx[:, 2]

    ri = xyz[i]
    rj = xyz[j]
    rk = xyz[k]

    v1 = ri - rj
    v2 = rk - rj

    dot = np.sum(v1 * v2, axis=1)
    norm1 = np.linalg.norm(v1, axis=1)
    norm2 = np.linalg.norm(v2, axis=1)

    cos_theta = dot / (norm1 * norm2)

    # Numerical safety
    cos_theta = np.clip(cos_theta, -1.0, 1.0)

    return np.arccos(cos_theta)

@njit(fastmath=True, cache=True)
def compute_angles_numba(xyz: np.ndarray, angle_idx: np.ndarray) -> np.ndarray:
    n_angles = angle_idx.shape[0]
    out = np.empty(n_angles, dtype=np.float64)

    for t in range(n_angles):
        i = angle_idx[t, 0]
        j = angle_idx[t, 1]
        k = angle_idx[t, 2]

        xi, yi, zi = xyz[i]
        xj, yj, zj = xyz[j]
        xk, yk, zk = xyz[k]

        v1x = xi - xj
        v1y = yi - yj
        v1z = zi - zj

        v2x = xk - xj
        v2y = yk - yj
        v2z = zk - zj

        dot = v1x*v2x + v1y*v2y + v1z*v2z

        n1 = np.sqrt(v1x*v1x + v1y*v1y + v1z*v1z)
        n2 = np.sqrt(v2x*v2x + v2y*v2y + v2z*v2z)

        cos_theta = dot / (n1 * n2)

        if cos_theta > 1.0:
            cos_theta = 1.0
        elif cos_theta < -1.0:
            cos_theta = -1.0

        out[t] = np.arccos(cos_theta)

    return out

def generate_bond_indices(sel: AtomGroup) -> np.ndarray:
    """
    Generate bond index array for a given AtomGroup.

    Only bonds fully contained within the selection are kept.

    Returns
    -------
    np.ndarray, shape (n_bonds, 2)
        Bond indices relative to sel.positions.
    """
    sel_atom_indices = sel.atoms.indices
    sel_atom_set = set(sel_atom_indices)

    atom_to_sel = {
        atom_idx: i
        for i, atom_idx in enumerate(sel_atom_indices)
    }

    bond_idx = cast(np.ndarray, sel.bonds.indices)

    bonds = []

    for i, j in bond_idx:
        if i in sel_atom_set and j in sel_atom_set:
            bonds.append([
                atom_to_sel[i],
                atom_to_sel[j],
            ])

    return np.asarray(bonds, dtype=np.int32)

def compute_bonds(xyz: np.ndarray, bond_idx: np.ndarray) -> np.ndarray:
    """
    Compute bond lengths for a set of atom pairs.

    Parameters
    ----------
    xyz : np.ndarray, shape (n_atoms, 3)
        Cartesian coordinates.
    bond_idx : np.ndarray, shape (n_bonds, 2)
        Atom indices (i, j) defining each bond.

    Returns
    -------
    np.ndarray, shape (n_bonds,)
        Bond lengths (same length units as xyz).
    """
    i = bond_idx[:, 0]
    j = bond_idx[:, 1]

    diff = xyz[i] - xyz[j]
    return np.linalg.norm(diff, axis=1)

@njit(fastmath=True, cache=True)
def compute_bonds_numba(xyz: np.ndarray, bond_idx: np.ndarray) -> np.ndarray:
    n_bonds = bond_idx.shape[0]
    out = np.empty(n_bonds, dtype=np.float64)

    for t in range(n_bonds):
        i = bond_idx[t, 0]
        j = bond_idx[t, 1]

        dx = xyz[i, 0] - xyz[j, 0]
        dy = xyz[i, 1] - xyz[j, 1]
        dz = xyz[i, 2] - xyz[j, 2]

        out[t] = np.sqrt(dx*dx + dy*dy + dz*dz)

    return out

