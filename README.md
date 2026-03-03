# Cartesian ↔ BAT Converter

Converts a protein PDB file to **Bond-Angle-Torsion (BAT)** internal coordinates, and reconstructs Cartesian coordinates from a BAT file. Written in pure C++17 with no external dependencies.

---

## Contents

- [Build](#build)
- [CLI Usage](#cli-usage)
- [API Reference](#api-reference)
  - [Data Types](#data-types)
  - [Geometry Functions](#geometry-functions)
  - [I/O](#io)
  - [Structure Analysis](#structure-analysis)
  - [BAT Coordinates](#bat-coordinates)
- [BAT File Format](#bat-file-format)
- [End-to-End Example](#end-to-end-example)
- [Error Handling](#error-handling)
- [Known Limitations](#known-limitations)

---

## Build

**Requirements:** C++17, CMake ≥ 3.15. No external libraries.

```bash
cd cpp
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

**Executables produced in `cpp/build/`:**

| Binary | Purpose |
|--------|---------|
| `bat_convert` | Main CLI driver |
| `test_bat_coords` | Unit test: BAT computation |
| `test_bat_roundtrip` | Unit test: compute → reconstruct round-trip |
| `test_bat_writer` | Unit test: file I/O round-trip |
| `test_pdb_roundtrip` | Integration test: full pipeline on a real PDB |
| `test_bfs_tree` | Unit test: BFS spanning forest |
| `test_backbone_frame` | Unit test: backbone frame selection |
| `test_bonds` | Unit test: bond inference |

**Run integration test:**
```bash
build/test_pdb_roundtrip path/to/protein.pdb [CHAIN]
```

---

## CLI Usage

```
bat_convert <input.pdb> <output.bat> [--chain <C>]
```

| Argument | Description |
|----------|-------------|
| `input.pdb` | Input PDB file (ATOM records only; HETATM ignored) |
| `output.bat` | Output BAT file (created or overwritten) |
| `--chain C` | Optional. Filter to single chain `C` (single character) |

**Exit codes:** `0` = success, `1` = error (message to stderr).

**Examples:**
```bash
# Single-chain protein
./build/bat_convert 1gia.pdb 1gia.bat

# Select chain A from multi-chain structure
./build/bat_convert 6crk.pdb 6crk_A.bat --chain A
```

---

## API Reference

All symbols are in `namespace geometry`. Include headers from `cpp/include/`.

Link against the `cbc` static library (built by CMake).

---

### Data Types

#### `Vec3` — `cartesian/vec3.hpp`

A 3D point/vector using `double` components.

```cpp
struct Vec3 {
    double x, y, z;

    Vec3();                             // default-constructible
    Vec3(double x_, double y_, double z_);
};
```

Coordinates are in **Ångströms**.

---

#### `Atom` — `structure/structure.hpp`

Represents a single atom parsed from a PDB file.

```cpp
struct Atom {
    int         index;         // PDB serial number (1-based)
    std::string name;          // Atom name, e.g. "CA", "CB"
    std::string residue_name;  // Residue name, e.g. "ALA"
    int         residue_id;    // Residue sequence number
    char        chain_id;      // Chain identifier
    std::string element;       // Element symbol: "H", "C", "N", "O", "S"
    Vec3        position;      // Cartesian coordinates (Å)
    char        alt_loc = ' '; // Alternate location indicator
    double      occupancy = 1.0;
};
```

---

#### `Bond` — `structure/structure.hpp`

An undirected covalent bond between two atoms.

```cpp
struct Bond {
    std::size_t i, j;  // Indices into Structure::atoms (unordered)
};
```

---

#### `Structure` — `structure/structure.hpp`

The complete molecular structure.

```cpp
struct Structure {
    std::vector<Atom> atoms;
    std::vector<Bond> bonds;

    // Returns a new Structure containing only atoms from the given chain.
    // If chain == '\0', returns a copy of the whole structure.
    Structure select_chain(char chain) const;
};
```

`bonds` is populated by `infer_bonds()`. Before that call, it is empty.

---

#### `ProteinFrame` — `structure/protein_root.hpp`

Backbone anchor indices for the first N/CA/C residue found.

```cpp
struct ProteinFrame {
    std::size_t N  = npos;   // Index of backbone N
    std::size_t CA = npos;   // Index of backbone CA
    std::size_t C  = npos;   // Index of backbone C

    bool valid() const;      // true iff all three indices are set
};
```

Indices refer to positions in `Structure::atoms`.

---

#### `BFSTree` — `structure/bfs_tree.hpp`

BFS spanning forest over the molecular graph.

```cpp
struct BFSTree {
    static constexpr std::size_t npos = static_cast<std::size_t>(-1);

    std::size_t n = 0;                               // Number of nodes

    std::vector<std::size_t> parent;                 // parent[v], or npos if root
    std::vector<int>         depth;                  // depth[v] from its root
    std::vector<std::size_t> root_of;                // component root of v
    std::vector<std::size_t> order;                  // BFS discovery order (whole forest)
    std::vector<std::vector<std::size_t>> children;  // tree edges: children[u]

    struct Ref3 { std::size_t p, g, h; };  // parent / grandparent / great-grandparent
    std::vector<Ref3> refs;                // refs[v]; entries are npos near roots
};
```

`refs[v]` stores the three BAT reference atoms for atom `v`:
- `p` = parent in spanning tree
- `g` = grandparent (`parent[p]`), or `npos` if depth < 2
- `h` = great-grandparent, or `npos` if depth < 3

---

#### `BATCoords` — `bat/bat_coords.hpp`

BAT representation of an entire protein structure.

```cpp
struct BATCoords {

    struct Anchor {
        std::size_t N_idx  = npos;   // Index of backbone N in Structure::atoms
        std::size_t CA_idx = npos;   // Index of backbone CA
        std::size_t C_idx  = npos;   // Index of backbone C
        Vec3 N_pos{};                // Absolute Cartesian position of N (Å)
        Vec3 CA_pos{};               // Absolute Cartesian position of CA (Å)
        Vec3 C_pos{};                // Absolute Cartesian position of C (Å)
    };
    Anchor anchor{};

    struct Entry {
        std::size_t atom_idx;  // Index into Structure::atoms
        double bond;           // Distance to BFS parent (Å)
        double angle;          // Bond angle at parent (radians, [0, π])
        double torsion;        // Dihedral angle (radians, (−π, π])
    };
    std::vector<Entry> entries;  // One entry per non-anchor atom, in BFS order
};
```

The three anchor atoms carry absolute Cartesian positions that define the global frame. All other atoms are described relative to their BFS parent/grandparent/great-grandparent.

---

#### `BATFile` — `io/bat_writer.hpp`

Return type of `read_bat()`. Holds everything needed for `reconstruct_cartesian()`.

```cpp
struct BATFile {
    BATCoords bat;   // Anchor positions + per-atom (bond, angle, torsion)
    BFSTree   tree;  // Only n and refs[] are populated
};
```

---

### Geometry Functions

#### `geometry.hpp` — `cartesian/geometry.hpp`

Low-level geometric primitives. All angles are in **radians**.

```cpp
// Euclidean distance between two points (Å).
double distance(const Vec3& a, const Vec3& b);

// Bond angle at vertex b (radians, [0, π]).
double angle(const Vec3& a, const Vec3& b, const Vec3& c);

// Dihedral angle around the b–c bond (radians, (−π, π]).
// Convention: b1 = b−a, b2 = c−b, b3 = d−c
//   n1 = cross(b1,b2), n2 = cross(b2,b3)
//   returns atan2( dot(cross(n1,n2), normalize(b2)), dot(n1,n2) )
double dihedral(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d);
```

---

### I/O

#### `PDBReader::read` — `io/pdb_reader.hpp`

Parses a PDB file (ATOM records only) into a `Structure`.

```cpp
static Structure PDBReader::read(const std::string& filename, char chain = '\0');
```

| Parameter | Description |
|-----------|-------------|
| `filename` | Path to PDB file |
| `chain` | If `'\0'`, reads all chains; otherwise reads only the named chain |

**Alternate locations:** the highest-occupancy conformation is kept; ties broken by preferring altLoc `' '` then `'A'`. If two conformations have equal occupancy and different positions, an error is thrown.

**Supported elements:** H, C, N, O, S (protein-only).

**Throws** `std::runtime_error` on: file not found, malformed coordinates, unsupported element, or ambiguous alternate locations.

---

#### `write_bat` — `io/bat_writer.hpp`

Writes BAT coordinates to a text file.

```cpp
void write_bat(const std::string& filename,
               const BATCoords&   bat,
               const BFSTree&     tree);
```

The file stores anchor Cartesian positions and, per non-anchor atom, the BFS reference indices and (bond, angle, torsion) triple. Doubles are written at full IEEE 754 precision (17 significant digits).

**Throws** `std::runtime_error` on I/O error.

---

#### `read_bat` — `io/bat_writer.hpp`

Reads a BAT file produced by `write_bat()`.

```cpp
BATFile read_bat(const std::string& filename);
```

Returns a `BATFile` with `bat` and `tree` populated. The `tree` only has `n` and `refs[]` set — sufficient for `reconstruct_cartesian()`. The original `Structure` is not required.

**Throws** `std::runtime_error` on: file not found, malformed input, or out-of-range indices.

---

### Structure Analysis

#### `infer_bonds` — `structure/bond_inference.hpp`

Detects covalent bonds from atomic positions using a spatial hash and covalent radii.

```cpp
struct BondInferOptions {
    double tolerance_angstrom      = 0.35;  // tolerance added to sum of covalent radii
    double min_distance            = 0.60;  // distance below this → error (atom overlap)
    double max_bond_search_distance = 2.5;  // distance above this → skip (neighbour cutoff)
    bool   include_hydrogen        = true;  // if false, H atoms are treated as non-bonding
    bool   enforce_valence         = true;  // hard error if max valence is exceeded
};

std::vector<Bond> infer_bonds(const Structure& st, const BondInferOptions& opt);
```

A bond is accepted when:
```
distance(i, j) ≤ covalent_radius(i) + covalent_radius(j) + tolerance_angstrom
```

**Covalent radii used (Å):** H 0.31, C 0.76, N 0.71, O 0.66, S 1.05.

**Max valence:** H 1, C 4, N 4, O 2, S 2.

**Throws** `std::runtime_error` on atom overlap, valence exceeded, or unsupported element.

---

#### `build_adjacency` — `structure/structure_graph.hpp`

Builds an undirected adjacency list from `Structure::bonds`.

```cpp
std::vector<std::vector<std::size_t>> build_adjacency(const Structure& st);
```

`adj[i]` contains the 0-based indices of all atoms bonded to atom `i`. Requires `st.bonds` to be populated (call `infer_bonds` first).

---

#### `connected_component_labels` — `structure/structure_graph.hpp`

Returns a component label for each atom.

```cpp
std::vector<int> connected_component_labels(const Structure& st);
```

`labels[i]` is the 0-based component index of atom `i`. Atoms in the same connected component share the same label. Requires `st.bonds` to be populated.

---

#### `select_first_backbone_frame` — `structure/protein_root.hpp`

Finds the first residue in the structure with all three backbone atoms (N, CA, C).

```cpp
ProteinFrame select_first_backbone_frame(
    const Structure&    st,
    std::optional<char> chain = std::nullopt
);
```

If `chain` is `std::nullopt`, the first chain in lexicographic order that yields a valid frame is used. Within a chain, residues are searched in increasing `residue_id` order.

**Throws** `std::runtime_error` if no valid backbone frame is found.

---

#### `build_bfs_forest` — `structure/bfs_tree.hpp`

Builds a deterministic BFS spanning forest over the molecular graph.

```cpp
BFSTree build_bfs_forest(
    const std::vector<std::vector<std::size_t>>& adj,
    std::size_t first_root,
    bool sort_neighbors = true
);
```

BFS starts at `first_root`. Any nodes not reached from `first_root` (disconnected components) become additional BFS roots in order of their index. Setting `sort_neighbors = true` ensures deterministic output across runs.

**Convenience alias** (same behaviour, semantically signals a connected graph):
```cpp
inline BFSTree build_bfs_tree(
    const std::vector<std::vector<std::size_t>>& adj,
    std::size_t root,
    bool sort_neighbors = true
);
```

---

### BAT Coordinates

#### `compute_bat` — `bat/bat_coords.hpp`

Computes BAT coordinates from a structure and its BFS spanning tree.

```cpp
BATCoords compute_bat(
    const Structure&    st,
    const BFSTree&      tree,
    const ProteinFrame& frame
);
```

The three anchor atoms (N, CA, C) are stored as absolute Cartesian positions. Every other atom gets a `(bond, angle, torsion)` triple relative to its BFS parent, grandparent, and great-grandparent.

**Virtual-reference convention** (applied when BFS refs are `npos`, i.e., for atoms near the BFS root):
- Missing grandparent → use `anchor.N_idx`
- Missing great-grandparent → use `anchor.C_idx`, unless `anchor.C_idx == parent`, in which case use `anchor.N_idx`

**Throws** `std::runtime_error` if: frame is invalid, tree/structure sizes mismatch, or a non-anchor disconnected-component root is encountered.

---

#### `place_atom` — `bat/bat_reconstruct.hpp`

Places one atom using the **Natural Extension Reference Frame (NeRF)** formula.

```cpp
Vec3 place_atom(const Vec3& h, const Vec3& g, const Vec3& p,
                double bond, double angle, double torsion);
```

Given great-grandparent `h`, grandparent `g`, parent `p`, and BAT parameters, returns position `v` satisfying:
- `distance(v, p) = bond`
- `angle(v, p, g)` at vertex `p` = `angle`
- `dihedral(v, p, g, h) = torsion`

**Algorithm:**
1. Build a right-handed local frame at `p`:
   - `e_x = normalize(p − g)`
   - `e_y = normalize((h − g) − proj of (h−g) onto e_x)` (fallback to arbitrary perp if collinear)
   - `e_z = cross(e_x, e_y)`
2. Compute displacement: `dx = −bond·cos(angle)`, `dy = bond·sin(angle)·cos(torsion)`, `dz = bond·sin(angle)·sin(torsion)`
3. Return `p + dx·e_x + dy·e_y + dz·e_z`

---

#### `reconstruct_cartesian` — `bat/bat_reconstruct.hpp`

Reconstructs all atomic Cartesian positions from BAT coordinates.

```cpp
std::vector<Vec3> reconstruct_cartesian(
    const BATCoords& bat,
    const BFSTree&   tree
);
```

Returns a vector indexed identically to `Structure::atoms` (i.e., `positions[i]` is the position of atom `i`). Applies the same virtual-reference convention as `compute_bat`, guaranteeing exact round-trip recovery (error < 10⁻¹² Å on real PDBs).

**Throws** `std::runtime_error` if anchor indices are out of range.

---

## BAT File Format

Text format, version `BAT 1`.

```
BAT 1
NATOMS <n>
ANCHOR_N  <N_idx>   <x>  <y>  <z>
ANCHOR_CA <CA_idx>  <x>  <y>  <z>
ANCHOR_C  <C_idx>   <x>  <y>  <z>
ENTRY <atom_idx>  <p>  <g>  <h>  <bond>  <angle>  <torsion>
ENTRY <atom_idx>  <p>  <g>  <h>  <bond>  <angle>  <torsion>
...
```

| Field | Description |
|-------|-------------|
| `NATOMS` | Total number of atoms in the original structure |
| `ANCHOR_*` | 0-based atom index and absolute Cartesian coordinates (Å) |
| `atom_idx` | 0-based index of this atom |
| `p`, `g`, `h` | 0-based BFS parent / grandparent / great-grandparent indices; `-1` = `npos` |
| `bond` | Distance to parent (Å) |
| `angle` | Bond angle at parent (radians) |
| `torsion` | Dihedral angle (radians) |

`ENTRY` lines are written in BFS traversal order. Doubles are at full IEEE 754 precision (17 significant digits).

---

## End-to-End Example

```cpp
#include "io/pdb_reader.hpp"
#include "io/bat_writer.hpp"
#include "structure/bond_inference.hpp"
#include "structure/structure_graph.hpp"
#include "structure/protein_root.hpp"
#include "structure/bfs_tree.hpp"
#include "bat/bat_coords.hpp"
#include "bat/bat_reconstruct.hpp"

// --- Forward conversion: PDB → BAT ---

geometry::Structure st = geometry::PDBReader::read("protein.pdb", 'A');

geometry::BondInferOptions opts;
st.bonds = geometry::infer_bonds(st, opts);

auto adj   = geometry::build_adjacency(st);
auto frame = geometry::select_first_backbone_frame(st, 'A');
auto tree  = geometry::build_bfs_forest(adj, frame.CA, /*sort_neighbors=*/true);
auto bat   = geometry::compute_bat(st, tree, frame);

geometry::write_bat("output.bat", bat, tree);

// --- Reverse conversion: BAT → Cartesian ---

geometry::BATFile loaded  = geometry::read_bat("output.bat");
auto positions = geometry::reconstruct_cartesian(loaded.bat, loaded.tree);
// positions[i] is the Cartesian position of atom i (Å)
```

---

## Error Handling

All errors are thrown as `std::runtime_error` with a descriptive message. There are no silent failures. Common causes:

| Situation | Thrown by |
|-----------|-----------|
| PDB file not found or malformed | `PDBReader::read` |
| Unsupported element (not H/C/N/O/S) | `PDBReader::read`, `infer_bonds` |
| Atom overlap (distance < 0.6 Å) | `infer_bonds` |
| Valence exceeded | `infer_bonds` |
| No backbone N/CA/C found | `select_first_backbone_frame` |
| Disconnected chain (gap in structure) | `compute_bat` |
| BAT file not found or malformed | `read_bat` |

---

## Known Limitations

- **Protein-only:** Elements H, C, N, O, S are supported. Ligands, metals, and non-standard residues (HETATM) are ignored.
- **Single structure per file:** No trajectory support.
- **Disconnected chains:** A chain with missing residues (gap in electron density) produces a disconnected molecular graph and is rejected by `compute_bat`. This is by design.
- **Binary format:** Not yet implemented. Only text BAT files are currently supported.
