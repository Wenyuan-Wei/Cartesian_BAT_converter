from MDAnalysis.lib.distances import calc_bonds
import numpy as np
import argparse
import util 
import cartesian_bat_converter as cbc

def parse_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments for BAM conversion.

    Defines all CLI options required for reading topology and trajectory
    files, atom selection, memory handling, and BAM output mode.

    Returns
    -------
    argparse.Namespace
        Parsed command-line arguments with the following attributes:

        top (str)  
            Path to the topology file.  
        traj (str)  
            Path to the trajectory file.  
        out (str)  
            Output BAM file path.  
        sel (str)  
            MDAnalysis atom selection string.  
        in_memory (bool)  
            Whether to load the entire trajectory into memory.  
        bam_mode ({"full", "torsion"})  
            BAM storage mode.
    """
    parser = argparse.ArgumentParser(description="Convert Cartesian coordinates to BAM format.")
    parser.add_argument("-s","--top", type=str, required=True, help="Path to the topology file.")
    parser.add_argument("-f","--traj", type=str, required=True, help="Path to the trajectory file.")
    parser.add_argument("--start", type=int, default=0, help="Starting frame index (default: 0).")
    parser.add_argument("--stop", type=int, default=None, help="Stopping frame index (default: None, meaning last frame).")
    parser.add_argument("--sel", type=str, default="chainid PROA", help="MDAnalysis atom selection string (default: all).")
    parser.add_argument("--no-in-memory", dest="in_memory", action="store_false", help="Do not load trajectory into memory.")
    parser.set_defaults(in_memory=True)

    return parser.parse_args()

def main():
    args = parse_arguments()

    U, sel = util.load_universe(args.top, args.traj, args.sel, args.in_memory)
    
    assert U.atoms is not None  # PyLance: narrow Optional[AtomGroup]

    bond_idx = cbc.generate_bond_indices(sel)
    i=bond_idx[:, 0]
    j=bond_idx[:, 1]

    ref_bonds=[]
    custom_bonds=[]
    custom_bonds_numba=[]

    for ts in U.trajectory[args.start:args.stop]:

        xyz = sel.positions.astype(np.float64, copy=False)
        
        r_i=xyz[i]
        r_j=xyz[j]
    
        ref_bond = calc_bonds(r_i, r_j)
        custom_bond = cbc.compute_bonds(xyz, bond_idx)
        custom_bond_numba = cbc.compute_bonds_numba(xyz, bond_idx)

        ref_bonds.append(ref_bond)
        custom_bonds.append(custom_bond)
        custom_bonds_numba.append(custom_bond_numba)

    ref_bonds = np.asarray(ref_bonds)
    custom_bonds = np.asarray(custom_bonds)
    custom_bonds_numba = np.asarray(custom_bonds_numba)

    np.testing.assert_allclose(
        ref_bonds,
        custom_bonds,
        atol=1e-6
    )

    np.testing.assert_allclose(
        custom_bonds,
        custom_bonds_numba,
        atol=1e-6
    )

if __name__=="__main__":
    main()