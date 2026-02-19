import numpy as np
import argparse
import util 
import cartesian_bat_converter as cbc
import tester

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
    
    dihedral_indices = cbc.generate_torsion_indices(sel)

    ref_torsions = tester.generate_torsion_ground_truth(U, sel, start=args.start, stop=args.stop, step=1)
    ref_torsions = cbc.ensure_radians(ref_torsions, name="BAT torsions")

    custom_torsions = []
    custom_torsions_numba = []
    for ts in U.trajectory[args.start:args.stop]:
        positions = sel.positions.copy()
        torsions = cbc.compute_torsions(positions, dihedral_indices)
        custom_torsions.append(torsions)
        torsions = cbc.compute_torsions_numba(positions, dihedral_indices)
        custom_torsions_numba.append(torsions)

    custom_torsions = np.array(custom_torsions)
    custom_torsions = cbc.ensure_radians(custom_torsions, name="custom torsions")
    custom_torsions_numba = np.array(custom_torsions_numba)
    custom_torsions_numba = cbc.ensure_radians(custom_torsions_numba, name="custom torsions numba")


    # Verify that the torsions match (up to periodicity)
    util.validate_torsions(ref_torsions, custom_torsions_numba)
    util.validate_torsions(custom_torsions, custom_torsions_numba)

if __name__=="__main__":
    main()