from typing import Tuple, Optional, cast

import MDAnalysis as mda 
from MDAnalysis.core.universe import Universe
from MDAnalysis.core.groups import AtomGroup
import numpy as np 
import argparse 


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
    parser.add_argument("-o","--out", type=str, required=True, help="Output BAM file path.")
    parser.add_argument("--sel", type=str, default="chainid PROA", help="MDAnalysis atom selection string (default: all).")
    parser.add_argument("--no-in-memory", dest="in_memory", action="store_false", help="Do not load trajectory into memory.")
    parser.add_argument("--bam-mode", choices=["full", "torsion"], default="full", help="BAM storage mode (default: full).")
    parser.set_defaults(in_memory=True)

    return parser.parse_args()

def load_universe(top:str, traj:str, selection:str, in_memory:bool) -> Tuple[Universe, AtomGroup]:
    """
    Load an MDAnalysis Universe and apply an atom selection.

    This function constructs an MDAnalysis Universe from the provided
    topology and trajectory files and returns both the Universe and an
    AtomGroup corresponding to the specified atom selection. The trajectory
    can either be streamed from disk or loaded entirely into memory,
    depending on the value of ``in_memory``.

    Parameters
    ----------
    top (str)  
        Path to the topology file (e.g., PDB, PSF, GRO, TPR).  
    traj (str)  
        Path to the trajectory file (e.g., XTC, DCD, TRR).  
    selection (str)  
        MDAnalysis atom selection string used to define the AtomGroup.  
    in_memory (bool)  
        Whether to load the entire trajectory into memory. If False, the trajectory is streamed from disk.

    Returns
    -------
    Universe  
        The loaded MDAnalysis Universe object.  
    AtomGroup  
        AtomGroup corresponding to the specified atom selection.  

    Raises
    ------
    ValueError  
        If the trajectory contains zero frames or if no atoms are found in the topology.
    """
    U = mda.Universe(top, traj, in_memory=in_memory)
    sel = U.select_atoms(selection)

    # immediate sanity check 
    if U.trajectory.n_frames == 0:
        raise ValueError("Trajectory contains zero frames.")

    assert U.atoms is not None # PyLance: narrow Optional[AtomGroup]
    if U.atoms.n_atoms == 0:  
        raise ValueError("No atoms found in topology.")
    
    # log the trajectory information: number of atoms and frames
    print(f"Number of atoms: {U.atoms.n_atoms}")
    print(f"Number of frames: {U.trajectory.n_frames}")

    return U, sel

def validate_torsions(ref: np.ndarray, custom: np.ndarray, atol: float = 1e-6) -> None:
    """
    Validate torsional angles against a reference, allowing for
    known convention differences.

    Accepted equivalences:
    1) ref ≈ custom
    2) |ref + custom| ≈ π   (supplementary convention)
    3) ref ≈ custom (mod 2π)

    Raises
    ------
    AssertionError
        If none of the equivalence checks pass.
    """
    # Case 1: direct equality
    try:
        np.testing.assert_allclose(ref, custom, atol=atol)
        return
    except AssertionError:
        pass

    # Case 2: supplementary (π shift)
    try:
        np.testing.assert_allclose(
            np.abs(ref + custom),
            np.pi,
            atol=atol
        )
        return
    except AssertionError:
        pass

    # Case 3: circular equivalence (mod 2π)
    diff = np.arctan2(
        np.sin(ref - custom),
        np.cos(ref - custom)
    )
    try:
        np.testing.assert_allclose(diff, 0.0, atol=atol)
        return
    except AssertionError:
        pass

    # If we reach here, none matched
    raise AssertionError(
        "Torsion validation failed: no known equivalence "
        "(direct, supplementary, or circular) matched."
    )