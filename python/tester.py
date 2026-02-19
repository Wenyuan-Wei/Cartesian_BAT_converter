from typing import Optional
from MDAnalysis.core.universe import Universe
from MDAnalysis.core.groups import AtomGroup
from MDAnalysis.analysis.bat import BAT
from MDAnalysis.analysis.dihedrals import Dihedral
import numpy as np 
import cartesian_bat_converter as cbc 

def compute_bat_ground_truth(U: Universe, sel: AtomGroup, start: int = 0, stop: Optional[int] = None, step: int = 1) -> BAT:
    """
    Compute BAM (bond-angle-torsion) coordinates using MDAnalysis BAT.

    This function uses MDAnalysis's BAT implementation as a ground-truth
    reference for internal coordinate extraction. It is intended for
    validation and testing purposes, not for performance-critical use.

    Parameters
    ----------
    U : Universe
        MDAnalysis Universe containing the trajectory.
    sel : AtomGroup
        AtomGroup defining the atoms used for BAT conversion.
    start : int, optional
        Starting frame index (default: 0).
    stop : int or None, optional
        Stopping frame index (default: None, meaning last frame).
    step : int, optional
        Step size for frame iteration (default: 1).

    Returns
    -------
    np.ndarray
        Array of BAT coordinates with shape (n_frames, n_internal_coords).

    Raises
    ------
    ValueError  
        If the atom selection results in zero atoms.
    ValueError  
        If the atom selection contains multiple molecules.
    """

    stop_frame: int = (
        U.trajectory.n_frames if stop is None else stop
    )

    # quick sanity check for selection validty
    if sel.n_atoms == 0:
        raise ValueError("Atom selection contains zero atoms.")
    
    # check if the selection has multiple segments
    if len(sel.fragments) != 1:
        raise ValueError(
            f"BAT requires a single molecule, but selection contains "
            f"{len(sel.fragments)} molecules."
    )

    # MDAnalysis BAT expects a static AtomGroup
    bat = BAT(sel)

    bat.run(start=start, stop=stop_frame, step=step)

    return bat 

def generate_torsion_ground_truth(U: Universe, sel: AtomGroup, start: int = 0, stop: Optional[int] = None, step: int = 1) -> np.ndarray:
    """
    Compute torsional angles using MDAnalysis Dihedral analysis.

    This function uses MDAnalysis's Dihedral implementation as a ground-truth
    reference for torsional angle extraction. It is intended for
    validation and testing purposes, not for performance-critical use.

    Parameters
    ----------
    U : Universe
        MDAnalysis Universe containing the trajectory.
    sel : AtomGroup
        AtomGroup defining the atoms used for torsion calculation.
    start : int, optional
        Starting frame index (default: 0).
    stop : int or None, optional
        Stopping frame index (default: None, meaning last frame).
    step : int, optional
        Step size for frame iteration (default: 1).

    Returns
    -------
    np.ndarray
        Array of torsional angles with shape (n_frames, n_torsions).

    Raises
    ------
    ValueError  
        If the atom selection results in zero atoms.
    ValueError  
        If the atom selection contains multiple molecules.
    """

    assert U.atoms is not None  # PyLance: narrow Optional[AtomGroup]

    stop_frame: int = (
        U.trajectory.n_frames if stop is None else stop
    )

    # quick sanity check for selection validty
    if sel.n_atoms == 0:
        raise ValueError("Atom selection contains zero atoms.")
    
    # check if the selection has multiple segments
    if len(sel.fragments) != 1:
        raise ValueError(
            f"Dihedral analysis requires a single molecule, but selection contains "
            f"{len(sel.fragments)} molecules."
    )

    # Generate torsion indices
    dihedral_indices = cbc.generate_torsion_indices(sel)
    sel_to_atom = sel.atoms.indices
    dih_groups = [
        U.atoms[sel_to_atom[idx_row]]
        for idx_row in dihedral_indices
    ]

    dih_analysis = Dihedral(dih_groups)
    dih_analysis.run(start=start, stop=stop_frame, step=step)

    return dih_analysis.results.angles

