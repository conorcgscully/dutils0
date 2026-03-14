import numpy as np
import numpy.typing as npt
from openeye import oechem


class InvalidSMARTSError(ValueError):
    """Error thrown when an invalid SMARTS string is parsed."""

    pass


def get_atom_mapping_from_match(
    match: oechem.OEMatchBase, /
) -> dict[oechem.OEAtomBase, oechem.OEAtomBase]:
    """Convert a match into a dictionary mapping source atoms to target atoms."""
    return dict(zip(match.GetPatternAtoms(), match.GetTargetAtoms(), strict=True))


def get_target_mask_from_match(match: oechem.OEMatchBase, /) -> npt.NDArray[np.bool_]:
    """Get a boolean mask over the target atoms in an OpenEye match, indicating which atoms are matched."""
    target = next(iter(match.GetTargetAtoms())).GetParent()
    n_atoms = target.NumAtoms()
    mask = np.zeros(n_atoms, dtype=bool)
    idxs = [atom.GetIdx() for atom in match.GetTargetAtoms()]
    mask[idxs] = True
    return mask
