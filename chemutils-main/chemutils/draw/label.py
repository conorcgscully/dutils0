from collections.abc import Callable
from typing import Any

from openeye import oechem

LabelFunction = Callable[[oechem.OEAtomBase], Any]


def get_enhanced_stereo_label(atom: oechem.OEAtomBase) -> str:
    """
    Get enhanced stereochemistry label for an atom.

    Absolute stereo is indicated as a `abs`, 'and' stereo groups as `&` and 'or'
    stereo groups as `or`.

    Args:
        atom: OpenEye atom.

    Returns:
        Stereochemistry label.
    """
    labels = []
    index_and_group = 1
    index_or_group = 1
    for group in list(atom.GetParent().GetGroups(oechem.OEIsMDLStereoGroup())):
        match group.GetGroupType():
            case oechem.OEGroupType_MDLAbsStereo:
                if group.HasAtom(atom):
                    labels.append("abs")
            case oechem.OEGroupType_MDLAndStereo:
                if group.HasAtom(atom):
                    labels.append(f"&{index_and_group}")
                index_and_group += 1
            case oechem.OEGroupType_MDLOrStereo:
                if group.HasAtom(atom):
                    labels.append(f"or{index_or_group}")
                index_or_group += 1
    return ",".join(labels)
