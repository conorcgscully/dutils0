from collections.abc import Generator

from openeye import oechem


def split_molecule(oemol: oechem.OEMolBase) -> Generator[oechem.OEMolBase, None, None]:
    """
    Split an OpenEye molecule into it's non-covalently bound parts.

    Args:
        oemol: OpenEye molecule to consider.

    Yields:
        Constituent OpenEye molecules.
    """
    num_parts, part_list = oechem.OEDetermineComponents(oemol)
    pred = oechem.OEPartPredAtom(part_list)

    for i in range(1, num_parts + 1):
        pred.SelectPart(i)
        subset_oemol = oechem.OEGraphMol()
        oechem.OESubsetMol(subset_oemol, oemol, pred)
        yield subset_oemol
