from openeye import oechem

OUTER_ELECTRON_COUNT = {7: 5, 8: 6, 9: 7, 16: 6}
"""Mapping of atomic number to outer electron count."""

HYBRIDISATION_TO_SP_ORBITALS: dict[int, int] = {
    oechem.OEHybridization_sp: 2,
    oechem.OEHybridization_sp2: 3,
    oechem.OEHybridization_sp3: 4,
}
"""Mapping of OpenEye hybridisation to number of sp orbitals."""


def get_hybridization(atom: oechem.OEAtomBase, /) -> int:
    """
    Calculate the hybridization of an atom.

    This uses OpenEye's `OEGetHybridization`, with some special cases to catch examples
    where the hybridization is Unknown.

    Args:
        atom: OpenEye atom.

    Returns:
        Hybridization, as a value in the `OEHybridization` enum.

    Raises:
        ValueError: Could not determine hybridization for given atom.
    """
    hybridisation = oechem.OEGetHybridization(atom)
    if hybridisation != oechem.OEHybridization_Unknown:
        return hybridisation  # type: ignore
    # -[N+]#[N-]
    if (
        atom.GetAtomicNum() == 7
        and atom.GetFormalCharge() == -1
        and atom.GetDegree() == 2
        and atom.GetValence() == 4
    ):
        return oechem.OEHybridization_sp  # type: ignore
    # -[O-] (hydroxyl)
    if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == -1 and atom.GetValence() == 1:
        nbrs = list(atom.GetAtoms())
        if len(nbrs) == 0:
            return oechem.OEHybridization_sp3  # type: ignore
        neigh_hyb = oechem.OEGetHybridization(nbrs[0])
        return (  # type: ignore
            oechem.OEHybridization_sp2
            if neigh_hyb == oechem.OEHybridization_sp2
            else oechem.OEHybridization_sp3
        )
    raise ValueError(f"Failed to calculate hybridization for atom {atom}.")


def get_num_sp_lone_pairs_on_atom(atom: oechem.OEAtomBase, /) -> int:
    """
    Calculate the number of lone pairs on an atom that lie in s or sp orbitals.

    Args:
        atom: oechem.OEAtomBase

    Returns:
        Number of lone pairs.
    """
    lone_pairs = (
        OUTER_ELECTRON_COUNT[atom.GetAtomicNum()] - atom.GetFormalCharge() - atom.GetValence()
    ) // 2
    hybridisation = get_hybridization(atom)
    sp_orbitals = HYBRIDISATION_TO_SP_ORBITALS[hybridisation]
    if sp_orbitals > 0:
        # Hydrogens/bonds may take up sp orbitals
        free_sp_orbitals = sp_orbitals - atom.GetDegree()
        lone_pairs = min(free_sp_orbitals, lone_pairs)

    return lone_pairs  # type: ignore
