from openeye import oechem

VALENCE_ELECTRONS = {
    1: 1,
    5: 3,
    6: 4,
    7: 5,
    8: 6,
    9: 7,
    14: 4,
    15: 5,
    16: 6,
    17: 7,
    34: 6,
    35: 7,
    53: 7,
}

AZIDO_CHARGE_SEPARATION = oechem.OEUniMolecularRxn(
    "[*:1][N:2]=[N+:3]=[NH:4]>>[*:1][N:2]=[N+:3]=[N-:4]"
)
# Nitro and nitrone
NITRO_CHARGE_SEPARATION = oechem.OEUniMolecularRxn(
    "[*+0:1]=[N+:2]([*+0:3])[OH:4]>>[*+0:1]=[N+:2]([*+0:3])[O-:4]"
)
AROMATIC_NITRO_CHARGE_SEPARATION = oechem.OEUniMolecularRxn(
    "[*+0:1]:[n+:2](:[*+0:3])[OH:4]>>[*+0:1]:[n+:2](:[*+0:3])[O-:4]"
)
ISONITRILE_CHARGE_SEPARATION = oechem.OEUniMolecularRxn("[*:1][N+:2]#[CH:3]>>[*:1][N+:2]#[C-:3]")


def _set_charge_and_h_count(*, atom: oechem.OEAtomBase, charge: int, h_count: int) -> None:
    total_h_to_remove = atom.GetTotalHCount() - h_count
    implicit_h_count = atom.GetImplicitHCount()
    explicit_h_count = atom.GetExplicitHCount()
    h_count = implicit_h_count + explicit_h_count
    if implicit_h_count >= total_h_to_remove:
        atom.SetImplicitHCount(atom.GetImplicitHCount() - total_h_to_remove)
    else:
        atom.SetImplicitHCount(0)
        total_h_to_remove -= implicit_h_count
        nbrs = [nbr for nbr in atom.GetAtoms() if nbr.GetAtomicNum() == 1]
        for to_remove in nbrs[:total_h_to_remove]:
            atom.GetParent().DeleteAtom(to_remove)
    atom.SetFormalCharge(charge)


def _neutralize_atom(atom: oechem.OEAtomBase, /) -> None:
    atomic_number = atom.GetAtomicNum()
    if atomic_number not in VALENCE_ELECTRONS:
        return

    valence_electrons = VALENCE_ELECTRONS[atomic_number]

    bonded_valence_electrons = atom.GetHvyValence()

    if bonded_valence_electrons >= valence_electrons:
        _set_charge_and_h_count(
            atom=atom, charge=valence_electrons - bonded_valence_electrons, h_count=0
        )
        return

    lowest_valence = 4 - abs(4 - valence_electrons)

    if bonded_valence_electrons < lowest_valence:
        _set_charge_and_h_count(
            atom=atom, charge=0, h_count=lowest_valence - bonded_valence_electrons
        )
        return

    unbonded_valence_electrons = bonded_valence_electrons - lowest_valence

    _set_charge_and_h_count(atom=atom, charge=unbonded_valence_electrons % 2, h_count=0)


def neutralize_molecule(mol: oechem.OEMolBase, /) -> None:
    """
    Neutralize a molecule to a standard form, adjusting charges and hydrogen counts.

    This technique effectively strips all hydrogens and charges, and then reassigns them
    based on the following rules:

    * If an atom has a standard valency, then it is left uncharged with no hydrogens.
    * If adding hydrogens would bring it up to a standard valency, then it is left uncharged and hydrogens are added.
    * Otherwise, a positive charge is assigned to the atom.

    Finally, a few standard charge-separation functional groups are corrected.

    * Azido groups are converted from `N=[N+]=[NH]` to `N=[N+]=[N-]`
    * Nitro/nitrone groups are converted from `O=[N+][OH]` to `O=[N+][O-]`
    * Isonitrile groups are converted from `-[N+]#[CH]` to `-[N+]#[C-]`

    Args:
        mol (oechem.OEMolBase): Molecule to neutralize.
    """
    for atom in mol.GetAtoms():
        _neutralize_atom(atom)
    AZIDO_CHARGE_SEPARATION(mol)
    NITRO_CHARGE_SEPARATION(mol)
    AROMATIC_NITRO_CHARGE_SEPARATION(mol)
    ISONITRILE_CHARGE_SEPARATION(mol)
