from openeye import oechem

IMPLICIT_H_COUNTS = {
    5: 3,  # Boron
    6: 4,  # Carbon
    7: 3,  # Nitrogen
    8: 2,  # Oxygen
    9: 1,  # Fluorine
    15: 3,  # Phosphorus
    16: 2,  # Sulfur
    17: 1,  # Chlorine
    35: 1,  # Bromine
    53: 1,  # Iodine
}


def modify_residue(
    residue: oechem.OEMolBase | oechem.OEAtomBase | oechem.OEResidue | oechem.OEAtomBondSet,
    *,
    residue_name: str | None = None,
    occupancy: float | None = None,
    b_factor: float | None = None,
    residue_number: int | None = None,
    sequence_number: int | None = None,
    serial_number: int | None = None,
    alternate_location: str | None = None,
    chain_id: str | None = None,
    asymmetric_unit_id: str | None = None,
    insertion_code: str | None = None,
    is_hetatm: bool | None = None,
) -> None:
    """Modify the residue information of an OpenEye object in place."""
    if isinstance(residue, oechem.OEMolBase | oechem.OEAtomBondSet):
        for atom in residue.GetAtoms():
            modify_residue(
                atom,
                residue_name=residue_name,
                occupancy=occupancy,
                b_factor=b_factor,
                residue_number=residue_number,
                sequence_number=sequence_number,
                serial_number=serial_number,
                alternate_location=alternate_location,
                chain_id=chain_id,
                asymmetric_unit_id=asymmetric_unit_id,
                insertion_code=insertion_code,
                is_hetatm=is_hetatm,
            )
    elif isinstance(residue, oechem.OEAtomBase):
        res = oechem.OEAtomGetResidue(residue)
        modify_residue(
            res,
            residue_name=residue_name,
            occupancy=occupancy,
            b_factor=b_factor,
            residue_number=residue_number,
            sequence_number=sequence_number,
            serial_number=serial_number,
            alternate_location=alternate_location,
            chain_id=chain_id,
            asymmetric_unit_id=asymmetric_unit_id,
            insertion_code=insertion_code,
            is_hetatm=is_hetatm,
        )
        oechem.OEAtomSetResidue(residue, res)
    else:
        if residue_name is not None:
            residue.SetName(residue_name)
        if occupancy is not None:
            residue.SetOccupancy(occupancy)
        if b_factor is not None:
            residue.SetBFactor(b_factor)
        if residue_number is not None:
            residue.SetResidueNumber(residue_number)
        if sequence_number is not None:
            residue.SetSequenceID(sequence_number)
        if serial_number is not None:
            residue.SetSerialNumber(serial_number)
        if alternate_location is not None:
            residue.SetAlternateLocation(alternate_location)
        if chain_id is not None:
            residue.SetExtChainID(chain_id)
        if asymmetric_unit_id is not None:
            residue.SetSubChainID(asymmetric_unit_id)
        if insertion_code is not None:
            residue.SetInsertCode(insertion_code)
        if is_hetatm is not None:
            residue.SetHetAtom(is_hetatm)


def create_residue(
    residue: oechem.OEResidue | None = None,
    residue_name: str | None = None,
    occupancy: float | None = None,
    b_factor: float | None = None,
    residue_number: int | None = None,
    sequence_number: int | None = None,
    serial_number: int | None = None,
    alternate_location: str | None = None,
    chain_id: str | None = None,
    asymmetric_unit_id: str | None = None,
    insertion_code: str | None = None,
    is_hetatm: bool | None = None,
) -> oechem.OEResidue:
    """
    Create a new `oechem.OEResidue` object with the specified properties.

    Args:
        residue: Preexisting residue to use as a base.
        residue_name: Residue name, such as `ALA`.
        occupancy: Occupancy of the atom described by this residue.
        b_factor: B-factor of the atom described by this residue.
        residue_number: Residue number as found in a PDB.
        sequence_number: Sequence number as found in an mmCIF.
        serial_number: Serial number of the atom described by this residue, as found in a PDB.
        alternate_location: Alternate location of the atom described by this residue, as found in a PDB.
        chain_id: Chain ID as found in a PDB. Corresponds to OpenEye's ExtChainID.
        asymmetric_unit_id: Asymmetric unit ID as found in an mmCIF. Corresponds to OpenEye's SubChainID.
        insertion_code: Insertion code of the residue, as found in a PDB.
        is_hetatm: Whether the atom described by this residue is a HETATM.

    Returns:
        `OEResidue` object with this properties.
    """
    residue = oechem.OEResidue(residue) if residue else oechem.OEResidue()
    modify_residue(
        residue,
        residue_name=residue_name,
        occupancy=occupancy,
        b_factor=b_factor,
        residue_number=residue_number,
        sequence_number=sequence_number,
        serial_number=serial_number,
        alternate_location=alternate_location,
        chain_id=chain_id,
        asymmetric_unit_id=asymmetric_unit_id,
        insertion_code=insertion_code,
        is_hetatm=is_hetatm,
    )
    return residue


def create_atom(
    mol: oechem.OEMolBase,
    /,
    *,
    atomic_number: int | None = None,
    atom: oechem.OEAtomBase | None = None,
    coordinates: tuple[float, float, float] | list[float] | None = None,
    name: str | None = None,
    isotope: int | None = None,
    formal_charge: int | None = None,
    partial_charge: float | None = None,
    implicit_h_count: int | None = None,
    residue: oechem.OEResidue | None = None,
    residue_name: str | None = None,
    occupancy: float | None = None,
    b_factor: float | None = None,
    residue_number: int | None = None,
    sequence_number: int | None = None,
    serial_number: int | None = None,
    alternate_location: str | None = None,
    chain_id: str | None = None,
    asymmetric_unit_id: str | None = None,
    insertion_code: str | None = None,
    is_hetatm: bool | None = None,
) -> oechem.OEAtomBase:
    """
    Create a new atom in a molecule with the specified properties.

    Args:
        mol: Molecule to add the atom to.
        atomic_number: Atomic number, such as 6 for carbon.
        atom: Prexisting atom to use as a base.
        coordinates: Coordinates of the atom.
        name: Name of the atom.
        isotope: Isotope number.
        formal_charge: Formal charge.
        partial_charge: Partial charge.
        implicit_h_count: Number of implicit hydrogen atoms on this atom.
        residue: Preexisting residue to use as a base.
        residue_name: Residue name, such as `ALA`.
        occupancy: Occupancy of the atom.
        b_factor: B-factor of the atom.
        residue_number: Residue number as found in a PDB.
        sequence_number: Sequence number as found in an mmCIF.
        serial_number: Serial number of the atom, as found in a PDB.
        alternate_location: Alternate location of the atom, as found in a PDB or mmCIF.
        chain_id: Chain ID as found in a PDB. Corresponds to OpenEye's ExtChainID.
        asymmetric_unit_id: Asymmetric unit ID as found in an mmCIF. Corresponds to OpenEye's SubChainID.
        insertion_code: Insertion code of the residue, as found in a PDB.
        is_hetatm: Whether the atom is a HETATM.

    Returns:
        `OEAtomBase` that has been added to the molecule.
    """
    if atomic_number is None and atom is None:
        raise ValueError("Either `atom` or `atomic_number` must be specified.")
    if atom is not None:
        atom = mol.NewAtom(atom)
        if atomic_number is not None:
            atom.SetAtomicNum(atomic_number)
    else:
        atom = mol.NewAtom(atomic_number)
    if name is not None:
        atom.SetName(name)
    if isotope is not None:
        atom.SetIsotope(isotope)
    if coordinates:
        if atom.GetParent().GetDimension() == 0:
            atom.GetParent().SetDimension(len(coordinates))
        atom.GetParent().SetCoords(atom, coordinates)
    if formal_charge is not None:
        atom.SetFormalCharge(formal_charge)
    if partial_charge is not None:
        atom.SetPartialCharge(partial_charge)
    if implicit_h_count is not None:
        atom.SetImplicitHCount(implicit_h_count)
    else:
        # Default number of implicit hydrogens
        if atomic_number in IMPLICIT_H_COUNTS:
            atom.SetImplicitHCount(IMPLICIT_H_COUNTS[atomic_number])
    residue = create_residue(
        residue=residue,
        residue_name=residue_name,
        occupancy=occupancy,
        b_factor=b_factor,
        residue_number=residue_number,
        sequence_number=sequence_number,
        serial_number=serial_number,
        alternate_location=alternate_location,
        chain_id=chain_id,
        asymmetric_unit_id=asymmetric_unit_id,
        insertion_code=insertion_code,
        is_hetatm=is_hetatm,
    )
    oechem.OEAtomSetResidue(atom, residue)
    return atom


def create_bond(
    mol: oechem.OEMolBase,
    /,
    *,
    atom1: oechem.OEAtomBase,
    atom2: oechem.OEAtomBase,
    order: int | None = None,
    adjust_h_count: bool = True,
    adjust_charge: bool = True,
) -> oechem.OEBondBase:
    """
    Create a new bond in a molecule with the specified properties.

    By default, hydrogen counts and charges can be adjusted (both `adjust_h_count` and `adjust_charge` are `True`).
    The behaviour of this function depends on these two:

    * If both `adjust_h_count` and `adjust_charge` are `True`, then for each atom, hydrogens are removed to satisfy the
      bond orders. If hydrogens run out, then the formal charge is increased.
    * If just `adjust_h_count` is `True`, then hydrogens are removed to satisfy the bond orders. If hydrogens run out,
      an error is raised.
    * If just `adjust_charge` is `True`, then the charges of each atom are increased to satisfy the bond orders.
    * If neither is `True`, the bond is added without modifying either atom.

    Args:
        mol: Molecule to add the bond to.
        atom1: First atom in the bond.
        atom2: Second atom in the bond.
        order: Bond order.
        adjust_h_count: Whether to adjust hydrogen counts to account for the new bond.
        adjust_charge: Whether to adjust formal charges to account for the new bond.

    Returns:
        `OEBondBase` that has been added to the molecule.
    """
    if adjust_h_count and (atom1.GetExplicitHCount() > 0 or atom2.GetExplicitHCount() > 0):
        raise ValueError("Cannot create a bond between atoms with explicit hydrogens.")

    order = order or 1

    hydrogens1_to_remove = min(order, atom1.GetImplicitHCount()) if adjust_h_count else 0
    hydrogens2_to_remove = min(order, atom2.GetImplicitHCount()) if adjust_h_count else 0

    if (
        adjust_h_count
        and not adjust_charge
        and ((hydrogens1_to_remove != order) or (hydrogens2_to_remove != order))
    ):
        raise ValueError("Atoms must have enough implicit hydrogens to create a bond.")

    charge1_to_increase = order - hydrogens1_to_remove if adjust_charge else 0
    charge2_to_increase = order - hydrogens2_to_remove if adjust_charge else 0

    bond = mol.NewBond(atom1, atom2, order)
    if hydrogens1_to_remove:
        atom1.SetImplicitHCount(atom1.GetImplicitHCount() - hydrogens1_to_remove)
    if hydrogens2_to_remove:
        atom2.SetImplicitHCount(atom2.GetImplicitHCount() - hydrogens2_to_remove)

    if charge1_to_increase:
        atom1.SetFormalCharge(atom1.GetFormalCharge() + charge1_to_increase)
    if charge2_to_increase:
        atom2.SetFormalCharge(atom2.GetFormalCharge() + charge2_to_increase)

    return bond
