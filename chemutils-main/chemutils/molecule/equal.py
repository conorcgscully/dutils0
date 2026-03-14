import numpy as np
from openeye import oechem

# By default, assert coordinate equality to 0.001 angstroms
DEFAULT_COORD_ATOL = 0.001
# By default, assert occupancy equality to 0.001
DEFAULT_OCCUPANCY_ATOL = 0.001

from enum import Flag, auto
from typing import Final

from chemutils.molecule.stereo import (
    get_atom_tetrahedral_stereo,
    get_bond_cistrans_stereo,
)


class AtomProperty(Flag):
    Nothing = 0
    Index = auto()
    AtomicNumber = auto()
    HeavyDegree = auto()
    Degree = auto()
    HCount = auto()
    ImplicitHCount = auto()
    ExplicitHCount = auto()
    FormalCharge = auto()
    PartialCharge = auto()
    Hybridisation = auto()
    Isotope = auto()
    Name = auto()
    Valence = auto()
    Aromatic = auto()
    Chiral = auto()
    IsInRing = auto()
    Coordinates = auto()
    Stereo = auto()
    ResidueName = auto()
    AlternateLocation = auto()
    BFactor = auto()
    ChainID = auto()
    AsymmetricUnitID = auto()
    InsertionCode = auto()
    Occupancy = auto()
    ResidueNumber = auto()
    SerialNumber = auto()
    SequenceNumber = auto()
    HETATM = auto()


class BondProperty(Flag):
    Nothing = 0
    Index = auto()
    BeginIndex = auto()
    EndIndex = auto()
    Aromatic = auto()
    Chiral = auto()
    IsInRing = auto()
    Stereo = auto()
    Type = auto()
    Order = auto()
    IsRotor = auto()


AllBondProperties: Final = (
    BondProperty.Index
    | BondProperty.BeginIndex
    | BondProperty.EndIndex
    | BondProperty.Aromatic
    | BondProperty.Chiral
    | BondProperty.IsInRing
    | BondProperty.Stereo
    | BondProperty.Type
    | BondProperty.IsRotor
    | BondProperty.Order
)


AllAtomProperties: Final = (
    AtomProperty.Index
    | AtomProperty.AtomicNumber
    | AtomProperty.Degree
    | AtomProperty.ImplicitHCount
    | AtomProperty.ExplicitHCount
    | AtomProperty.FormalCharge
    | AtomProperty.PartialCharge
    | AtomProperty.Hybridisation
    | AtomProperty.Isotope
    | AtomProperty.Name
    | AtomProperty.Valence
    | AtomProperty.Aromatic
    | AtomProperty.Chiral
    | AtomProperty.IsInRing
    | AtomProperty.Coordinates
    | AtomProperty.Stereo
    | AtomProperty.ResidueName
    | AtomProperty.AlternateLocation
    | AtomProperty.BFactor
    | AtomProperty.ChainID
    | AtomProperty.AsymmetricUnitID
    | AtomProperty.InsertionCode
    | AtomProperty.Occupancy
    | AtomProperty.ResidueNumber
    | AtomProperty.SerialNumber
    | AtomProperty.SequenceNumber
    | AtomProperty.HETATM
)


def assert_atom_equal(
    a: oechem.OEAtomBase,
    b: oechem.OEAtomBase,
    /,
    *,
    coord_atol: float = DEFAULT_COORD_ATOL,
    occupancy_atol: float = DEFAULT_OCCUPANCY_ATOL,
    properties: AtomProperty = AllAtomProperties,
) -> None:
    """
    Assert that two atoms are equal.

    Args:
        a: The first atom to compare.
        b: The second atom to compare.
        coord_atol: Absolute tolerance for coordinate equality.
        occupancy_atol: Absolute tolerance for occupancy equality.
        properties: Flags for which properties to compare (by default, all).

    Raises:
        AssertionError: The two atoms are not equal.
    """
    if AtomProperty.Index in properties:
        assert a.GetIdx() == b.GetIdx(), f"Atom index unequal: {a.GetIdx()} != {b.GetIdx()}"
    if AtomProperty.AtomicNumber in properties:
        assert a.GetAtomicNum() == b.GetAtomicNum(), (
            f"Atomic number unequal: {a.GetAtomicNum()} != {b.GetAtomicNum()}"
        )
    if AtomProperty.HeavyDegree in properties:
        assert a.GetHvyDegree() == b.GetHvyDegree(), (
            f"Atom heavy degree unequal: {a.GetHvyDegree()} != {b.GetHvyDegree()}"
        )
    if AtomProperty.Degree in properties:
        assert a.GetDegree() == b.GetDegree(), (
            f"Atom degree unequal: {a.GetDegree()} != {b.GetDegree()}"
        )
    if AtomProperty.HCount in properties:
        assert a.GetHCount() == b.GetHCount(), (
            f"Atom H count unequal: {a.GetHCount()} != {b.GetHCount()}"
        )
    if AtomProperty.ImplicitHCount in properties:
        assert a.GetImplicitHCount() == b.GetImplicitHCount(), (
            f"Atom implicit H count unequal: {a.GetImplicitHCount()} != {b.GetImplicitHCount()}"
        )
    if AtomProperty.ExplicitHCount in properties:
        assert a.GetExplicitHCount() == b.GetExplicitHCount(), (
            f"Atom explicit H count unequal: {a.GetExplicitHCount()} != {b.GetExplicitHCount()}"
        )
    if AtomProperty.FormalCharge in properties:
        assert a.GetFormalCharge() == b.GetFormalCharge(), (
            f"Atom formal charge unequal: {a.GetFormalCharge()} != {b.GetFormalCharge()}"
        )
    if AtomProperty.PartialCharge in properties:
        assert a.GetPartialCharge() == b.GetPartialCharge(), (
            f"Atom partial charge unequal: {a.GetPartialCharge()} != {b.GetPartialCharge()}"
        )
    if AtomProperty.Hybridisation in properties:
        assert a.GetHyb() == b.GetHyb(), f"Atom hybridisation unequal: {a.GetHyb()} != {b.GetHyb()}"
    if AtomProperty.Isotope in properties:
        assert a.GetIsotope() == b.GetIsotope(), (
            f"Atom isotope unequal: {a.GetIsotope()} != {b.GetIsotope()}"
        )
    if AtomProperty.Name in properties:
        assert a.GetName() == b.GetName(), f"Atom name unequal: {a.GetName()} != {b.GetName()}"
    if AtomProperty.Valence in properties:
        assert a.GetValence() == b.GetValence(), (
            f"Atom valence unequal: {a.GetValence()} != {b.GetValence()}"
        )
    if AtomProperty.Aromatic in properties:
        assert a.IsAromatic() == b.IsAromatic(), (
            f"Atom aromaticity unequal: {a.IsAromatic()} != {b.IsAromatic()}"
        )
    if AtomProperty.Chiral in properties:
        assert a.IsChiral() == b.IsChiral(), (
            f"Atom chirality unequal: {a.IsChiral()} != {b.IsChiral()}"
        )
    if AtomProperty.IsInRing in properties:
        assert a.IsInRing() == b.IsInRing(), (
            f"Atom ring membership unequal: {a.IsInRing()} != {b.IsInRing()}"
        )
    if AtomProperty.Stereo in properties:
        stereo_a = get_atom_tetrahedral_stereo(a)
        stereo_b = get_atom_tetrahedral_stereo(b)
        assert stereo_a == stereo_b, f"Atom tetrahedral stereo unequal: {stereo_a} != {stereo_b}"
    if AtomProperty.Coordinates in properties:
        a_coords = np.array(a.GetParent().GetCoords(a))
        b_coords = np.array(b.GetParent().GetCoords(b))
        # Check coordinate equality to 0.01 angstroms
        assert np.allclose(a_coords, b_coords, atol=coord_atol), (
            f"Atom coordinates unequal: {a_coords} != {b_coords}"
        )

    # Residue-related properties
    residue1: oechem.OEResidue = oechem.OEAtomGetResidue(a)
    residue2: oechem.OEResidue = oechem.OEAtomGetResidue(b)
    if AtomProperty.ResidueName in properties:
        assert residue1.GetName() == residue2.GetName(), (
            f"Residue name unequal: {residue1.GetName()} != {residue2.GetName()}"
        )
    if AtomProperty.AlternateLocation in properties:
        assert residue1.GetAlternateLocation() == residue2.GetAlternateLocation(), (
            f"Alternate location unequal: {residue1.GetAlternateLocation()} != {residue2.GetAlternateLocation()}"
        )
    if AtomProperty.BFactor in properties:
        assert residue1.GetBFactor() == residue2.GetBFactor(), (
            f"B factor unequal: {residue1.GetBFactor()} != {residue2.GetBFactor()}"
        )
    if AtomProperty.ChainID in properties:
        assert residue1.GetExtChainID() == residue2.GetExtChainID(), (
            f"Chain ID unequal: {residue1.GetExtChainID()} != {residue2.GetExtChainID()}"
        )
    if AtomProperty.AsymmetricUnitID in properties:
        assert residue1.GetSubChainID() == residue2.GetSubChainID(), (
            f"Asymmetric Unit ID unequal: {residue1.GetSubChainID()} != {residue2.GetSubChainID()}"
        )
    if AtomProperty.InsertionCode in properties:
        assert residue1.GetInsertCode() == residue2.GetInsertCode(), (
            f"Insert code unequal: {residue1.GetInsertCode()} != {residue2.GetInsertCode()}"
        )
    if AtomProperty.Occupancy in properties:
        try:
            np.testing.assert_allclose(
                residue1.GetOccupancy(), residue2.GetOccupancy(), atol=occupancy_atol
            )
        except AssertionError:
            raise AssertionError(
                f"Occupancy unequal: {residue1.GetOccupancy()} != {residue2.GetOccupancy()}"
            ) from None
    if AtomProperty.ResidueNumber in properties:
        assert residue1.GetResidueNumber() == residue2.GetResidueNumber(), (
            f"Residue number unequal: {residue1.GetResidueNumber()} != {residue2.GetResidueNumber()}"
        )
    if AtomProperty.SerialNumber in properties:
        assert residue1.GetSerialNumber() == residue2.GetSerialNumber(), (
            f"Serial number unequal: {residue1.GetSerialNumber()} != {residue2.GetSerialNumber()}"
        )
    if AtomProperty.SequenceNumber in properties:
        assert residue1.GetSequenceID() == residue2.GetSequenceID(), (
            f"Sequence ID unequal: {residue1.GetSequenceID()} != {residue2.GetSequenceID()}"
        )
    if AtomProperty.HETATM in properties:
        assert residue1.IsHetAtom() == residue2.IsHetAtom(), (
            f"HETATM unequal: {residue1.IsHetAtom()} != {residue2.IsHetAtom()}"
        )


def assert_bond_equal(
    a: oechem.OEBondBase,
    b: oechem.OEBondBase,
    /,
    *,
    ignore_kekulization: bool = True,
    properties: BondProperty = AllBondProperties,
) -> None:
    """
    Assert that two bonds are equal.

    Args:
        a: The first bond to compare.
        b: The second bond to compare.
        ignore_kekulization: Should bond order be ignored when comparing aromatic bonds?
        properties: Flags for which properties to compare (by default, all).

    Raises:
        AssertionError: The two bonds are not equal.
    """
    if BondProperty.Index in properties:
        assert a.GetIdx() == b.GetIdx(), f"Bond index unequal: {a.GetIdx()} != {b.GetIdx()}"
    if BondProperty.BeginIndex in properties:
        assert a.GetBgnIdx() == b.GetBgnIdx(), (
            f"Bond start index unequal: {a.GetBgnIdx()} != {b.GetBgnIdx()}"
        )
    if BondProperty.EndIndex in properties:
        assert a.GetEndIdx() == b.GetEndIdx(), (
            f"Bond end index unequal: {a.GetEndIdx()} != {b.GetEndIdx()}"
        )
    if BondProperty.Type in properties:
        assert a.GetType() == b.GetType(), f"Bond type unequal: {a.GetType()} != {b.GetType()}"
    if BondProperty.Aromatic in properties:
        assert a.IsAromatic() == b.IsAromatic(), (
            f"Bond aromaticity unequal: {a.IsAromatic()} != {b.IsAromatic()}"
        )
    if BondProperty.Order in properties and (
        not ignore_kekulization or (not a.IsAromatic() and not b.IsAromatic())
    ):
        assert a.GetOrder() == b.GetOrder(), f"Bond order unequal: {a.GetOrder()} != {b.GetOrder()}"
    if BondProperty.Chiral in properties:
        assert a.IsChiral() == b.IsChiral(), (
            f"Bond chirality unequal: {a.IsChiral()} != {b.IsChiral()}"
        )
    if BondProperty.IsInRing in properties:
        assert a.IsInRing() == b.IsInRing(), (
            f"Bond ring membership unequal: {a.IsInRing()} != {b.IsInRing()}"
        )
    if BondProperty.IsRotor in properties:
        assert a.IsRotor() == b.IsRotor(), f"Bond rotatable unequal: {a.IsRotor()} != {b.IsRotor()}"

    stereo_a = get_bond_cistrans_stereo(a)
    stereo_b = get_bond_cistrans_stereo(b)
    if BondProperty.Stereo in properties:
        assert stereo_a == stereo_b, f"Bond stereo unequal: {stereo_a} != {stereo_b}"


def assert_mol_equal(
    a: oechem.OEMolBase,
    b: oechem.OEMolBase,
    /,
    *,
    ignore_kekulization: bool = True,
    coord_atol: float = DEFAULT_COORD_ATOL,
    occupancy_atol: float = DEFAULT_OCCUPANCY_ATOL,
    atom_properties: AtomProperty = AllAtomProperties,
    bond_properties: BondProperty = AllBondProperties,
) -> None:
    """
    Assert that two molecules are equal.

    Args:
        a: The first molecule to compare.
        b: The second molecule to compare.
        ignore_kekulization: Should bond order be ignored when comparing aromatic bonds?
        coord_atol: Absolute tolerance for coordinate equality.
        occupancy_atol: Absolute tolerance for occupancy equality.
        atom_properties: Flags indicating which properties to compare atoms on.
        bond_properties: Flags indicating which properties to compare bonds on.

    Raises:
        AssertionError: The two molecules are not equal.
    """
    assert a.NumAtoms() == b.NumAtoms(), (
        f"Molecule atom count unequal: {a.NumAtoms()} != {b.NumAtoms()}"
    )
    assert a.NumBonds() == b.NumBonds(), (
        f"Molecule atom count unequal: {a.NumAtoms()} != {b.NumAtoms()}"
    )
    for i, atom_a, atom_b in zip(range(a.NumAtoms()), a.GetAtoms(), b.GetAtoms(), strict=True):
        try:
            assert_atom_equal(
                atom_a,
                atom_b,
                coord_atol=coord_atol,
                occupancy_atol=occupancy_atol,
                properties=atom_properties,
            )
        except AssertionError as e:
            raise AssertionError(f"Atoms with index {i} are not equal") from e
    for i, bond_a, bond_b in zip(range(a.NumBonds()), a.GetBonds(), b.GetBonds(), strict=True):
        try:
            assert_bond_equal(
                bond_a, bond_b, ignore_kekulization=ignore_kekulization, properties=bond_properties
            )
        except AssertionError as e:
            raise AssertionError(f"Bonds with index {i} are not equal") from e


def atom_equal(
    a: oechem.OEAtomBase,
    b: oechem.OEAtomBase,
    /,
    coord_atol: float = DEFAULT_COORD_ATOL,
    occupancy_atol: float = DEFAULT_OCCUPANCY_ATOL,
    properties: AtomProperty = AllAtomProperties,
) -> bool:
    """
    Check if two OpenEye atoms are equal.

    This checks value equality, which allows comparison between atoms belonging to two
    different `OEMolBase` objects. The atoms must be exactly the same, including their
    index and atomic properties.

    Args:
        a: The first atom to compare.
        b: The second atom to compare.
        coord_atol: Absolute tolerance for coordinate equality.
        occupancy_atol: Absolute tolerance for occupancy equality.
        properties: Flags for which properties to compare (by default, all).

    Returns:
        True if the two atoms are equal, and False otherwise.
    """
    try:
        assert_atom_equal(
            a, b, coord_atol=coord_atol, occupancy_atol=occupancy_atol, properties=properties
        )
        return True
    except AssertionError:
        return False


def bond_equal(
    a: oechem.OEBondBase,
    b: oechem.OEBondBase,
    /,
    *,
    ignore_kekulization: bool = True,
    properties: BondProperty = AllBondProperties,
) -> bool:
    """
    Check if two OpenEye bonds are equal.

    This checks value equality, which allows comparison between bonds belonging to two
    different `OEMolBase` objects. The bonds must be exactly the same, including their
    index, atom indices and bond properties.

    Args:
        a: The first bond to compare.
        b: The second bond to compare.
        ignore_kekulization: Should bond order be ignored when comparing aromatic bonds?
        properties: Flags for which properties to compare (by default, all).

    Returns:
        True if the two bonds are equal, and False otherwise.
    """
    try:
        assert_bond_equal(a, b, ignore_kekulization=ignore_kekulization, properties=properties)
        return True
    except AssertionError:
        return False


def mol_equal(
    a: oechem.OEMolBase,
    b: oechem.OEMolBase,
    /,
    *,
    ignore_kekulization: bool = True,
    coord_atol: float = DEFAULT_COORD_ATOL,
    occupancy_atol: float = DEFAULT_OCCUPANCY_ATOL,
    atom_properties: AtomProperty = AllAtomProperties,
    bond_properties: BondProperty = AllBondProperties,
) -> bool:
    """
    Check if two OpenEye molecules are equal.

    This checks that the number of atoms and bonds agree, are in the same order and
    match exactly.

    Args:
        a: The first molecule to compare.
        b: The second molecule to compare.
        ignore_kekulization: Should bond order be ignored when comparing aromatic bonds?
        coord_atol: Absolute tolerance for coordinate equality.
        occupancy_atol: Absolute tolerance for occupancy equality.
        atom_properties: Flags indicating which properties to compare atoms on.
        bond_properties: Flags indicating which properties to compare bonds on.

    Returns:
        True if the two molecules are equal, and False otherwise.
    """
    try:
        assert_mol_equal(
            a,
            b,
            ignore_kekulization=ignore_kekulization,
            coord_atol=coord_atol,
            occupancy_atol=occupancy_atol,
            atom_properties=atom_properties,
            bond_properties=bond_properties,
        )
        return True
    except AssertionError:
        return False
