import itertools
from collections import defaultdict
from collections.abc import Generator

import numpy as np
import numpy.typing as npt
from openeye import oechem


class UnderspecifiedStereoError(Exception):
    def __init__(self, oemol: oechem.OEMolBase):
        super().__init__(f"Stereochemistry is underspecified for {oechem.OEMolToSmiles(oemol)}")


def is_specific_stereoisomer(oemol: oechem.OEMolBase, /, include_nitrogens: bool = True) -> bool:
    """
    Does the given molecule have stereochemistry defined for all of its chiral atoms.

    Args:
        oemol: OpenEye molecule to consider.
        include_nitrogens: Wether nitrogen stereocenters count. default = True

    Returns:
        True if every chiral atom is assigned either clockwise or anticlockwise, False otherwise.
    """
    stereoatoms, stereobonds = get_stereocenters_and_stereobonds(oemol)
    for stereocenter in stereoatoms:
        if not include_nitrogens and can_invert(stereocenter):
            continue
        if get_atom_tetrahedral_stereo(stereocenter) == oechem.OEAtomStereo_Undefined:
            return False
    for stereobond in stereobonds:
        if get_bond_cistrans_stereo(stereobond) == oechem.OEBondStereo_Undefined:
            return False
    return True


def reflect_tetrahedral_stereo(stereo_label: int, /) -> int:
    """
    Get the reflected version of an OpenEye tetrahedral stereochemistry label.

    This flips right-handed to left-handed and vice-versa, and does not flip undefined centers.

    Args:
        stereo_label: Tetrahedral stereolabel from the OEAtomStereo namespace.

    Returns:
        Flipped tetrahdral stereolabel from the OEAtomStereo namespace.

    Raises:
        ValueError: Unknown stereo label.
    """
    match stereo_label:
        case oechem.OEAtomStereo_LeftHanded:
            return oechem.OEAtomStereo_RightHanded  # type: ignore
        case oechem.OEAtomStereo_RightHanded:
            return oechem.OEAtomStereo_LeftHanded  # type: ignore
        case oechem.OEAtomStereo_Undefined:
            return oechem.OEAtomStereo_Undefined  # type: ignore
        case _:
            raise ValueError(f"Unknown Tetrahedral Stereo Label {stereo_label}")


def reflect_cistrans_stereo(stereo_label: int, /) -> int:
    """
    Get the reflected version of an OpenEye cis-trans stereochemistry label.

    This flips cis to trans and vice-versa, and does not flip undefined centers.

    Args:
        stereo_label: Cis-trans stereolabel from the OEBondStereo namespace.

    Returns:
        Flipped cis-trans stereolabel from the OEBondStereo namespace.

    Raises:
        ValueError: Unknown stereo label.
    """
    match stereo_label:
        case oechem.OEBondStereo_Undefined:
            return oechem.OEBondStereo_Undefined  # type: ignore
        case oechem.OEBondStereo_Cis:
            return oechem.OEBondStereo_Trans  # type: ignore
        case oechem.OEBondStereo_Trans:
            return oechem.OEBondStereo_Cis  # type: ignore
        case _:
            raise ValueError(f"Unknown Tetrahedral Stereo Label {stereo_label}")


def reflect_stereocenter(atom: oechem.OEAtomBase, /) -> None:
    """
    Reflect a given atom's stereochemistry label.

    If the atom lacks stereo information, this is a no-op.

    Args:
        atom: OpenEye atom to reflect.
    """
    if atom.HasStereoSpecified(oechem.OEAtomStereo_Tetrahedral):
        stereo_label = get_atom_tetrahedral_stereo(atom)
        stereo_label = reflect_tetrahedral_stereo(stereo_label)
        set_atom_tetrahedral_stereo(atom, stereo_label=stereo_label)


def reflect_stereobond(bond: oechem.OEBondBase, /) -> None:
    """
    Reflect a given bond's cis-trans label.

    If the bond lacks stereo information, this is a no-op.

    Args:
        bond: OpenEye bond to reflect.
    """
    if bond.HasStereoSpecified(oechem.OEBondStereo_CisTrans):
        stereo_label = get_bond_cistrans_stereo(bond)
        stereo_label = reflect_cistrans_stereo(stereo_label)
        set_bond_cistrans_stereo(bond, stereo_label=stereo_label)


def reflect_molecule_stereo(oemol: oechem.OEMolBase, /) -> oechem.OEMolBase:
    """
    Returns a new molecule by reflecting all chiral tetrahdral atoms.

    Chiral atoms with unspecified stereochemistry are unaffected.

    The returned molecule will have any coordinates that the original molecule had removed.

    Args:
        oemol: OpenEye molecule to reflect.

    Returns:
        Reflected OpenEye molecule.
    """
    copy = oemol.CreateCopy()
    copy.ClearCoords()
    for atom in copy.GetAtoms():
        reflect_stereocenter(atom)
    for bond in copy.GetBonds():
        reflect_stereobond(bond)
    return copy


def are_stereo_related(oemol1: oechem.OEMolBase, oemol2: oechem.OEMolBase, /) -> bool:
    """
    Are two OpenEye molecules equivalent ignoring stereochemistry (same structure, but different stereo labels).

    Identical molecules are considered to be stereo related. Molecules can also be stereo related even
    when they don't have stereochemistry specified.

    Args:
        oemol1: First OpenEye molecule to consider.
        oemol2: Second OpenEye molecule to consider.

    Returns:
        True if the two molecules are stereoisomers, and False otherwise.
    """
    return oechem.OECreateCanSmiString(oemol1) == oechem.OECreateCanSmiString(oemol2)  # type: ignore


def are_stereoisomers(oemol1: oechem.OEMolBase, oemol2: oechem.OEMolBase, /) -> bool:
    """
    Are two OpenEye molecules stereoisomers (same structure but differing stereo labels).

    By the IUPAC definition, two identical structures are not stereoisomers to each other.

    If either molecule does not have all its stereochemistry specified, then an `UnderspecifiedStereoError`
    is raised.

    Args:
        oemol1: First OpenEye molecule to consider.
        oemol2: Second OpenEye molecule to consider.

    Returns:
        True if the two molecules are stereoisomers, and False otherwise.

    Raises:
        UnderspecifiedStereoError: Either of the molecules does not have all its stereochemistry specified.
    """
    if not is_specific_stereoisomer(oemol1):
        raise UnderspecifiedStereoError(oemol1)
    if not is_specific_stereoisomer(oemol2):
        raise UnderspecifiedStereoError(oemol2)
    return are_stereo_related(oemol1, oemol2) and oechem.OEMolToSmiles(
        oemol1
    ) != oechem.OEMolToSmiles(oemol2)


def are_enantiomers(oemol1: oechem.OEMolBase, oemol2: oechem.OEMolBase, /) -> bool:
    """
    Are two OpenEye molecules enantiomers (stereoisomers which are mirror images).

    By the IUPAC definition, two identical molecules cannot be considered to be enantiomers, as they must both
    be mirror images and non-superimposable.

    If either molecule does not have all its stereochemistry specified, then an `UnderspecifiedStereoError`
    is raised.

    Args:
        oemol1: First OpenEye molecule to consider.
        oemol2: Second OpenEye molecule to consider.

    Returns:
        True if the two molecules are enantiomers, and False otherwise.

    Raises:
        UnderspecifiedStereoError: Either of the molecules does not have all its stereochemistry specified.
    """
    return are_stereoisomers(oemol1, oemol2) and are_mirror_images(oemol1, oemol2)


def are_diastereomers(oemol1: oechem.OEMolBase, oemol2: oechem.OEMolBase, /) -> bool:
    """
    Are two OpenEye molecules diastereomers (stereoisomers which are not mirror images).

    By the IUPAC definition, two identical molecules cannot be considered to be diastereomers, as they
    must be stereoisomers which must be unique.

    If either molecule does not have all its stereochemistry specified, then an `UnderspecifiedStereoError`
    is raised.

    Args:
        oemol1: First OpenEye molecule to consider.
        oemol2: Second OpenEye molecule to consider.

    Returns:
        True if the two molecules are diastereomers, and False otherwise.

    Raises:
        UnderspecifiedStereoError: Either of the molecules does not have all its stereochemistry specified.
    """
    return are_stereoisomers(oemol1, oemol2) and not are_mirror_images(oemol1, oemol2)


def are_mirror_images(oemol1: oechem.OEMolBase, oemol2: oechem.OEMolBase, /) -> bool:
    """
    Are two OpenEye molecules mirror images.

    Two identical molecules can be considered to be mirror images. This indicates that the molecule is achiral.

    If either molecule does not have all its stereochemistry specified, then an `UnderspecifiedStereoError`
    is raised.

    Args:
        oemol1: First OpenEye molecule to consider.
        oemol2: Second OpenEye molecule to consider.

    Returns:
        True if the two molecules are mirror images, and False otherwise.

    Raises:
        UnderspecifiedStereoError: Either of the molecules does not have all its stereochemistry specified.
    """
    if not is_specific_stereoisomer(oemol1):
        raise UnderspecifiedStereoError(oemol1)
    if not is_specific_stereoisomer(oemol2):
        raise UnderspecifiedStereoError(oemol2)
    reflected_oemol1 = reflect_molecule_stereo(oemol1)
    return oechem.OEMolToSmiles(oemol2) == oechem.OEMolToSmiles(reflected_oemol1)  # type: ignore


def is_meso_isomer(oemol: oechem.OEMolBase, /) -> bool:
    """
    Is the given molecule a meso (stereo)isomer, i.e. superimposable to its mirror image.

    If the molecule does not have all its stereochemistry specified, then an `UnderspecifiedStereoError`
    is raised.

    If the molecule has less than two stereocenters, returns False.

    Args:
        oemol: OpenEye molecule to consider.

    Returns:
        True if the molecule is a meso isomer, and False otherwise.

    Raises:
        UnderspecifiedStereoError:The molecule does not have all its stereochemistry specified.
    """
    if sum(atom.IsChiral() for atom in oemol.GetAtoms()) < 2:
        return False
    return not is_chiral_molecule(oemol)


def is_chiral_molecule(oemol: oechem.OEMolBase, /) -> bool:
    """
    Is the given molecule not superimposable to its mirror image.

    If the molecule does not have all its stereochemistry specified, then an `UnderspecifiedStereoError`
    is raised.

    Args:
        oemol: OpenEye molecule to consider.

    Returns:
        True if the molecule is a chiral, and False otherwise.

    Raises:
        UnderspecifiedStereoError:The molecule does not have all its stereochemistry specified.
    """
    if not is_specific_stereoisomer(oemol):
        raise UnderspecifiedStereoError(oemol)
    reflected_oemol = reflect_molecule_stereo(oemol)
    return oechem.OEMolToSmiles(oemol) != oechem.OEMolToSmiles(reflected_oemol)  # type: ignore


def check_tetrahedral_coordinates(
    *, atom_coordinates: npt.NDArray[np.float32], neighbour_coordinates: npt.NDArray[np.float32]
) -> None:
    """
    Verify that `atom_coordinates` and `neighbour_coordinates` are the right shape for tetrahedral stereochemistry.

    Shapes should be i.e. (3,) and (3/4, 3) respectively.
    """
    if (
        atom_coordinates.shape != (3,)
        or neighbour_coordinates.shape[1] != 3
        or neighbour_coordinates.shape[0] not in (3, 4)
    ):
        raise ValueError(
            f"Expected shapes: (3,) for atom_coordinates and (3, 3) or (4, 3) for neighbour_coordinates. Got: {atom_coordinates.shape} and {neighbour_coordinates.shape} respectively."
        )


def get_atom_tetrahedral_stereo_from_coordinates(
    *, atom_coordinates: npt.NDArray[np.float32], neighbour_coordinates: npt.NDArray[np.float32]
) -> int:
    """
    Get the OpenEye stereo label for an atom, given its neighbour coordinates.

    This is equivalent to `oechem.GetStereo`, but purely based on coordinates of the central atom
    and its neighbours. It is for use in chirality metrics for DLF.

    Args:
        atom_coordinates: Coordinates of the central atom, as a numpy array of shape `(3, )`
        neighbour_coordinates: Coordinates of the neighbouring atoms, as a numpy array of shape
            `(3, 3)` or `(4, 3)`.

    Returns:
        OpenEye stereo label (either `oechem.OEAtomStereo_Right` or `oechem.OEAtomStereo_Left`)
        based on the order of the neighbouring atoms.

    Raises:
        ValueError: Coordinates are wrong shape.
    """
    check_tetrahedral_coordinates(
        atom_coordinates=atom_coordinates, neighbour_coordinates=neighbour_coordinates
    )
    # Only need three bond directions for parallelepiped
    neighbour_bond_dirs = neighbour_coordinates - atom_coordinates[None, :]
    # Take the (orientated) volume of the parallelepiped from the three vectors leaving the center
    # This is equivalent to the determinant of the 3 vectors stacked, or the scalar triple product
    # Think a.(b x c), b x c describes the plane of two of the arms, and the dot product says which
    # side of the plane it lies
    oriented_volume = np.linalg.det(neighbour_bond_dirs[:3])
    if oriented_volume < 0:
        return oechem.OEAtomStereo_Right  # type: ignore
    else:
        return oechem.OEAtomStereo_Left  # type: ignore


def get_all_tetrahedral_stereo_from_coordinates(
    *,
    coordinates: npt.NDArray[np.float32],
    adjacency_matrix: npt.NDArray[np.int32],
    chiral_mask: npt.NDArray[np.bool_],
) -> list[int]:
    """
    Get a list of chiral labels for the atoms in a molecule, based on coordinates.

    Only atoms in the `chiral_mask` will be given a stereo label. The other atoms will
    be labelled `OEAtomStereo_Undefined`.

    Args:
        coordinates: Coordinates of the atoms.
        adjacency_matrix: Adjacency matrix in the same order as the coordinates.
        chiral_mask: Boolean mask indicating which atoms are chiral centers.

    Returns:
        List of stereo labels, which will be `oechem.OEAtomStereo_Undefined` for nonchiral atoms
        and either `oechem.OEAtomStereo_Right` or `oechem.OEAtomStereo_Left` for chiral atoms.
    """
    chiral_orientations = [oechem.OEAtomStereo_Undefined] * chiral_mask.shape[0]
    for index in chiral_mask.nonzero()[0]:
        atom_coordinates = coordinates[index]
        neighbour_indices = adjacency_matrix[index].nonzero()[0]
        neighbour_coordinates = coordinates[neighbour_indices]
        chiral_orientations[index] = get_atom_tetrahedral_stereo_from_coordinates(
            atom_coordinates=atom_coordinates, neighbour_coordinates=neighbour_coordinates
        )
    return chiral_orientations


def get_atom_tetrahedral_stereo(atom: oechem.OEAtomBase, /) -> int:
    """
    Get the tetrahedral stereo label for an atom, based on the current ordering of atoms.

    Args:
        atom: OpenEye atom.

    Returns:
        One of `oechem.OEAtomStereo_Undefined`, `oechem.OEAtomStereo_Right` or
        `oechem.OEAtomStereo_Left`.
    """
    if not atom.HasStereoSpecified(oechem.OEAtomStereo_Tetrahedral):
        return oechem.OEAtomStereo_Undefined  # type: ignore
    nbrs = list(atom.GetAtoms())
    return atom.GetStereo(nbrs, oechem.OEAtomStereo_Tetrahedral)  # type: ignore


def set_atom_tetrahedral_stereo(atom: oechem.OEAtomBase, /, stereo_label: int) -> None:
    """
    Set the tetrahedral stereo label for an atom, based on the current ordering of atoms.

    Args:
        atom: OpenEye atom.
        stereo_label: OpenEye stereo label, one of `oechem.OEAtomStereo_Undefined`,
            `oechem.OEAtomStereo_Right` or `oechem.OEAtomStereo_Left`.

    Raises:
        ValueError: Failed to set stereo label of atom.
    """
    if not atom.SetStereo(list(atom.GetAtoms()), oechem.OEAtomStereo_Tetrahedral, stereo_label):
        raise ValueError(f"Failed to set stereo for atom {atom} to {stereo_label}")


def get_bond_cistrans_stereo(bond: oechem.OEBondBase, /) -> int:
    """
    Get the cis-trans stereo label for a bond, based on the current ordering of atoms.

    Args:
        bond: OpenEye bond.

    Returns:
        One of `oechem.OEBondStereo_Undefined`, `oechem.OEBondStereo_Cis` or
        `oechem.OEBondStereo_Trans`.
    """
    if not bond.HasStereoSpecified(oechem.OEBondStereo_CisTrans):
        return oechem.OEBondStereo_Undefined  # type: ignore
    start = bond.GetBgn()
    end = bond.GetEnd()
    neigh1 = next(atom for atom in start.GetAtoms() if atom != end)
    neigh2 = next(atom for atom in end.GetAtoms() if atom != start)
    return bond.GetStereo([neigh1, neigh2], oechem.OEBondStereo_CisTrans)  # type: ignore


def set_bond_cistrans_stereo(bond: oechem.OEBondBase, /, stereo_label: int) -> None:
    """
    Set the cis-trans stereo label for a bond, based on the current ordering of atoms.

    Args:
        bond: OpenEye bond.
        stereo_label: OpenEye stereo label, one of `oechem.OEBondStereo_Undefined`,
            `oechem.OEBondStereo_Cis` or `oechem.OEBondStereo_Trans`.

    Raises:
        ValueError: Failed to set stereo label of bond.
    """
    start = bond.GetBgn()
    end = bond.GetEnd()
    neigh1 = next(atom for atom in start.GetAtoms() if atom != end)
    neigh2 = next(atom for atom in end.GetAtoms() if atom != start)
    if not bond.SetStereo([neigh1, neigh2], oechem.OEBondStereo_CisTrans, stereo_label):
        raise ValueError("Failed to set bond stereo")


def get_stereocenters_and_stereobonds(
    oemol: oechem.OEGraphMol, /
) -> tuple[list[oechem.OEAtomBase], list[oechem.OEBondBase]]:
    _, ring_systems = oechem.OEDetermineRingSystems(oemol)

    ring_system_to_possible_stereocenters = defaultdict(list)

    atoms = []
    bonds = []

    for atom in oemol.GetAtoms():
        if atom.IsChiral():
            atoms.append(atom)
        else:
            ring = ring_systems[atom.GetIdx()]
            if ring == 0:
                continue
            nbrs = list(atom.GetAtoms())
            if len(nbrs) < 3:
                continue
            num_ring_nbrs = sum(1 for nbr in nbrs if ring_systems[nbr.GetIdx()] == ring)
            if num_ring_nbrs != 2:
                continue
            ring_symmetries = {
                atom.GetSymmetryClass() for atom in nbrs if ring_systems[atom.GetIdx()] == ring
            }
            if len(ring_symmetries) != 1:
                continue
            ring_system_to_possible_stereocenters[ring].append(atom)
    for bond in oemol.GetBonds():
        if bond.IsChiral():
            bonds.append(bond)

    for possible_stereocenters in ring_system_to_possible_stereocenters.values():
        if len(possible_stereocenters) > 1:
            for stereocenter in possible_stereocenters:
                double_bonds = [bond for bond in stereocenter.GetBonds() if bond.GetOrder() == 2]
                if len(double_bonds) == 1:
                    bonds.append(double_bonds[0])
                else:
                    atoms.append(stereocenter)

    return atoms, bonds


def enumerate_stereoisomers(oemol: oechem.OEMolBase, /) -> Generator[None, None, oechem.OEMolBase]:
    """
    Enumerate all stereoisomers of a molecule.

    Atoms with stereochemistry assigned are not changed. Only atoms which are either unspecified or
    part of enhanced stereochemistry AND or OR groups are flipped.

    Unspecified atoms which can invert are not enumerated. This means tertiary amines that have
    no specified stereochemistry are not enumerated.

    Args:
        oemol: OpenEye molecule to enumerate.

    Yields:
        OpenEye molecule with all stereoisomers.
    """
    oemol = oemol.CreateCopy()

    groups = list(oemol.GetGroups(oechem.OEIsMDLStereoGroup()))

    unspecified_atoms: list[oechem.OEAtomBase] = []
    unspecified_bonds: list[oechem.OEBondBase] = []
    enhanced_stereo_groups: list[list[oechem.OEAtomBase]] = []

    for group in list(groups):
        if group.GetGroupType() == oechem.OEGroupType_MDLAbsStereo:
            oemol.DeleteGroup(group)
            continue
        if group.GetGroupType() not in [
            oechem.OEGroupType_MDLAndStereo,
            oechem.OEGroupType_MDLOrStereo,
        ]:
            continue
        atoms = [
            atom
            for atom in group.GetAtoms()
            if atom.HasStereoSpecified(oechem.OEAtomStereo_Tetrahedral)
        ]
        if len(atoms) >= 2:
            enhanced_stereo_groups.append(atoms)
        else:
            unspecified_atoms.append(atoms[0])
        oemol.DeleteGroup(group)

    stereoatoms, stereobonds = get_stereocenters_and_stereobonds(oemol)

    for atom in stereoatoms:
        if atom.HasStereoSpecified(oechem.OEAtomStereo_Tetrahedral) or can_invert(atom):
            continue
        unspecified_atoms.append(atom)
        set_atom_tetrahedral_stereo(atom, oechem.OEAtomStereo_Left)

    for bond in stereobonds:
        if bond.HasStereoSpecified(oechem.OEBondStereo_CisTrans):
            continue
        unspecified_bonds.append(bond)
        set_bond_cistrans_stereo(bond, oechem.OEBondStereo_Cis)

    all_smiles = set()

    stereo_atom_groups = []
    stereo_bond_groups = []

    for group in enhanced_stereo_groups:
        stereo_atom_groups.append(group)
    for atom in unspecified_atoms:
        stereo_atom_groups.append([atom])

    for bond in unspecified_bonds:
        stereo_bond_groups.append([bond])

    if len(stereo_atom_groups) == 0 and len(stereo_bond_groups) == 0:
        yield oemol
        return

    for stereo_atom_group_list in itertools.product(*([[], group] for group in stereo_atom_groups)):
        for stereo_bond_group_list in itertools.product(
            *([[], group] for group in stereo_bond_groups)
        ):
            atoms = [
                atom for stereo_atom_group in stereo_atom_group_list for atom in stereo_atom_group
            ]
            bonds = [
                bond for stereo_bond_group in stereo_bond_group_list for bond in stereo_bond_group
            ]
            new_oemol = oemol.CreateCopy()
            for atom in atoms:
                reflect_stereocenter(new_oemol.GetAtom(oechem.OEHasAtomIdx(atom.GetIdx())))
            for bond in bonds:
                reflect_stereobond(new_oemol.GetBond(oechem.OEHasBondIdx(bond.GetIdx())))
            smiles = oechem.OEMolToSmiles(new_oemol)
            if smiles in all_smiles:
                continue
            all_smiles.add(smiles)
            yield new_oemol


def can_invert(atom: oechem.OEAtomBase, /) -> bool:
    """
    Can the given atom possibly invert at room temperature?

    Args:
        atom: OpenEye atom to consider.

    Returns:
        True if the atom can be inverted, False otherwise.
    """
    return atom.GetAtomicNum() == 7 and atom.GetDegree() == 3  # type: ignore[no-any-return]
