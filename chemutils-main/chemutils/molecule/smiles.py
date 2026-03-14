import re

from openeye import oechem

from .standardise import standardise_oemol
from .stereo import set_bond_cistrans_stereo
from .write import write_molecule_str

CXSMILES_REGEX = re.compile(r"^(?P<smiles>.+) \|(?P<tags>.+)\|$")

ENHANCED_STEREO_REGEX = re.compile(r"(&\d+|o\d+|a):\d+(,\d+)*")


def clean_cxsmiles_tags(tags: str, /) -> str:
    return ",".join(match.group() for match in ENHANCED_STEREO_REGEX.finditer(tags))


def oemol_from_smiles(smiles: str, /) -> oechem.OEGraphMol:
    """
    Convert a SMILES or CXSMILES string to an OpenEye molecule.

    No standardisation is applied (such as removing hydrogens or charges), this can be
    done using a subsequent call to `standardise_oemol`, or in one step using the
    function `std_oemol_from_smiles`.

    Args:
        smiles: SMILES or CXSMILES string.

    Returns:
        OpenEye molecule parsed from SMILES or CXSMILES string.

    Raises:
        ValueError: Failed to parse SMILES or CXSMILES string.
    """
    # If OpenEye sees CXSMILES tags it can't understand, it fails to parse anything
    is_cxsmiles = "|" in smiles
    has_bond_stereo = "\\" in smiles or "/" in smiles
    if is_cxsmiles and (match := CXSMILES_REGEX.match(smiles)):
        smiles = f"{match.group('smiles')} |{clean_cxsmiles_tags(match.group('tags'))}|"
    mol = oechem.OEGraphMol()
    opts = oechem.OEParseSmilesOptions(cxsmiles=is_cxsmiles)
    if not oechem.OEParseSmiles(mol, smiles, opts):
        raise ValueError(f"Failed to parse SMILES string {smiles}")
    mol.SetTitle("")  # Remove title (which gets set to the CXSMILES tags)
    oechem.OEFindRingAtomsAndBonds(mol)
    oechem.OEAssignAromaticFlags(mol, oechem.OEAroModel_OpenEye)
    oechem.OEPerceiveChiral(mol)
    if has_bond_stereo:
        # OEMolToSmiles remove cis/trans from double bonds in rings of size 6 or less.
        for bond in mol.GetBonds():
            if not bond.HasStereoSpecified(oechem.OEBondStereo_CisTrans):
                continue
            if (
                oechem.OEBondIsInRingSize(bond, 6)
                or oechem.OEBondIsInRingSize(bond, 5)
                or oechem.OEBondIsInRingSize(bond, 4)
            ):
                set_bond_cistrans_stereo(bond, oechem.OEBondStereo_Undefined)
    return mol


def std_oemol_from_smiles(smiles: str, /) -> oechem.OEGraphMol:
    """
    Convert a SMILES string to a standardised OpenEye molecule.

    The various standardisation procedures include the removal of charges and hydrogens,
    atom reordering and discarding all but the largest component are performed, as described
    by `standardise_oemol`.

    Args:
        smiles: SMILES string.

    Returns:
        OpenEye molecule parsed from SMILES string, after standardisation.

    Raises:
        ValueError: Failed to parse SMILES string.
    """
    mol = oemol_from_smiles(smiles)
    standardise_oemol(mol)
    return mol


def cxsmiles_from_oemol(oemol: oechem.OEMolBase, /) -> str:
    """
    Convert an OpenEye molecule to its CXSMILES representation.

    Args:
        oemol: OpenEye molecule.

    Returns:
        CXSMILES string.
    """
    return write_molecule_str(oemol, format="cxsmiles").strip()


def smiles_from_oemol(
    oemol: oechem.OEMolBase, /, *, atom_order: list[oechem.OEAtomBase] | None = None
) -> str:
    """
    Convert an OpenEye molecule to its SMILES representation.

    Args:
        oemol: OpenEye molecule.
        atom_order: Optional list to populate with the OpenEye atoms in the order they appear in the SMILES string.

    Returns:
        SMILES string.

    Raises:
        ValueError: Invalid argument for `atom_order`.
    """
    if atom_order is not None:
        if len(atom_order) != 0:
            raise ValueError("`atom_order` argument must be an empty list.")
        atom_order.extend(
            oechem.OEGetSmiStringOrder(
                oemol, oechem.OESMILESFlag_Canonical | oechem.OESMILESFlag_ISOMERIC
            )
        )
    return oechem.OEMolToSmiles(oemol)  # type: ignore


def canonicalize_smiles(smiles: str, /) -> str:
    """Canonicalize a SMILES string to a canonical isomeric SMILES."""
    return cxsmiles_from_oemol(oemol_from_smiles(smiles))
