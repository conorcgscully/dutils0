from __future__ import annotations

import warnings
from typing import Any

from openeye import oechem
from rdkit import Chem, Geometry
from rdkit.Chem.rdchem import AtomPDBResidueInfo, HybridizationType

ORDER_TO_RDKIT_BONDTYPE = {
    1: Chem.BondType.SINGLE,
    2: Chem.BondType.DOUBLE,
    3: Chem.BondType.TRIPLE,
    4: Chem.BondType.QUADRUPLE,
    5: Chem.BondType.QUINTUPLE,
    6: Chem.BondType.HEXTUPLE,
}

OE_TO_RDKIT_HYBRIDISATION = {
    oechem.OEHybridization_Unknown: HybridizationType.UNSPECIFIED,
    oechem.OEHybridization_sp: HybridizationType.SP,
    oechem.OEHybridization_sp2: HybridizationType.SP2,
    oechem.OEHybridization_sp3: HybridizationType.SP3,
    oechem.OEHybridization_sp3d: HybridizationType.SP3D,
    oechem.OEHybridization_sp3d2: HybridizationType.SP3D2,
}


def clean_rdmol_hydrogens(mol: Chem.Mol, /) -> Chem.Mol:
    """
    Remove hydrogens from an RDKit molecule and recalculate the ring information.

    OpenEye keeps certain hydrogens when importing SMILES where RDKit would delete them.
    When `rdmol_from_oemol` is used, it hence results in additional hydrogens that would
    not be present if using RDKit directly. This causes issues with certain functionality,
    such as Morgan fingerprints.

    Args:
        mol: RDKit molecule with possible hydrogen atoms.

    Returns:
        New RDKit molecule with hydrogens removed.
    """
    if any(atom.GetSymbol() == "H" for atom in mol.GetAtoms()):
        # Morgan fingerprint is affected by explicit hydrogens.
        mol = Chem.RemoveHs(mol, updateExplicitCount=True, sanitize=False)
        Chem.GetSSSR(mol)
    return mol


def as_rdmol(mol: oechem.OEGraphMol | Chem.Mol, /) -> Chem.Mol:
    """
    Ensure a molecule is an RDKit molecule, converting from OpenEye if applicable.

    Args:
        mol: Either an OpenEye or RDKit molecule.

    Returns:
        RDKit molecule.
    """
    if not isinstance(mol, Chem.Mol):
        return rdmol_from_oemol(mol)
    return mol


def rdmol_from_oemol(
    oemol: oechem.OEGraphMol,
    *,
    atom_map: dict[oechem.OEAtomBase, Chem.Atom] | None = None,
    bond_map: dict[oechem.OEBondBase, Chem.Bond] | None = None,
) -> Chem.Mol:
    """
    Convert an OEMol molecule to an RDKit molecule.

    Based on https://gist.github.com/bannanc/810ccc4636b930a4522636baab1965a6, and used under
    the MIT license.

    Args:
        oemol: OpenEye molecule to convert to RDKit.
        atom_map: Optional map that will be populated by mappings from OpenEye atoms to RDKit atoms
            of the resultant atoms.
        bond_map: Optional map that will be populated by mappings from OpenEye bonds to RDKit bonds
            of the resultant bonds.

    Returns:
        Converted RDKit molecule.
    """
    rdmol = Chem.RWMol()

    map_oeatom_idx_to_rdatom_idx = {}
    map_rdatom_idx_to_oeatom_idx = {}

    for oea in oemol.GetAtoms():
        oe_idx = oea.GetIdx()
        rda = Chem.Atom(oea.GetAtomicNum())
        rda.SetFormalCharge(oea.GetFormalCharge())
        rda.SetIsAromatic(oea.IsAromatic())
        rda.SetHybridization(OE_TO_RDKIT_HYBRIDISATION[oea.GetHyb()])
        rda.SetIsotope(oea.GetIsotope())

        rda_idx = rdmol.AddAtom(rda)
        map_oeatom_idx_to_rdatom_idx[oe_idx] = rda_idx
        map_rdatom_idx_to_oeatom_idx[rda_idx] = oe_idx

    map_oebond_idx_to_rdbond_idx = {}

    for oeb in oemol.GetBonds():
        # get neighboring rd atoms
        rd_a1 = map_oeatom_idx_to_rdatom_idx[oeb.GetBgnIdx()]
        rd_a2 = map_oeatom_idx_to_rdatom_idx[oeb.GetEndIdx()]

        rdmol.AddBond(rd_a1, rd_a2)
        rdbond = rdmol.GetBondBetweenAtoms(rd_a1, rd_a2)
        map_oebond_idx_to_rdbond_idx[oeb.GetIdx()] = rdbond.GetIdx()

        # Assign bond type, which is based on order unless it is aromatic
        order = oeb.GetOrder()
        if oeb.IsAromatic():
            rdbond.SetBondType(Chem.BondType.AROMATIC)
            rdbond.SetIsAromatic(True)
        else:
            rdbond.SetBondType(ORDER_TO_RDKIT_BONDTYPE[order])
            rdbond.SetIsAromatic(False)

    # add atom stereochemistry
    for oea in oemol.GetAtoms():
        if oea.HasStereoSpecified(oechem.OEAtomStereo_Tetrahedral):
            rda = rdmol.GetAtomWithIdx(map_oeatom_idx_to_rdatom_idx[oea.GetIdx()])
            # Order neighbours by 'rdkit' reckoning
            # This may not be in the order we created the bonds
            rd_neighbours_idx = [bond.GetOtherAtomIdx(rda.GetIdx()) for bond in rda.GetBonds()]
            oe_neighbours = [
                oemol.GetAtom(oechem.OEHasAtomIdx(map_rdatom_idx_to_oeatom_idx[rd_idx]))
                for rd_idx in rd_neighbours_idx
            ]
            stereo = oea.GetStereo(oe_neighbours, oechem.OEAtomStereo_Tetrahedral)
            if stereo == oechem.OEAtomStereo_Left:
                rda.SetChiralTag(Chem.CHI_TETRAHEDRAL_CCW)
            if stereo == oechem.OEAtomStereo_Right:
                rda.SetChiralTag(Chem.CHI_TETRAHEDRAL_CW)

    # add bond stereochemistry:
    for oebond in oemol.GetBonds():
        if oebond.HasStereoSpecified(oechem.OEBondStereo_CisTrans):
            stereo, oe_idx1, oe_idx2 = get_oemol_cistrans(oebond)
            rd_neighbours_idx = [
                map_oeatom_idx_to_rdatom_idx[oe_idx1],
                map_oeatom_idx_to_rdatom_idx[oe_idx2],
            ]
            oe_neighbours = [
                oemol.GetAtom(oechem.OEHasAtomIdx(map_rdatom_idx_to_oeatom_idx[rd_idx]))
                for rd_idx in rd_neighbours_idx
            ]
            stereo = oebond.GetStereo(oe_neighbours, oechem.OEBondStereo_CisTrans)
            rdbond = rdmol.GetBondWithIdx(map_oebond_idx_to_rdbond_idx[oebond.GetIdx()])
            if stereo == oechem.OEBondStereo_Cis:
                rdbond.SetStereoAtoms(*rd_neighbours_idx)
                rdbond.SetStereo(Chem.rdchem.BondStereo.STEREOCIS)
            if stereo == oechem.OEBondStereo_Trans:
                rdbond.SetStereoAtoms(*rd_neighbours_idx)
                rdbond.SetStereo(Chem.rdchem.BondStereo.STEREOTRANS)

    # if oemol has coordinates (The dimension is non-zero)
    # add those coordinates to the rdmol
    if oemol.GetDimension() > 0:
        conformer = Chem.Conformer()
        oecoords = oemol.GetCoords()
        for oe_idx, rd_idx in map_oeatom_idx_to_rdatom_idx.items():
            (x, y, z) = oecoords[oe_idx]
            conformer.SetAtomPosition(rd_idx, Geometry.Point3D(x, y, z))
        rdmol.AddConformer(conformer)

    # Save the molecule title
    rdmol.SetProp("_Name", oemol.GetTitle())

    rdmol.UpdatePropertyCache(strict=False)
    Chem.GetSSSR(rdmol)
    Chem.AssignStereochemistry(rdmol)
    try:
        Chem.rdCIPLabeler.AssignCIPLabels(rdmol)
    except RuntimeError:
        warnings.warn("CIP labelling failed when converting to RDKit molecule.", stacklevel=1)

    # Handle hydrogens disappearing by making them explicit
    for oea in oemol.GetAtoms():
        oe_idx = oea.GetIdx()
        rd_idx = map_oeatom_idx_to_rdatom_idx[oe_idx]

        rda = rdmol.GetAtomWithIdx(rd_idx)

        oe_implicit_hs = oea.GetImplicitHCount()
        rd_implicit_hs = rda.GetNumImplicitHs()
        if oe_implicit_hs > rd_implicit_hs:
            rda.SetNumExplicitHs(rda.GetNumExplicitHs() + (oe_implicit_hs - rd_implicit_hs))

    # Copy PDB Information
    for atom in oemol.GetAtoms():
        rd_idx = map_oeatom_idx_to_rdatom_idx[atom.GetIdx()]
        oe_data = oechem.OEAtomGetResidue(atom)
        rd_data = AtomPDBResidueInfo(atom.GetName())
        rd_data.SetAltLoc(oe_data.GetAlternateLocation())
        rd_data.SetChainId(oe_data.GetChainID())
        rd_data.SetInsertionCode(oe_data.GetInsertCode())
        rd_data.SetOccupancy(oe_data.GetOccupancy())
        rd_data.SetResidueName(oe_data.GetName())
        rd_data.SetResidueNumber(oe_data.GetResidueNumber())
        rd_data.SetSerialNumber(oe_data.GetSerialNumber())
        rd_data.SetTempFactor(oe_data.GetBFactor())
        rd_data.SetIsHeteroAtom(oe_data.IsHetAtom())
        rdmol.GetAtomWithIdx(rd_idx).SetPDBResidueInfo(rd_data)

    # Copy SD Tags
    for pair in oechem.OEGetSDDataPairs(oemol):
        rdmol.SetProp(pair.GetTag(), pair.GetValue())

    if atom_map is not None:
        for oeatom_idx, rdatom_idx in map_oeatom_idx_to_rdatom_idx.items():
            atom_map[oemol.GetAtom(oechem.OEHasAtomIdx(oeatom_idx))] = rdmol.GetAtomWithIdx(
                rdatom_idx
            )

    if bond_map is not None:
        for oebond_idx, rdbond_idx in map_oebond_idx_to_rdbond_idx.items():
            bond_map[oemol.GetBond(oechem.OEHasBondIdx(oebond_idx))] = rdmol.GetBondWithIdx(
                rdbond_idx
            )

    return rdmol


def get_oemol_cistrans(bond: oechem.OEBondBase) -> tuple[Any, int, int]:
    """
    Get the cis/trans isomerisation of an OpenEye bond, as well as the pair of atom indices to define it.

    Cis/trans isomerisation has to be defined with respect to two other atoms (one bound to each end of
    the double bond). OpenEye doesn't store this, instead it tells you the stereoisomerisation when you
    ask. Hence, to find if a bond is cis/trans, you must iterate over all pairs of neighbours of each
    bond end, and try getting the stereoisomerisation with respect to these two.

    Args:
        bond: OpenEye bond.

    Returns:
        Tuple of the stereoisomisation and the two atom indices used to define it.
    """
    if bond.HasStereoSpecified(oechem.OEBondStereo_CisTrans):
        for atomB in bond.GetBgn().GetAtoms():
            if atomB == bond.GetEnd():
                continue
            for atomE in bond.GetEnd().GetAtoms():
                if atomE == bond.GetBgn():
                    continue
                v = [atomB, atomE]
                stereo = bond.GetStereo(v, oechem.OEBondStereo_CisTrans)
                if stereo != oechem.OEBondStereo_Undefined:
                    return stereo, atomB.GetIdx(), atomE.GetIdx()
    return oechem.OEBondStereo_Undefined, -1, -1
