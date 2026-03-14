from openeye import oechem
from rdkit import Chem
from rdkit.Chem import QED, Crippen, MolSurf, rdMolDescriptors, rdmolops

from chemutils.rdmol_from_oemol import as_rdmol

from .property import MolecularProperty


def get_qed(mol: Chem.Mol | oechem.OEMolBase) -> float:
    """
    Get the Quantitative Estimation of Drug-Likeness (QED) score.

    The QED score is described by Bickerton et al (2012) and compares various molecular
    properties to the distribution of known drug molecules. RDKit makes some minor
    modifications as described at https://www.rdkit.org/docs/source/rdkit.Chem.QED.html.

    Args:
        mol: RDKit molecule.

    Returns:
        QED score, from 0 to 1.

    Raises:
        AtomValenceException: Valencies that are not handled by RDKit, such as hypervalent Nitrogen.
    """
    mol = as_rdmol(mol)
    mol = rdmolops.RemoveAllHs(mol)
    return QED.qed(mol, qedProperties=_properties(mol))  # type: ignore


def _properties(mol: Chem.Mol) -> QED.QEDproperties:
    # Override standard RDKit QED properties, to remove `mol = Chem.RemoveHs(mol)`.
    # We assume the hydrogens have already been stripped, and calling this causes
    # issues with certain structures with hypervalency.
    qedProperties = QED.QEDproperties(
        MW=rdMolDescriptors._CalcMolWt(mol),
        ALOGP=Crippen.MolLogP(mol),
        HBA=sum(
            len(mol.GetSubstructMatches(pattern))
            for pattern in QED.Acceptors
            if mol.HasSubstructMatch(pattern)
        ),
        HBD=rdMolDescriptors.CalcNumHBD(mol),
        PSA=MolSurf.TPSA(mol),
        ROTB=rdMolDescriptors.CalcNumRotatableBonds(
            mol, rdMolDescriptors.NumRotatableBondsOptions.Strict
        ),
        AROM=len(Chem.GetSSSR(Chem.DeleteSubstructs(Chem.Mol(mol), QED.AliphaticRings))),
        ALERTS=sum(1 for alert in QED.StructuralAlerts if mol.HasSubstructMatch(alert)),
    )
    return qedProperties


QuantitativeEstimationDrugLikeness = MolecularProperty(id="qed", rdkit_func=get_qed)
