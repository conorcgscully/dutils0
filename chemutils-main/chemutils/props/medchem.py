from .alerts import NumAmberAlerts, NumPAINSAlerts, NumRedAlerts
from .atomic import (
    FractionCarbonSP3,
    MolecularWeight,
    NumberChiralAtoms,
    NumberHeavyAtoms,
    NumberSpiroAtoms,
    SpacialScore,
)
from .blood_brain_barrier import BloodBrainBarrierScore
from .bonds_rings import IsMacrocycle, NumberAromaticRings, NumberRotatableBonds
from .cns_mpo import CentralNervousSystemMPO
from .druglike import QuantitativeEstimationDrugLikeness
from .hydrogen_bonds import (
    NumberAcceptorLonePairs,
    NumberDonorHydrogens,
    NumberHydrogenAcceptors,
    NumberHydrogenDonors,
    NumberLipinskiHydrogenAcceptors,
    NumberLipinskiHydrogenDonors,
)
from .identifier import CXSMILES, SMILES
from .logd import LogD
from .logp import LogP, OEXLogP, WildmanCrippenLogP
from .pka import PKa, PKaConjugateAcid
from .synthetic_accessibility import SyntheticAccessibilityScore
from .tpsa import TwoDimensionalTPSA
from .uncharged import PercentagedUncharged74, PercentageUncharged6

MEDCHEM_PROPERTIES = {
    "BBB Score": BloodBrainBarrierScore,
    "Num Rotatable Bonds": NumberRotatableBonds,
    "Num Aromatic Rings": NumberAromaticRings,
    "2D TPSA": TwoDimensionalTPSA,
    "QED": QuantitativeEstimationDrugLikeness,
    "Wildman-Crippen ALogP": WildmanCrippenLogP,
    "Molecular Weight": MolecularWeight,
    "Num Heavy Atoms": NumberHeavyAtoms,
    "Spacial Score": SpacialScore,
    "Num Chiral Atoms": NumberChiralAtoms,
    "Num Spiro Atoms": NumberSpiroAtoms,
    "Num Lipinski HBA": NumberLipinskiHydrogenAcceptors,
    "Num Lipinski HBD": NumberLipinskiHydrogenDonors,
    "SA Score": SyntheticAccessibilityScore,
    "Is Macrocycle": IsMacrocycle,
    "SMILES": SMILES,
    "CXSMILES": CXSMILES,
    "Fraction C sp3": FractionCarbonSP3,
    "Chemaxon LogP": LogP,
    "OpenEye XLogP": OEXLogP,
    "Chemaxon LogD": LogD,
    "pKa": PKa,
    "Conjugate Acid pKa": PKaConjugateAcid,
    "Num HBD": NumberHydrogenDonors,
    "Num HBA": NumberHydrogenAcceptors,
    "Num HBA Lone Pairs": NumberAcceptorLonePairs,
    "Num HBD Hydrogens": NumberDonorHydrogens,
    "CNS MPO": CentralNervousSystemMPO,
    "% Uncharged at pH 6": PercentageUncharged6,
    "% Uncharged at pH 7.4": PercentagedUncharged74,
    "Num Red Alerts": NumRedAlerts,
    "Num Amber Alerts": NumAmberAlerts,
    "Num PAINS Alerts": NumPAINSAlerts,
}
"""
Subset of molecular properties relevant to Medicinal Chemistry.

These properties are mapped using human-readable names, and omit
properties such as Morgan Fingerprints which are not relevant.

These are designed to be used with `calculate_properties_for_representative_smiles`,
and to be inserted into software such as Torx & Signals.
"""
