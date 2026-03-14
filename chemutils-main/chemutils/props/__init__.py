from __future__ import annotations

from collections.abc import Mapping
from typing import Any, TypeVar

from openeye import oechem
from rdkit import Chem

from chemutils.chemaxon import cxmol_from_oemol
from chemutils.chemaxon.jvm.chemaxon.struc import Molecule
from chemutils.molecule.representative import representative_oemol_from_smiles
from chemutils.molecule.smiles import oemol_from_smiles
from chemutils.molecule.symmetry import get_symmetry_classes
from chemutils.props.hydrogen_bonds import (
    NumberAcceptorLonePairs,
    NumberDonorHydrogens,
    NumberHydrogenAcceptors,
    NumberHydrogenDonors,
    NumberLipinskiHydrogenAcceptors,
    NumberLipinskiHydrogenDonors,
)
from chemutils.rdmol_from_oemol import rdmol_from_oemol

from .atomic import (
    FractionCarbonSP3,
    MolecularWeight,
    NumberChiralAtoms,
    NumberHeavyAtoms,
    NumberSpiroAtoms,
    SpacialScore,
)
from .blood_brain_barrier import BloodBrainBarrierScore
from .bonds_rings import IsMacrocycle, LargestRingSize, NumberAromaticRings, NumberRotatableBonds
from .cns_mpo import (
    CentralNervousSystemMPO,
    cns_mpo,  # noqa: F401
)
from .druglike import QuantitativeEstimationDrugLikeness
from .filter import PAINS
from .fingerprint import MorganFingerprint
from .identifier import (
    CXSMILES,
    SMILES,
    InChICharge,
    InChIConnections,
    InChIHydrogens,
    InChIKey,
    InChIMolecularFormula,
    InChIStereo,
)
from .logd import LogD
from .logp import LogP, OEXLogP, WildmanCrippenLogP
from .medchem import MEDCHEM_PROPERTIES  # noqa: F401
from .pka import NumberAcidicAtoms, NumberBasicAtoms, PKa, PKaConjugateAcid
from .property import MolecularProperty
from .synthetic_accessibility import SyntheticAccessibilityScore
from .tpsa import TwoDimensionalTPSA
from .uncharged import PercentagedUncharged74, PercentageUncharged6

SymmetryClasses = MolecularProperty("symmetry_classes", oemol_func=get_symmetry_classes)

ALL_PROPERTIES = frozenset(
    [
        BloodBrainBarrierScore,
        NumberRotatableBonds,
        NumberAromaticRings,
        TwoDimensionalTPSA,
        QuantitativeEstimationDrugLikeness,
        OEXLogP,
        WildmanCrippenLogP,
        MolecularWeight,
        NumberHeavyAtoms,
        SpacialScore,
        NumberChiralAtoms,
        NumberSpiroAtoms,
        NumberLipinskiHydrogenAcceptors,
        NumberLipinskiHydrogenDonors,
        MorganFingerprint,
        SyntheticAccessibilityScore,
        IsMacrocycle,
        LargestRingSize,
        SMILES,
        CXSMILES,
        InChIKey,
        InChIHydrogens,
        InChICharge,
        InChIConnections,
        InChIMolecularFormula,
        InChIStereo,
        FractionCarbonSP3,
        LogP,
        LogD,
        PKa,
        PKaConjugateAcid,
        NumberAcidicAtoms,
        NumberBasicAtoms,
        NumberHydrogenDonors,
        NumberHydrogenAcceptors,
        NumberAcceptorLonePairs,
        NumberDonorHydrogens,
        CentralNervousSystemMPO,
        PercentageUncharged6,
        PercentagedUncharged74,
        PAINS,
    ]
)
"""All molecular properties that can be calculated."""

_T = TypeVar("_T")


def calculate_properties_for_smiles(
    smiles: str,
    properties: dict[str, MolecularProperty[_T]] | set[MolecularProperty[_T]],
) -> dict[str, _T]:
    """
    Calculate a set of molecular properties for a given SMILES string.

    Note that no standardisation or neutralisation is applied to the SMILES. You may want to consider
    using `calculate_properties_for_representative_smiles`.

    Args:
        smiles: SMILES string.
        properties: Either a set of MolecularProperty objects (in which case the default key is used for
            each property) or a dictionary mapping names to MolecularProperty values.

    Returns:
        Dictionary, mapping names (either provided or the defaults) to their calculated values.

    """
    oemol = oemol_from_smiles(smiles)
    return calculate_properties_for_oemol(oemol, properties)


def calculate_properties_for_representative_smiles(
    smiles: str,
    properties: dict[str, MolecularProperty[_T]] | set[MolecularProperty[_T]],
) -> dict[str, _T]:
    """
    Calculate a set of molecular properties for a given SMILES string, converting first to a representative form.

    Args:
        smiles: SMILES string.
        properties: Either a set of MolecularProperty objects (in which case the default key is used for
            each property) or a dictionary mapping names to MolecularProperty values.

    Returns:
        Dictionary, mapping names (either provided or the defaults) to their calculated values.

    """
    oemol = representative_oemol_from_smiles(smiles)
    return calculate_properties_for_oemol(oemol, properties)


def get_properties_to_calculate(
    *,
    properties: dict[str, MolecularProperty[_T]] | set[MolecularProperty[_T]],
) -> list[MolecularProperty[Any]]:
    desired_properties: list[MolecularProperty[Any]]

    if isinstance(properties, Mapping):
        desired_properties = list(properties.values())
    else:
        desired_properties = list(properties)

    for property in desired_properties:
        if property.composite_func is not None and property.composite_properties is not None:
            desired_properties.extend(list(property.composite_properties.values()))

    unique_properties = []
    for prop in reversed(desired_properties):
        if prop not in unique_properties:
            unique_properties.append(prop)

    return unique_properties


def calculate_properties_for_oemol(
    oemol: oechem.OEGraphMol,
    properties: dict[str, MolecularProperty[_T]] | set[MolecularProperty[_T]],
) -> dict[str, _T]:
    """
    Calculate a set of molecular properties for a given OpenEye molecule.

    Args:
        oemol: OpenEye molecule.
        properties: Either a set of MolecularProperty objects (in which case the default key is used for
            each property) or a dictionary mapping names to MolecularProperty values.

    Returns:
        Dictionary, mapping names (either provided or the defaults) to their calculated values.

    Raises:
        ValueError: Failed to calculate molecular property.
    """
    calculated_properties = get_properties_to_calculate(properties=properties)

    results: dict[MolecularProperty[Any], Any] = {}
    rdmol: Chem.Mol | None = None
    cxmol: Molecule | None = None

    for property in calculated_properties:
        try:
            if property.composite_func is not None and property.composite_properties is not None:
                results[property] = property.get_composite(results)
            elif property.oemol_func is not None:
                results[property] = property.get_openeye(oemol)
            elif property.rdkit_func is not None:
                rdmol = rdmol if rdmol is not None else rdmol_from_oemol(oemol)
                results[property] = property.get_rdkit(rdmol)
            elif property.cxmol_func is not None:
                cxmol = cxmol if cxmol is not None else cxmol_from_oemol(oemol)
                results[property] = property.get_chemaxon(cxmol)
            else:
                raise ValueError(f"Don't know how to calculate property {property}")
        except Exception as e:
            raise ValueError(f"Failed to calculate property {property}") from e

    if isinstance(properties, Mapping):
        return {name: results[property] for name, property in properties.items()}
    else:
        return {property.id: results[property] for property in properties}
