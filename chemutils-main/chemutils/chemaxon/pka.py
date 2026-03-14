from collections import defaultdict
from collections.abc import MutableMapping
from typing import Any, Literal

from chemutils.molecule import oemol_from_smiles

from .ion_class import get_ion_class
from .jvm.chemaxon.marvin.calculations import pKaPlugin
from .molecule import ConvertibleToCXMol, as_cxmol
from .tautomer import canonical_smiles_from_cxmol


def get_microspecies_distribution(
    mol: ConvertibleToCXMol,
    /,
    *,
    pH: float = 7.4,
    temperature: float | None = None,
    consider_tautomomerization: bool = False,
) -> dict[str, float]:
    """
    Get the distribution of microspecies at a given pH.

    Args:
        mol: Either a SMILES or an OpenEye molecule.
        pH: pH to calculate microspecies distribution at.
        temperature: Temperature to calculate at, in Kelvin.
        consider_tautomomerization: Should tautomers and resonance be considered?

    Returns:
        Mapping of microspecies (either as SMILES or OpenEye molecules, depending on input)
        to the calculated fraction that this microspecies exists.
    """
    cxmol = as_cxmol(mol)

    plugin = pKaPlugin()
    if temperature:
        plugin.setTemperature(float(temperature))
    if consider_tautomomerization:
        plugin.setConsiderTautomerization(True)

    plugin.setpH(float(pH))
    plugin.setMolecule(cxmol)

    assert plugin.isMsCalc()
    plugin.run()

    result: MutableMapping[str, float] = defaultdict(lambda: 0)
    smiles_to_cxmol = {}

    for i in range(plugin.getMsCount()):
        # Deduplicate by SMILES
        cxmol = plugin.getMsMolecule(i)
        smiles = canonical_smiles_from_cxmol(cxmol)
        smiles_to_cxmol[smiles] = cxmol
        result[smiles] += plugin.getSingleMsDistribution(i) / 100

    return dict(result)


def get_fraction_microspecies_uncharged(
    mol: ConvertibleToCXMol, /, *, pH: float = 7.4, temperature: float | None = None
) -> float:
    """
    Get the overall fraction of microspecies at a given pH that have no charged atoms.

    This excludes both ions and zwitterions with a overall neutral charge.

    Args:
        mol: Either a SMILES or an OpenEye molecule.
        pH: pH to calculate microspecies distribution at.
        temperature: Temperature to calculate at, in Kelvin.

    Returns:
        Fraction of microspecies with uncharged atoms.
    """
    distr: dict[str, float] = get_microspecies_distribution(mol, pH=pH, temperature=temperature)

    total = 0.0

    for smiles, fraction in distr.items():
        if get_ion_class(oemol_from_smiles(smiles)) == "uncharged":
            total += fraction

    return total


def _get_pkas(
    mol: ConvertibleToCXMol,
    /,
    *,
    temperature: float | None = None,
    mode: Literal["micro", "macro"],
    type: Any,
) -> list[float]:
    cxmol = as_cxmol(mol)

    plugin = pKaPlugin()

    if temperature:
        plugin.setTemperature(float(temperature))

    plugin.setMolecule(cxmol)
    if mode == "micro":
        plugin.setMicropKaCalc(True)

    plugin.run()

    pkas = plugin.getMacropKaValues(type)

    return list(pkas) if pkas else []


def get_macro_pka(mol: ConvertibleToCXMol, /, *, temperature: float | None = None) -> float | None:
    """
    Calculate the macro-pKa of a molecule.

    The provided molecule is converted to a canonical form, so this functions will return the same
    result for anionic, cationic, neutral and zwitterionic versions of the same molecule.

    Args:
        mol: Either a SMILES or an OpenEye molecule.
        temperature: Temperature to calculate pKa at.

    Returns:
        pKa of the molecule.
    """
    pkas = _get_pkas(mol, temperature=temperature, mode="macro", type=pKaPlugin.ACIDIC)
    return pkas[0] if len(pkas) > 0 else None


def get_macro_pka_conjugate_acid(
    mol: ConvertibleToCXMol, /, *, temperature: float | None = None
) -> float | None:
    """
    Calculate the macro-pKa of the conjugate acid of a molecule.

    The provided molecule is converted to a canonical form, so this functions will return the same
    result for anionic, cationic, neutral and zwitterionic versions of the same molecule.

    Args:
        mol: Either a SMILES or an OpenEye molecule.
        temperature: Temperature to calculate pKa at.

    Returns:
        pKa of the molecule.
    """
    pkas = _get_pkas(mol, temperature=temperature, mode="macro", type=pKaPlugin.BASIC)
    return pkas[0] if len(pkas) > 0 else None


def get_num_acidic_atoms(mol: ConvertibleToCXMol, /, *, temperature: float | None = None) -> int:
    """
    Get the number of acidic atoms in a molecule.

    An acidic atom is defined as any atom with a micro-pKa of 7 or less, as calculated by ChemAxon.

    Args:
        mol: Either a SMILES or an OpenEye molecule.
        temperature: Temperature to calculate pKa at.

    Returns:
        Number of acidic atoms in the molecule.
    """
    return len(
        [
            pka
            for pka in _get_pkas(mol, temperature=temperature, mode="micro", type=pKaPlugin.ACIDIC)
            if pka <= 7
        ]
    )


def get_num_basic_atoms(mol: ConvertibleToCXMol, /, *, temperature: float | None = None) -> int:
    """
    Get the number of basic atoms in a molecule.

    A basic atom is defined as any atom with a micro-pKb of 7 or less (or equivalently, a micro-pKa of
    7 or more for the conjugate acid) as calculated by ChemAxon.

    Args:
        mol: Either a SMILES or an OpenEye molecule.
        temperature: Temperature to calculate pKa at.

    Returns:
        Number of basic atoms in the molecule.
    """
    return len(
        [
            pka
            for pka in _get_pkas(mol, temperature=temperature, mode="micro", type=pKaPlugin.BASIC)
            if pka >= 7
        ]
    )
