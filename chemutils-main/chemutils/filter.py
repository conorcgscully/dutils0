"""
Molecular filters which discard molecules based on a set of rules.

Both OpenEye and RDKit provides various filters (and both provide a PAINS filter, though they differ
in implementation).

To create a filter function, use either `openeye_filter` or `rdkit_filter`

```
filter = openeye_filter('BLOCKBUSTER')
```

This filter is a function which accepts a molecule (correspondingly either an `oechem.OEGraphMol` or `Chem.Mol`)
and returns a `FilterResult`.

This result has a `passed` parameter that indicates if the molecule passed the filter, and a `failed_rules`
parameter that lists rules that it failed. Note that for OpenEye, a molecule may fail a filter even though
`failed_rules` is empty. This is because the filters also check properties such as molecular weight.

To create a molecular property that returns True or False if a filter passed, use:

```
filter = openeye_filter('BLOCKBUSTER')
property = MolecularProperty(id='filter_blockbuster', oemol_func=filter.passed)
```

To create a molecular property that returns the set of failed rules, use:

```
filter = rdkit_filter('PAINS')
property = MolecularProperty(id='filter_blockbuster', rdkit_func=filter.failed_rules)
```

"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import NamedTuple

from openeye import oechem, oemolprop
from rdkit import Chem
from rdkit.Chem.rdfiltercatalog import FilterCatalog, FilterCatalogParams


class FilterResult(NamedTuple):
    """Result of a filter on a molecule."""

    passed: bool
    """Has the molecule passed the filter."""
    failed_rules: set[str]
    """List of internal names of rules that have failed for this molecule."""


OE_FILTER_TYPES: dict[str, int] = {}
"""Mapping of string names to OpenEye filter types (integers)."""

for attr in oemolprop.__dict__:
    if attr.startswith("OEFilterType_"):
        OE_FILTER_TYPES[attr.split("OEFilterType_")[1].upper()] = getattr(oemolprop, attr)

CHARMTX_OE_FILTER_NAME_TO_FILE: dict[str, Path] = {
    # sarah_skerratt_virtual_screening_filter.txt
    # Taken from https://www.notion.so/charmtx/Medchem-Filters-1172f7d9efb9478ca25ddae2f542f2d8
    "CHARM_20221101": (Path(__file__).parent / "CHARM_20221101.txt"),
}

OE_FILTEROUTPUT_NONRULES = {
    "SMILES",
    "Name",
    "atom count",
    "carbons",
    "heteroatoms",
    "heteroatom to carbon ratio",
    "chiral centers",
    "hydrogen-bond acceptors",
    "hydrogen-bond donors",
    "lipinski h-bond acceptors",
    "lipinski h-bond donors",
    "molecular weight",
    "molecular weight halide fraction",
    "formal charge count",
    "sum of formal charges",
    "Non-ring size",
    "Unbranched chain",
    "total fcn ct",
    "Veber",
    "Lipinski violations",
    "ABS",
    "Pharmacopia",
    "Aggregator",
    "Predicted Aggregator",
    "Filter",
    "2d PSA",
    "Solubility",
    "XLogP",
    "rigid bonds",
    "rotatable bonds",
    "maximum size of ring system",
    "number of ring systems",
}
"""OpenEye filter output contains various properties as well as rules.

This set is the standard OpenEye headers which may appear in a filter output.
"""


def charmtx_openeye_filter(name: str, /) -> Callable[[oechem.OEGraphMol], FilterResult]:
    """
    Create an OpenEye molecular filter from a CHARMTx filter file.

    Args:
        name: Name of a pre defined filter (i.e. 'CHARM_20221101')

    Returns:
        A function which maps an OpenEye molecule to a FilterResult.

    Raises:
        ValueError: The filter is unrecognized.
    """
    if name in CHARMTX_OE_FILTER_NAME_TO_FILE:
        ifs = oechem.oeifstream()
        filter_text_file = CHARMTX_OE_FILTER_NAME_TO_FILE[name]
        if not ifs.open(str(filter_text_file)):
            raise ValueError(
                f"Could not open oechem filter file {filter_text_file} that corresponds to {name}"
            )
        oe_filter = oemolprop.OEFilter(ifs)
        return openeye_filter(oe_filter)
    else:
        raise ValueError(f"Unrecognized filter {name}")


def openeye_filter(
    filter: str | int | oemolprop.OEFilter,
) -> Callable[[oechem.OEGraphMol], FilterResult]:
    """
    Create an OpenEye molecular filter.

    Args:
        filter: Either a string name of a filter (i.e. 'BLOCKBUSTER'), a constant (i.e.
            `oemolprop.OEFilterType_BlockBuster`) or a already created OpenEye filter.

    Returns:
        A function which maps an OpenEye molecule to a FilterResult.

    Raises:
        ValueError: The filter is unrecognized.
    """
    if isinstance(filter, str):
        try:
            filter = OE_FILTER_TYPES[filter.upper()]
        except KeyError as e:
            raise ValueError(f"Unrecognized filter type: {filter}") from e

    if isinstance(filter, int):
        filter = oemolprop.OEFilter(filter)

    def filter_mol(mol: oechem.OEGraphMol) -> FilterResult:
        mol = mol.CreateCopy()
        ostr = oechem.oeosstream()
        filter.SetTable(ostr, False)
        filter.SetFlagTableFailures(True)
        filter(mol)
        [headers, fields] = [
            line.split("\t") for line in ostr.str().decode().rstrip("\n").split("\n")
        ]
        results = {header.strip(): value for header, value in zip(headers, fields, strict=True)}
        return FilterResult(
            passed=results["Filter"] == "Pass",
            failed_rules={name for name, value in results.items() if value.endswith("*")},
        )

    filter_mol.passed = lambda mol: filter_mol(mol).passed  # type: ignore
    filter_mol.failed_rules = lambda mol: filter_mol(mol).failed_rules  # type: ignore

    return filter_mol


def rdkit_filter(
    filter: str | FilterCatalogParams.FilterCatalogs,
) -> Callable[[Chem.Mol], FilterResult]:
    """
    Create an RDKit molecular filter.

    Args:
        filter: Either a string name of a filter (i.e. 'PAINS') or a constant (i.e.
            `FilterCatalogs.PAINS`)

    Returns:
        A function which maps an RDKit molecule to a FilterResult.

    Raises:
        ValueError: The filter is unrecognized.
    """
    if isinstance(filter, str):
        try:
            filter = getattr(FilterCatalogParams.FilterCatalogs, filter.upper())
        except AttributeError as e:
            raise ValueError(f"Unrecognized filter {filter}") from e

    filter = FilterCatalog(filter)

    def filter_mol(mol: Chem.Mol) -> FilterResult:
        # RDkit filter consists of many SMARTS matches for substructures that should be filtered out
        # A success occurs when no matches are found
        matches = filter.GetMatches(mol)
        failed_rules = {match.GetDescription() for match in matches}
        return FilterResult(passed=failed_rules == set(), failed_rules=failed_rules)

    filter_mol.passed = lambda mol: filter_mol(mol).passed  # type: ignore
    filter_mol.failed_rules = lambda mol: filter_mol(mol).failed_rules  # type: ignore

    return filter_mol
