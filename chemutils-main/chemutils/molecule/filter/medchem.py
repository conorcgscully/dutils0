from functools import cache
from pathlib import Path
from typing import Any, Literal, TypedDict

import fsutils as fs
from openeye import oechem

from chemutils.molecule import oemol_from_smiles
from chemutils.molecule.graph import iterate_substructure_matches


@cache
def _get_medchem_filter(version: str, /) -> Any:
    return fs.read_yaml(Path(__file__).parent / f"medchem_filter_{version}.yaml")


class MedChemFilterAlert(TypedDict):
    category: Literal["red", "amber", "pains"]
    rule: str
    matches: list[oechem.OEAtomBondSet]


class MedChemFilterResult(TypedDict):
    num_red_alerts: int
    num_amber_alerts: int
    num_pains_alerts: int
    alerts: list[MedChemFilterAlert]


def get_medchem_filter_alerts(
    molecule: oechem.OEMolBase, /, *, version: str = "2025-06-09"
) -> MedChemFilterResult:
    molecule = oemol_from_smiles(molecule) if isinstance(molecule, str) else molecule
    medchem_filter = _get_medchem_filter(version)
    alerts = []
    for filter_group in medchem_filter.values():
        matches: set[tuple[frozenset[oechem.OEAtomBase], frozenset[oechem.OEBondBase]]] = set()
        for _, smarts in filter_group["filter_smarts"].items():
            for match in iterate_substructure_matches(pattern=smarts, target=molecule):
                atoms = frozenset(atom for atom in match.GetTargetAtoms())
                bonds = frozenset(bond for bond in match.GetTargetBonds())
                # merge into existing matches
                new_matches = set()
                for existing_atoms, existing_bonds in matches:
                    if existing_atoms & atoms:
                        atoms = existing_atoms | atoms
                        bonds = existing_bonds | bonds
                    else:
                        new_matches.add((existing_atoms, existing_bonds))
                new_matches.add((atoms, bonds))
                matches = new_matches
        if len(matches) > 0:
            oematches = []
            for atoms, bonds in matches:
                oematch = oechem.OEAtomBondSet()
                for atom in atoms:
                    oematch.AddAtom(atom)
                for bond in bonds:
                    oematch.AddBond(bond)
                oematches.append(oematch)
            alerts.append(
                MedChemFilterAlert(
                    rule=filter_group["filter_group"],
                    category=filter_group["category"],
                    matches=oematches,
                )
            )
    return MedChemFilterResult(
        num_red_alerts=len([alert for alert in alerts if alert["category"] == "red"]),
        num_amber_alerts=len([alert for alert in alerts if alert["category"] == "amber"]),
        num_pains_alerts=len([alert for alert in alerts if alert["category"] == "pains"]),
        alerts=alerts,
    )
