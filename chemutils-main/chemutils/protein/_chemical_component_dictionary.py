from functools import cache
from typing import TypedDict, cast

import fsutils as fs
from openeye import oechem
from typing_extensions import NotRequired

from chemutils.molecule import make_hydrogens_implicit, modify_residue, oemol_from_smiles

CHEMICAL_COMPONENT_DICTIONARY_PATH = "s3://charmtx-datalake/rcsb/chemical_components.cbor.gz"


class ChemicalComponent(TypedDict):
    smiles: str
    atom_names: list[str]
    monomer_parent_id: NotRequired[str]


@cache
def _get_chemical_component_dictionary() -> dict[str, ChemicalComponent]:
    return cast(dict[str, ChemicalComponent], fs.read_cbor(CHEMICAL_COMPONENT_DICTIONARY_PATH))


def get_chemical_component(chemical_component_id: str, /) -> oechem.OEGraphMol:
    chemical_component = _get_chemical_component_dictionary().get(chemical_component_id)
    if chemical_component is None:
        raise KeyError(f"Chemical component `{chemical_component_id}` not found.")
    oemol = oemol_from_smiles(chemical_component["smiles"])
    make_hydrogens_implicit(oemol, remove_isotopic_hydrogens=True)
    for atom, name in zip(oemol.GetAtoms(), chemical_component["atom_names"], strict=True):
        atom.SetName(name)
    modify_residue(oemol, residue_name=chemical_component_id)
    return oemol
