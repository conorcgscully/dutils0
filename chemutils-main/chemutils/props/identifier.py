from chemutils.molecule import (
    representative_cxsmiles_from_oemol,
    representative_smiles_from_oemol,
)
from chemutils.molecule.inchi import get_inchi_breakdown

from .property import MolecularProperty

SMILES = MolecularProperty(
    "smiles",
    oemol_func=lambda oemol: representative_smiles_from_oemol(oemol, already_representative=True),
)

CXSMILES = MolecularProperty(
    "cxsmiles",
    oemol_func=lambda oemol: representative_cxsmiles_from_oemol(oemol, already_representative=True),
)

InChIKey = MolecularProperty(
    "inchi_key", oemol_func=lambda oemol: get_inchi_breakdown(oemol).inchi_key
)
InChIMolecularFormula = MolecularProperty(
    "inchi_molformula", oemol_func=lambda oemol: get_inchi_breakdown(oemol).molformula
)
InChIConnections = MolecularProperty(
    "inchi_connections", oemol_func=lambda oemol: get_inchi_breakdown(oemol).connections
)
InChICharge = MolecularProperty(
    "inchi_charge", oemol_func=lambda oemol: get_inchi_breakdown(oemol).charge
)
InChIStereo = MolecularProperty(
    "inchi_stereo", oemol_func=lambda oemol: get_inchi_breakdown(oemol).stereo
)
InChIHydrogens = MolecularProperty(
    "inchi_hydrogens", oemol_func=lambda oemol: get_inchi_breakdown(oemol).hydrogens
)
