from chemutils.filter import openeye_filter

from .property import MolecularProperty

pains_filter = openeye_filter("pains")

PAINS = MolecularProperty("pains_passed", oemol_func=lambda mol: pains_filter(mol).passed)
