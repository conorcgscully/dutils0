from chemutils.chemaxon.pka import get_fraction_microspecies_uncharged

from .property import MolecularProperty

PercentageUncharged6 = MolecularProperty(
    "percent_uncharged_ph_6",
    cxmol_func=lambda oemol: 100 * get_fraction_microspecies_uncharged(oemol, pH=6),
)

PercentagedUncharged74 = MolecularProperty(
    "percent_uncharged_ph_7_4",
    cxmol_func=lambda oemol: 100 * get_fraction_microspecies_uncharged(oemol, pH=7.4),
)
