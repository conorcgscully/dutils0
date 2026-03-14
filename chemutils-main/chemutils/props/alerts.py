from chemutils.molecule.filter.medchem import get_medchem_filter_alerts

from .property import MolecularProperty

NumRedAlerts = MolecularProperty(
    id="num_red_alerts", oemol_func=lambda oemol: get_medchem_filter_alerts(oemol)["num_red_alerts"]
)
NumAmberAlerts = MolecularProperty(
    id="num_amber_alerts",
    oemol_func=lambda oemol: get_medchem_filter_alerts(oemol)["num_amber_alerts"],
)
NumPAINSAlerts = MolecularProperty(
    id="num_pains_alerts",
    oemol_func=lambda oemol: get_medchem_filter_alerts(oemol)["num_pains_alerts"],
)
