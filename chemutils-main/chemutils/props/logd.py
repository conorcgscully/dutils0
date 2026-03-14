from chemutils.chemaxon.logd import get_logd

from .property import MolecularProperty

LogD = MolecularProperty("logd", cxmol_func=get_logd)
