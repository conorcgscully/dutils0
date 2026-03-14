from chemutils.chemaxon.pka import (
    get_macro_pka,
    get_macro_pka_conjugate_acid,
    get_num_acidic_atoms,
    get_num_basic_atoms,
)

from .property import MolecularProperty

PKa = MolecularProperty("pka", cxmol_func=get_macro_pka)
PKaConjugateAcid = MolecularProperty("pka_conjugate_acid", cxmol_func=get_macro_pka_conjugate_acid)
NumberAcidicAtoms = MolecularProperty("num_acidic_atoms", cxmol_func=get_num_acidic_atoms)
NumberBasicAtoms = MolecularProperty("num_basic_atoms", cxmol_func=get_num_basic_atoms)
