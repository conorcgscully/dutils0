import pytest

from chemutils.molecule import oemol_from_smiles
from chemutils.molecule.inchi import InChI, get_inchi_breakdown


@pytest.mark.parametrize(
    ["smiles", "expected"],
    [
        # Charged ion
        (
            "[Fe+2]",
            InChI(
                inchi_key="CWYNVVGOOAEACU-UHFFFAOYSA-N",
                molformula="Fe",
                charge="/q+2",
            ),
        ),
        # Unspecified tetrahedral stereochemistry
        (
            "[C](F)(Cl)(Br)I",
            InChI(
                inchi_key="XEGUVFFZWHRVAV-UHFFFAOYSA-N",
                molformula="CBrClFI",
                connections="/c2-1(3,4)5",
            ),
        ),
        # Clockwise tetrahedral stereochemistry
        (
            "[C@](F)(Cl)(Br)I",
            InChI(
                inchi_key="XEGUVFFZWHRVAV-PVQJCKRUSA-N",
                molformula="CBrClFI",
                connections="/c2-1(3,4)5",
                stereo="/t1-/m1/s1",
            ),
        ),
        # Anticlockwise tetrahedral stereochemistry
        (
            "[C@@](F)(Cl)(Br)I",
            InChI(
                inchi_key="XEGUVFFZWHRVAV-SFOWXEAESA-N",
                molformula="CBrClFI",
                connections="/c2-1(3,4)5",
                stereo="/t1-/m0/s1",
            ),
        ),
        # Unspecified double bond stereochemistry
        (
            "FC(Cl)=C(Br)I",
            InChI(
                inchi_key="RMZBCOYELGKOKR-UHFFFAOYSA-N",
                molformula="C2BrClFI",
                connections="/c3-1(6)2(4)5",
            ),
        ),
        # Trans double bond stereochemistry
        (
            "F/C(Cl)=C/(Br)I",
            InChI(
                inchi_key="RMZBCOYELGKOKR-UPHRSURJSA-N",
                molformula="C2BrClFI",
                connections="/c3-1(6)2(4)5",
                stereo="/b2-1-",
            ),
        ),
        # Cis double bond stereochemistry
        (
            "F/C(Cl)=C\\(Br)I",
            InChI(
                inchi_key="RMZBCOYELGKOKR-OWOJBTEDSA-N",
                molformula="C2BrClFI",
                connections="/c3-1(6)2(4)5",
                stereo="/b2-1+",
            ),
        ),
        # Pair of structural isomers (same molformula, different connections/hydrogens)
        (
            "CCCC",
            InChI(
                inchi_key="IJDNQMDRQITEOD-UHFFFAOYSA-N",
                molformula="C4H10",
                connections="/c1-3-4-2",
                hydrogens="/h3-4H2,1-2H3",
            ),
        ),
        (
            "CC(C)C",
            InChI(
                inchi_key="NNPPMTNAJDCUHE-UHFFFAOYSA-N",
                molformula="C4H10",
                connections="/c1-4(2)3",
                hydrogens="/h4H,1-3H3",
            ),
        ),
        # Pair of tautomers (same molformula/connections, different hydrogens)
        (
            "C=CO",
            InChI(
                inchi_key="IMROMDMJAWUWLK-UHFFFAOYSA-N",
                molformula="C2H4O",
                connections="/c1-2-3",
                hydrogens="/h2-3H,1H2",
            ),
        ),
        (
            "CC=O",
            InChI(
                inchi_key="IKHGUXGNUITLKF-UHFFFAOYSA-N",
                molformula="C2H4O",
                connections="/c1-2-3",
                hydrogens="/h2H,1H3",
            ),
        ),
        # Pair of microspecies (same molformula/connections/hydrogens, different charges)
        (
            "CC(=O)O",
            InChI(
                inchi_key="QTBSBXVTEAMEQO-UHFFFAOYSA-N",
                molformula="C2H4O2",
                connections="/c1-2(3)4",
                hydrogens="/h1H3,(H,3,4)",
            ),
        ),
        (
            "CC(=O)[O-]",
            InChI(
                inchi_key="QTBSBXVTEAMEQO-UHFFFAOYSA-M",
                molformula="C2H4O2",
                connections="/c1-2(3)4",
                hydrogens="/h1H3,(H,3,4)",
                charge="/p-1",
            ),
        ),
        # Isotopic hydrogen
        (
            "CC[2H]",
            InChI(
                inchi_key="OTMSDBZUPAUEDD-MICDWDOJSA-N",
                molformula="C2H6",
                connections="/c1-2",
                hydrogens="/h1-2H3",
                isotope="/i1D",
            ),
        ),
        # Isotopic carbon
        (
            "C[14CH2]C",
            InChI(
                inchi_key="ATUOYWHBWRKTHZ-YZRHJBSPSA-N",
                molformula="C3H8",
                connections="/c1-3-2",
                hydrogens="/h3H2,1-2H3",
                isotope="/i3+2",
            ),
        ),
    ],
)
def test_get_inchi_breakdown(smiles, expected):
    oemol = oemol_from_smiles(smiles)
    assert get_inchi_breakdown(oemol) == expected
