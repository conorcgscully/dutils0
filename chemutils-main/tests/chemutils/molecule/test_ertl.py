import pytest

from chemutils.molecule import oemol_from_smiles
from chemutils.molecule.ertl import get_ertl_functional_groups


@pytest.mark.parametrize(
    ["smiles", "pseudosmiles"],
    [
        # Examples of functional groups detected by Ertl's algorithm
        ### Carbon-only groups
        # Alkene
        ("C=C", {("C=C", (0, 1))}),
        ("C=CCC=C", {("C=C", (0, 1)), ("C=C", (3, 4))}),
        # Alkyne
        ("C#C", {("C#C", (0, 1))}),
        ### Aromatic atoms
        # Aromatic nitrogen
        ("c1ccncc1", {("[Nar]", (3,))}),
        ("c1cn[nH]c1", {("[Nar]", (2,)), ("[Nar]", (3,))}),
        # Aromatic oxygen
        ("c1cocc1", {("[Oar]", (2,))}),
        # Aromatic sulfur
        ("c1cscc1", {("[Sar]", (2,))}),
        ### Alcohols and ethers
        # Alcohol
        ("CO", {("[Cal]O", (1,))}),
        # Phenol
        ("c1ccccc1O", {("[Car]O", (6,))}),
        # Ether
        ("COC", {("[R]O[R]", (1,))}),
        ("COCCOC", {("[R]O[R]", (1,)), ("[R]O[R]", (4,))}),
        # [Acetal
        ("COCOC", {("[R]OCO[R]", (1, 2, 3))}),
        # Carbonyls
        # Ketone
        ("CC(=O)C", {("[R]C(=O)[R]", (1, 2))}),
        # Aldehyde
        ("CC(=O)", {("[R]C=O", (1, 2))}),
        ### Carboxylic acids and derivatives
        # Carboxylic acid
        ("CC(=O)O", {("[R]C(=O)O", (1, 2, 3))}),
        # Ester
        ("CC(=O)OC", {("[R]C(=O)O[R]", (1, 2, 3))}),
        ### Amides and derivatives
        # Amide
        ("CC(=O)N", {("[R]C(=O)N([R])[R]", (1, 2, 3))}),
        ("CC(=O)N(C)C", {("[R]C(=O)N([R])[R]", (1, 2, 3))}),
        ("NC(=O)CC(=O)N", {("[R]C(=O)N([R])[R]", (0, 1, 2)), ("[R]C(=O)N([R])[R]", (4, 5, 6))}),
        # a,b-unsaturated amide
        ("CC=CC(=O)NC", {("[R]N([R])C(=O)C=C", (1, 2, 3, 4, 5))}),
        # Imide
        ("CC(=O)N(C)C(=O)C", {("[R]C(=O)N([R])C(=O)[R]", (1, 2, 3, 5, 6))}),
        # Hydroxamic acid
        ("CC(=O)NO", {("[R]C(=O)N([R])O", (1, 2, 3, 4))}),
        ("CC(=O)N(C)O", {("[R]C(=O)N([R])O", (1, 2, 3, 5))}),
        ### Carbamic acid and derivatives
        # Carbamate ester
        ("CNC(=O)OC", {("[R]N([R])C(=O)O[R]", (1, 2, 3, 4))}),
        # Carbamic acid
        ("CNC(=O)O", {("[R]N([R])C(=O)O", (1, 2, 3, 4))}),
        ### Amines
        # Aliphatic primary amine
        ("CN", {("[Cal]N", (1,))}),
        # Aromatic primary amine
        ("c1ccccc1N", {("[Car]N", (6,))}),
        # Secondary amine
        ("CNC", {("[R]N[R]", (1,))}),
        # Tertiary amine
        ("CN(C)C", {("[R]N([R])[R]", (1,))}),
        ### Imines and derivatives
        # Iminic acid
        ("OC=NC", {("[R]N=CO", (0, 1, 2))}),
        ("CC(=N)O", {("[R]N=CO", (1, 2, 3))}),
        # Amidine
        ("CC(=NC)N(C)C", {("[R]N=CN([R])[R]", (1, 2, 4))}),
        # Guanidine
        ("CNC(=N)N", {("[R]N=C(N([R])[R])N([R])[R]", (1, 2, 3, 4))}),
        ### Other nitrogen groups
        # Nitrile
        ("CC#N", {("C#N", (1, 2))}),
        # Nitro
        ("C[N+](=O)[O-]", {("[R][N+](=O)[O-]", (1, 2, 3))}),
        ### Halides
        # Fluoride
        ("CF", {("[R]F", (1,))}),
        ("C(F)F", {("[R]F", (1,)), ("[R]F", (2,))}),
        # Chloride
        ("CCl", {("[R]Cl", (1,))}),
        ("C(Cl)Cl", {("[R]Cl", (1,)), ("[R]Cl", (2,))}),
        # Bromide
        ("CBr", {("[R]Br", (1,))}),
        ("C(Br)Br", {("[R]Br", (1,)), ("[R]Br", (2,))}),
        # Iodide
        ("CI", {("[R]I", (1,))}),
        ("C(I)I", {("[R]I", (1,)), ("[R]I", (2,))}),
        # Acyl fluoride
        ("CC(=O)F", {("[R]C(=O)F", (1, 2, 3))}),
        # Acyl chloride
        ("CC(=O)Cl", {("[R]C(=O)Cl", (1, 2, 3))}),
        ### Sulfur-containing
        # Thiol
        ("CS", {("[R]S", (1,))}),
        # Sulfide
        ("CSC", {("[R]S[R]", (1,))}),
        # Disulfide
        ("CSSC", {("[R]SS[R]", (1, 2))}),
        # Sulfone
        ("CS(=O)(=O)C", {("[R]S(=O)(=O)[R]", (1, 2, 3))}),
        # Sulfonamide
        ("CS(=O)(=O)N", {("[R]N([R])S(=O)(=O)[R]", (1, 2, 3, 4))}),
        ("CS(=O)(=O)N(C)C", {("[R]N([R])S(=O)(=O)[R]", (1, 2, 3, 4))}),
    ],
)
def test_get_ertl_functional_group_pseudosmiles(smiles, pseudosmiles):
    oemol = oemol_from_smiles(smiles)
    groups = get_ertl_functional_groups(oemol)
    result = {
        (
            next(iter(group.GetRoles())).GetName().removeprefix("ertl:pseudosmiles:"),
            tuple(sorted(atom.GetIdx() for atom in group.GetAtoms())),
        )
        for group in groups
    }
    assert result == pseudosmiles
