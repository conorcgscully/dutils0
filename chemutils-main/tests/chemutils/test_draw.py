import io

import pytest
from openeye import oechem
from PIL import Image
from syrupy.extensions.image import PNGImageSnapshotExtension, SVGImageSnapshotExtension

from chemutils.draw import Highlight
from chemutils.draw import draw_molecule as draw_molecule_actual
from chemutils.molecule import make_hydrogens_implicit, oemol_from_smiles
from chemutils.molecule.ertl import get_ertl_functional_groups
from tests.chemutils.props.drug_smiles import DRUG_SMILES

EXTENSION = {"svg": SVGImageSnapshotExtension, "png": PNGImageSnapshotExtension}


def clean_png(png: bytes) -> bytes:
    original = Image.open(io.BytesIO(png))
    stripped = Image.new(original.mode, original.size)
    stripped.putdata(list(original.getdata()))

    out = io.BytesIO()
    stripped.save(out, format="PNG")
    return out.getvalue()


# Wrapper around draw_molecule, to strip out metadata from PNG
# It makes image comparison difficult, due to RDKit serialising a whole
# molecule (with exact RDKit version) into the PNG metadata
def draw_molecule(*args, format, **kargs):
    result = draw_molecule_actual(*args, format=format, **kargs)

    if format != "png":
        return result

    return clean_png(result)


@pytest.mark.parametrize("format", ["png", "svg"])
def test_draw(snapshot, format):
    mol = oemol_from_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    assert draw_molecule(
        mol,
        format=format,
    ) == snapshot(extension_class=EXTENSION[format])


@pytest.mark.parametrize("format", ["png", "svg"])
def test_label(snapshot, format):
    mol = oemol_from_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    assert draw_molecule(
        mol,
        label="Caffeine",
        format=format,
    ) == snapshot(extension_class=EXTENSION[format])


@pytest.mark.parametrize("format", ["png", "svg"])
@pytest.mark.parametrize("shape", [(100, 100), (700, 300), (300, 700), (1024, 1024)])
def test_image_shape(snapshot, format, shape):
    mol = oemol_from_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    assert draw_molecule(
        mol,
        width=shape[0],
        height=shape[1],
        format=format,
    ) == snapshot(extension_class=EXTENSION[format])


@pytest.mark.parametrize("format", ["png", "svg"])
@pytest.mark.parametrize("background_color", ["red", "#ccffdd", None])
def test_background_color(snapshot, format, background_color):
    mol = oemol_from_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    assert draw_molecule(mol, format=format, background_color=background_color) == snapshot(
        extension_class=EXTENSION[format]
    )


@pytest.mark.parametrize("format", ["png", "svg"])
def test_grayscale(snapshot, format):
    mol = oemol_from_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    assert draw_molecule(mol, format=format, grayscale=True) == snapshot(
        extension_class=EXTENSION[format]
    )


@pytest.mark.parametrize("format", ["png", "svg"])
def test_atom_index_labels(snapshot, format):
    mol = oemol_from_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    assert draw_molecule(mol, format=format, atom_labels=lambda atom: atom.GetIdx()) == snapshot(
        extension_class=EXTENSION[format]
    )


@pytest.mark.parametrize("format", ["png", "svg"])
def test_highlight_atoms_single_atom_index(snapshot, format):
    mol = oemol_from_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")

    assert draw_molecule(mol, format=format, highlight_atoms=6) == snapshot(
        extension_class=EXTENSION[format]
    )


@pytest.mark.parametrize("format", ["png", "svg"])
def test_highlight_atom_index(snapshot, format):
    mol = oemol_from_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")

    assert draw_molecule(mol, format=format, highlight_atoms=6) == snapshot(
        extension_class=EXTENSION[format]
    )
    assert draw_molecule(mol, format=format, highlight=Highlight(atoms=6, color="red")) == snapshot(
        extension_class=EXTENSION[format]
    )


@pytest.mark.parametrize("format", ["png", "svg"])
def test_highlight_single_atom(snapshot, format):
    mol = oemol_from_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    atoms = list(mol.GetAtoms())

    assert draw_molecule(mol, format=format, highlight_atoms=atoms[2]) == snapshot(
        extension_class=EXTENSION[format]
    )
    assert draw_molecule(
        mol, format=format, highlight=Highlight(atoms=atoms[2], color="blue")
    ) == snapshot(extension_class=EXTENSION[format])


@pytest.mark.parametrize("format", ["png", "svg"])
def test_highlight_atom_index_list(snapshot, format):
    mol = oemol_from_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")

    assert draw_molecule(
        mol,
        format=format,
        highlight_atoms=[3, 7, 11],
    ) == snapshot(extension_class=EXTENSION[format])
    assert draw_molecule(
        mol,
        format=format,
        highlight=Highlight(atoms=[3, 7, 11], color="green"),
    ) == snapshot(extension_class=EXTENSION[format])


@pytest.mark.parametrize("format", ["png", "svg"])
def test_highlight_atom_list(snapshot, format):
    mol = oemol_from_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    atoms = list(mol.GetAtoms())

    assert draw_molecule(
        mol,
        format=format,
        highlight_atoms=[atoms[2], atoms[5], atoms[10]],
    ) == snapshot(extension_class=EXTENSION[format])
    assert draw_molecule(
        mol,
        format=format,
        highlight=Highlight(atoms=[atoms[2], atoms[5], atoms[10]], color="grey"),
    ) == snapshot(extension_class=EXTENSION[format])


@pytest.mark.parametrize("format", ["png", "svg"])
def test_highlight_atom_iter(snapshot, format):
    mol = oemol_from_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    atoms = list(mol.GetAtoms(oechem.OEHasHvyDegree(2)))

    assert draw_molecule(
        mol,
        format=format,
        highlight_atoms=atoms,
    ) == snapshot(extension_class=EXTENSION[format])
    assert draw_molecule(
        mol,
        format=format,
        highlight=Highlight(atoms=atoms, color="grey"),
    ) == snapshot(extension_class=EXTENSION[format])


@pytest.mark.parametrize("format", ["png", "svg"])
def test_highlight_atom_predicate(snapshot, format):
    mol = oemol_from_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")

    assert draw_molecule(
        mol,
        format=format,
        highlight_atoms=oechem.OEIsAromaticAtom(),
    ) == snapshot(extension_class=EXTENSION[format])
    assert draw_molecule(
        mol,
        format=format,
        highlight=Highlight(atoms=oechem.OEIsAromaticAtom(), color="grey"),
    ) == snapshot(extension_class=EXTENSION[format])


@pytest.mark.parametrize("format", ["png", "svg"])
def test_highlight_atom_callable(snapshot, format):
    mol = oemol_from_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")

    assert draw_molecule(
        mol,
        format=format,
        highlight_atoms=lambda atom: atom.GetDegree() == 3,
    ) == snapshot(extension_class=EXTENSION[format])
    assert draw_molecule(
        mol,
        format=format,
        highlight=Highlight(atoms=lambda atom: atom.GetDegree() == 3, color="grey"),
    ) == snapshot(extension_class=EXTENSION[format])


@pytest.mark.parametrize("format", ["png", "svg"])
def test_highlight_atom_bond_set_atoms(snapshot, format):
    mol = oemol_from_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    atoms = list(mol.GetAtoms())
    atom_bond_set = oechem.OEAtomBondSet()
    atom_bond_set.AddAtom(atoms[2])
    atom_bond_set.AddAtom(atoms[6])

    assert draw_molecule(
        mol,
        format=format,
        highlight_atoms=atom_bond_set,
    ) == snapshot(extension_class=EXTENSION[format])
    assert draw_molecule(
        mol,
        format=format,
        highlight=Highlight(atoms=atom_bond_set, color="grey"),
    ) == snapshot(extension_class=EXTENSION[format])


@pytest.mark.parametrize("format", ["png", "svg"])
def test_highlight_match_atoms(snapshot, format):
    mol = oemol_from_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    atoms = list(mol.GetAtoms())
    match = oechem.OEMatch()
    match.AddPair(atoms[2], atoms[6])

    assert draw_molecule(
        mol,
        format=format,
        highlight_atoms=match,
    ) == snapshot(extension_class=EXTENSION[format])
    assert draw_molecule(
        mol,
        format=format,
        highlight=Highlight(atoms=match, color="grey"),
    ) == snapshot(extension_class=EXTENSION[format])


@pytest.mark.parametrize("format", ["png", "svg"])
def test_highlight_single_bond_index(snapshot, format):
    mol = oemol_from_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")

    assert draw_molecule(mol, format=format, highlight_bonds=3) == snapshot(
        extension_class=EXTENSION[format]
    )
    assert draw_molecule(
        mol, format=format, highlight=Highlight(bonds=3, color="cyan")
    ) == snapshot(extension_class=EXTENSION[format])


@pytest.mark.parametrize("format", ["png", "svg"])
def test_highlight_single_bond(snapshot, format):
    mol = oemol_from_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    bonds = list(mol.GetBonds())

    assert draw_molecule(mol, format=format, highlight_bonds=bonds[2]) == snapshot(
        extension_class=EXTENSION[format]
    )
    assert draw_molecule(
        mol, format=format, highlight=Highlight(bonds=bonds[2], color="blue")
    ) == snapshot(extension_class=EXTENSION[format])


@pytest.mark.parametrize("format", ["png", "svg"])
def test_highlight_bond_index_list(snapshot, format):
    mol = oemol_from_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")

    assert draw_molecule(
        mol,
        format=format,
        highlight_bonds=[2, 5, 10],
    ) == snapshot(extension_class=EXTENSION[format])
    assert draw_molecule(
        mol,
        format=format,
        highlight=Highlight(bonds=[2, 5, 10], color="pink"),
    ) == snapshot(extension_class=EXTENSION[format])


@pytest.mark.parametrize("format", ["png", "svg"])
def test_highlight_bond_list(snapshot, format):
    mol = oemol_from_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    bonds = list(mol.GetBonds())

    assert draw_molecule(
        mol,
        format=format,
        highlight_bonds=[bonds[2], bonds[5], bonds[10]],
    ) == snapshot(extension_class=EXTENSION[format])
    assert draw_molecule(
        mol,
        format=format,
        highlight=Highlight(bonds=[bonds[2], bonds[5], bonds[10]], color="orange"),
    ) == snapshot(extension_class=EXTENSION[format])


@pytest.mark.parametrize("format", ["png", "svg"])
def test_highlight_bond_iter(snapshot, format):
    mol = oemol_from_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    bonds = list(mol.GetBonds(oechem.OEIsAromaticBond()))

    assert draw_molecule(
        mol,
        format=format,
        highlight_bonds=bonds,
    ) == snapshot(extension_class=EXTENSION[format])
    assert draw_molecule(
        mol,
        format=format,
        highlight=Highlight(bonds=bonds, color="red"),
    ) == snapshot(extension_class=EXTENSION[format])


@pytest.mark.parametrize("format", ["png", "svg"])
def test_highlight_bond_predicate(snapshot, format):
    mol = oemol_from_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")

    assert draw_molecule(
        mol,
        format=format,
        highlight_bonds=oechem.OEIsAromaticBond(),
    ) == snapshot(extension_class=EXTENSION[format])
    assert draw_molecule(
        mol,
        format=format,
        highlight=Highlight(bonds=oechem.OEIsAromaticBond(), color="blue"),
    ) == snapshot(extension_class=EXTENSION[format])


@pytest.mark.parametrize("format", ["png", "svg"])
def test_highlight_bond_callable(snapshot, format):
    mol = oemol_from_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")

    assert draw_molecule(
        mol,
        format=format,
        highlight_bonds=lambda bond: bond.GetOrder() == 2,
    ) == snapshot(extension_class=EXTENSION[format])
    assert draw_molecule(
        mol,
        format=format,
        highlight=Highlight(bonds=lambda bond: bond.GetOrder() == 2, color="green"),
    ) == snapshot(extension_class=EXTENSION[format])


@pytest.mark.parametrize("format", ["png", "svg"])
def test_highlight_bond_atom_bond_set(snapshot, format):
    mol = oemol_from_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    bonds = list(mol.GetBonds())
    atom_bond_set = oechem.OEAtomBondSet()
    atom_bond_set.AddBond(bonds[2])
    atom_bond_set.AddBond(bonds[6])

    assert draw_molecule(
        mol,
        format=format,
        highlight_bonds=atom_bond_set,
    ) == snapshot(extension_class=EXTENSION[format])
    assert draw_molecule(
        mol,
        format=format,
        highlight=Highlight(bonds=atom_bond_set, color="grey"),
    ) == snapshot(extension_class=EXTENSION[format])


@pytest.mark.parametrize("format", ["png", "svg"])
def test_highlight_bond_match(snapshot, format):
    mol = oemol_from_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    bonds = list(mol.GetBonds())
    match = oechem.OEMatch()
    match.AddPair(bonds[2], bonds[6])

    assert draw_molecule(
        mol,
        format=format,
        highlight_bonds=match,
    ) == snapshot(extension_class=EXTENSION[format])
    assert draw_molecule(
        mol,
        format=format,
        highlight=Highlight(bonds=match, color="grey"),
    ) == snapshot(extension_class=EXTENSION[format])


@pytest.mark.parametrize("format", ["png", "svg"])
def test_highlight_atom_bond_list(snapshot, format):
    mol = oemol_from_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    atoms = list(mol.GetAtoms())
    bonds = list(mol.GetBonds())

    assert draw_molecule(mol, format=format, highlight=[atoms[2], bonds[2]]) == snapshot(
        extension_class=EXTENSION[format]
    )
    assert draw_molecule(
        mol,
        format=format,
        highlight=[
            Highlight(atoms_and_bonds=[atoms[2], bonds[2]], color="red"),
        ],
    ) == snapshot(extension_class=EXTENSION[format])


@pytest.mark.parametrize("format", ["png", "svg"])
def test_highlight_atom_and_bonds_atom_bond_set(snapshot, format):
    mol = oemol_from_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    atoms = list(mol.GetAtoms())
    bonds = list(mol.GetBonds())
    atom_bond_set = oechem.OEAtomBondSet()
    atom_bond_set.AddAtom(atoms[2])
    atom_bond_set.AddBond(bonds[2])

    assert draw_molecule(mol, format=format, highlight=atom_bond_set) == snapshot(
        extension_class=EXTENSION[format]
    )
    assert draw_molecule(
        mol,
        format=format,
        highlight=[
            Highlight(atoms_and_bonds=atom_bond_set, color="red"),
        ],
    ) == snapshot(extension_class=EXTENSION[format])


@pytest.mark.parametrize("format", ["png", "svg"])
def test_highlight_atom_and_bonds_match(snapshot, format):
    mol = oemol_from_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    atoms = list(mol.GetAtoms())
    bonds = list(mol.GetBonds())
    match = oechem.OEMatch()
    match.AddPair(atoms[2], atoms[7])
    match.AddPair(bonds[2], bonds[5])

    assert draw_molecule(mol, format=format, highlight=match) == snapshot(
        extension_class=EXTENSION[format]
    )
    assert draw_molecule(
        mol,
        format=format,
        highlight=[
            Highlight(atoms_and_bonds=match, color="red"),
        ],
    ) == snapshot(extension_class=EXTENSION[format])


@pytest.mark.parametrize("format", ["png", "svg"])
def test_highlight_multiple(snapshot, format):
    mol = oemol_from_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    atoms = list(mol.GetAtoms())

    assert draw_molecule(
        mol,
        format=format,
        highlight=[
            Highlight(atoms=oechem.OEIsAromaticAtom(), color="grey"),
            Highlight(atoms=[atoms[2], atoms[3], atoms[4]], color="pink"),
        ],
    ) == snapshot(extension_class=EXTENSION[format])


@pytest.mark.parametrize("format", ["png", "svg"])
def test_highlight_multiple_bonds(snapshot, format):
    mol = oemol_from_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    bonds = list(mol.GetBonds())

    assert draw_molecule(
        mol,
        format=format,
        highlight=[
            Highlight(bonds=oechem.OEIsAromaticBond(), color="grey"),
            Highlight(bonds=[bonds[2], bonds[3], bonds[4]], color="pink"),
        ],
    ) == snapshot(extension_class=EXTENSION[format])


@pytest.mark.parametrize("format", ["png", "svg"])
def test_highlight_multiple_connected(snapshot, format):
    mol = oemol_from_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    atoms = list(mol.GetAtoms())

    assert draw_molecule(
        mol,
        format=format,
        highlight=[
            Highlight(atoms=oechem.OEIsAromaticAtom(), color="grey", connected=True),
            Highlight(atoms=[atoms[2], atoms[3], atoms[4]], color="pink", connected=True),
        ],
    ) == snapshot(extension_class=EXTENSION[format])


@pytest.mark.parametrize("format", ["png", "svg"])
def test_highlight_multiple_overwrite(snapshot, format):
    mol = oemol_from_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    atoms = list(mol.GetAtoms())

    assert draw_molecule(
        mol,
        format=format,
        highlight=[
            Highlight(atoms=oechem.OEIsAromaticAtom(), color="grey", connected=True),
            Highlight(
                atoms=[atoms[2], atoms[3], atoms[4]], color="pink", connected=True, overwrite=True
            ),
        ],
    ) == snapshot(extension_class=EXTENSION[format])


@pytest.mark.parametrize("format", ["png", "svg"])
def test_highlight_multiple_bonds_overwrite(snapshot, format):
    mol = oemol_from_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    bonds = list(mol.GetBonds())

    assert draw_molecule(
        mol,
        format=format,
        highlight=[
            Highlight(bonds=oechem.OEIsAromaticBond(), color="grey", connected=True),
            Highlight(
                bonds=[bonds[2], bonds[3], bonds[4]], color="pink", connected=True, overwrite=True
            ),
        ],
    ) == snapshot(extension_class=EXTENSION[format])


@pytest.mark.parametrize("format", ["png", "svg"])
def test_condense_abbreviations(snapshot, format):
    mol = oemol_from_smiles("O=C1N[C@@H]2[C@@H](SC[C@@H]2N1)CCCCC(=O)O")
    make_hydrogens_implicit(mol, remove_isotopic_hydrogens=True)
    assert draw_molecule(mol, format=format, condense_abbreviations=True) == snapshot(
        extension_class=EXTENSION[format]
    )


@pytest.mark.parametrize("format", ["png", "svg"])
def test_stereo_label_abs(snapshot, format):
    mol = oemol_from_smiles("[C@H](F)(Cl)Br |a:0|")

    assert draw_molecule(mol, format=format, condense_abbreviations=True) == snapshot(
        extension_class=EXTENSION[format]
    )


@pytest.mark.parametrize("format", ["png", "svg"])
def test_stereo_label_and(snapshot, format):
    mol = oemol_from_smiles("[C@H](F)(Cl)Br |&1:0|")

    assert draw_molecule(mol, format=format, condense_abbreviations=True) == snapshot(
        extension_class=EXTENSION[format]
    )


@pytest.mark.parametrize("format", ["png", "svg"])
def test_stereo_label_or(snapshot, format):
    mol = oemol_from_smiles("[C@H](F)(Cl)Br |o1:0|")

    assert draw_molecule(mol, format=format, condense_abbreviations=True) == snapshot(
        extension_class=EXTENSION[format]
    )


@pytest.mark.parametrize("format", ["png", "svg"])
def test_stereo_label_and_multiple(snapshot, format):
    mol = oemol_from_smiles("[C@H](F)(Cl)[C@H](Br)I |&1:0,&2:3|")

    assert draw_molecule(mol, format=format, condense_abbreviations=True) == snapshot(
        extension_class=EXTENSION[format]
    )


@pytest.mark.parametrize("format", ["png", "svg"])
def test_stereo_label_and_linked(snapshot, format):
    mol = oemol_from_smiles("[C@H](F)(Cl)[C@H](Br)I |&1:0,3|")

    assert draw_molecule(mol, format=format, condense_abbreviations=True) == snapshot(
        extension_class=EXTENSION[format]
    )


@pytest.mark.parametrize("format", ["png", "svg"])
def test_stereo_label_or_multiple(snapshot, format):
    mol = oemol_from_smiles("[C@H](F)(Cl)[C@H](Br)I |o1:0,o2:3|")

    assert draw_molecule(mol, format=format, condense_abbreviations=True) == snapshot(
        extension_class=EXTENSION[format]
    )


@pytest.mark.parametrize("format", ["png", "svg"])
def test_stereo_label_or_linked(snapshot, format):
    mol = oemol_from_smiles("[C@H](F)(Cl)[C@H](Br)I |o1:0,3|")

    assert draw_molecule(mol, format=format, condense_abbreviations=True) == snapshot(
        extension_class=EXTENSION[format]
    )


@pytest.mark.parametrize("format", ["png", "svg"])
def test_stereo_label_mixed(snapshot, format):
    mol = oemol_from_smiles("[C@H](F)(Cl)[C@H](Br)I |&1:0,o1:3|")

    assert draw_molecule(mol, format=format, condense_abbreviations=True) == snapshot(
        extension_class=EXTENSION[format]
    )


def test_oemol_repr_png(snapshot):
    mol = oemol_from_smiles("[C@H](F)(Cl)Br")

    assert clean_png(mol._repr_png_()) == snapshot(extension_class=EXTENSION["png"])


def test_oemol_repr_png_stereo(snapshot):
    mol = oemol_from_smiles("[C@H](F)(Cl)Br |&1:0|")

    assert clean_png(mol._repr_png_()) == snapshot(extension_class=EXTENSION["png"])


def test_oeatom_repr_png(snapshot):
    mol = oemol_from_smiles("[C@H](F)(Cl)Br")
    atom = list(mol.GetAtoms())[1]
    assert clean_png(atom._repr_png_()) == snapshot(extension_class=EXTENSION["png"])


def test_oeatomiter_repr_png(snapshot):
    mol = oemol_from_smiles("[C@H](F)(Cl)Br")
    atoms = mol.GetAtoms()
    assert clean_png(atoms._repr_png_()) == snapshot(extension_class=EXTENSION["png"])


def test_oebond_repr_png(snapshot):
    mol = oemol_from_smiles("[C@H](F)(Cl)Br")
    bond = list(mol.GetBonds())[1]
    assert clean_png(bond._repr_png_()) == snapshot(extension_class=EXTENSION["png"])


def test_oebonditer_repr_png(snapshot):
    mol = oemol_from_smiles("[C@H](F)(Cl)Br")
    bonds = mol.GetBonds()
    assert clean_png(bonds._repr_png_()) == snapshot(extension_class=EXTENSION["png"])


@pytest.mark.parametrize(
    "smiles", ["c1ccccc1O", "CC(=O)Nc1ccc(O)cc1", r"Cc1cnc(s1)NC(=O)C\3=C(/O)c2ccccc2S(=O)(=O)N/3C"]
)
def test_oeatombondsetvector_repr_png(snapshot, smiles):
    mol = oemol_from_smiles(smiles)
    groups = get_ertl_functional_groups(mol)
    assert clean_png(groups._repr_png_()) == snapshot(extension_class=EXTENSION["png"])


@pytest.mark.parametrize(
    "smiles", [pytest.param(smi, id=name) for name, smi in DRUG_SMILES.items()]
)
def test_drug(snapshot, smiles):
    mol = oemol_from_smiles(smiles)
    make_hydrogens_implicit(mol, remove_isotopic_hydrogens=True)
    assert draw_molecule(mol, format="svg") == snapshot(extension_class=EXTENSION["svg"])


@pytest.mark.parametrize("aromaticity", ["kekulize", "dashed", "circle"])
@pytest.mark.parametrize("format", ["png", "svg"])
def test_aromaticity(snapshot, aromaticity, format):
    mol = oemol_from_smiles("c1c2c(nc[nH]2)ncn1")
    assert draw_molecule(mol, format=format, aromaticity=aromaticity) == snapshot(
        extension_class=EXTENSION[format]
    )


def test_many_elements(snapshot):
    mol = oemol_from_smiles(
        "C1(F)C(Cl)=C([SeH])C(NC(C2C=C(P(O)(O)=O)C(Br)=C2C(I)B(O)O)=S)=CC=1C(=O)O"
    )
    assert draw_molecule(mol, format="svg") == snapshot(extension_class=EXTENSION["svg"])


@pytest.mark.parametrize("format", ["png", "svg"])
def test_unknown_bond_stereo(snapshot, format):
    mol = oemol_from_smiles("C(F)(Cl)=C(Br)I")
    assert draw_molecule(mol, format=format) == snapshot(extension_class=EXTENSION[format])
