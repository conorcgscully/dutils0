import pytest

from chemutils.molecule import oemol_from_smiles
from chemutils.molecule.cip.digraph import (
    HierarchicalDigraph,
    ImplicitHydrogenNode,
    PhantomNode,
)
from chemutils.molecule.cip.digraph.atom import make_subscript


@pytest.mark.parametrize(
    ["smiles", "expected"],
    [
        ("C", "C₀[H03,H03,H03,H03]"),
        ("CC", "C₀[C₁[H03,H03,H03],H03,H03,H03]"),
        ("C=C", "C₀[C₁[H03,H03]01,H03,H03]"),
        ("CC=C", "C₀[C₁[C₂[(C₁)03,H03,H03],(C₂)03,H03],H03,H03,H03]"),
        ("CC#C", "C₀[C₁[C₂[(C₁)03,(C₁)03,H03],(C₂)03,(C₂)03],H03,H03,H03]"),
        ("C1CC1", "C₀[C₂[C₁[(C₀)03,H03,H03],H03,H03],C₁[C₂[(C₀)03,H03,H03],H03,H03],H03,H03]"),
    ],
)
def test_smiles(smiles, expected):
    def make_node_name(node):
        name = f"{node.name}"
        children = node.children
        num_phantom = len([c for c in children if isinstance(c, PhantomNode)])
        children = [c for c in children if not isinstance(c, PhantomNode)]
        if len(children) > 0:
            name += "[" + ",".join(make_node_name(n) for n in children) + "]"
        if num_phantom > 0:
            name += "0" + str(num_phantom)
        return name

    oemol = oemol_from_smiles(smiles)
    digraph = HierarchicalDigraph(root=next(iter(oemol.GetAtoms())))
    assert make_node_name(digraph.root_node) == expected


@pytest.fixture
def digraph():
    oemol = oemol_from_smiles("C")
    return HierarchicalDigraph(root=next(iter(oemol.GetAtoms())))


EXAMPLE_ATOMS = [
    ("C", 6, 12.0107, "C₀"),
    ("[14C]", 6, 14.0032, "C₀"),
    ("O", 8, 15.9994, "O₀"),
    ("[H]", 1, 1.00794, "H"),
]


@pytest.mark.parametrize(["smiles", "atomic_number", "atomic_mass", "name"], EXAMPLE_ATOMS)
def test_atom_node(smiles, atomic_number, atomic_mass, name, digraph):
    oemol = oemol_from_smiles(smiles)
    atom_node = digraph.add_atom(atom=next(iter(oemol.GetAtoms())), parent=digraph.root_node)
    assert atom_node.atomic_number == atomic_number
    assert atom_node.atomic_mass == pytest.approx(atomic_mass, abs=0.01)
    assert atom_node.name == name


@pytest.mark.parametrize(["smiles", "atomic_number", "atomic_mass", "name"], EXAMPLE_ATOMS)
def test_duplicated_node(smiles, atomic_number, atomic_mass, name, digraph):
    oemol = oemol_from_smiles(smiles)
    atom_node = digraph.add_atom(atom=next(iter(oemol.GetAtoms())), parent=digraph.root_node)
    duplicated_node = digraph.add_duplicated(duplicated=atom_node, parent=atom_node)
    assert duplicated_node.atomic_number == atomic_number
    assert duplicated_node.atomic_mass == pytest.approx(atomic_mass, abs=0.01)
    assert duplicated_node.name == "(" + name + ")"


def test_implicit_hydrogen_node(digraph):
    node = ImplicitHydrogenNode(digraph=digraph)
    assert node.atomic_number == 1
    assert node.atomic_mass == pytest.approx(1.008, abs=0.01)
    assert node.name == "H"


def test_phantom_node(digraph):
    node = PhantomNode(digraph=digraph)
    assert node.atomic_number == 0
    assert node.atomic_mass == 0
    assert node.name == "0"


@pytest.mark.parametrize(
    ["value", "subscript"], [(0, "₀"), (12, "₁₂"), (543, "₅₄₃"), (6879, "₆₈₇₉")]
)
def test_make_subscript(value, subscript):
    assert make_subscript(value) == subscript


@pytest.mark.parametrize("value", [-1, -55])
def test_make_subscript_error(value):
    with pytest.raises(ValueError):
        _ = make_subscript(value)
