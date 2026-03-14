from functools import cmp_to_key

import pytest

from chemutils.molecule import oemol_from_smiles
from chemutils.molecule.cip.compare import (
    compare_atomic_mass,
    compare_atomic_numbers,
    compare_recursive,
)
from chemutils.molecule.cip.digraph import AtomNode, HierarchicalDigraph
from chemutils.molecule.cip.set import sorted_partition


class Node:
    def __init__(self, v):
        self.v = v
        self.children = []
        self.parent = None

    def __lt__(self, other):
        return self.v < other.v

    def __repr__(self):
        return f"{self.v}[{','.join(c.v for c in self.children)}]"

    @property
    def ancestors(self):
        if self.parent is None:
            return []
        return [*self.parent.ancestors, self.parent]


def make_graph(branches):
    def create_node(v):
        if isinstance(v, list):
            n = Node(v[0])
            for child in v[1]:
                n_c = create_node(child)
                n.children.append(n_c)
                n_c.parent = n
        else:
            n = Node(v)
        return n

    return [create_node(branch) for branch in branches]


@pytest.mark.parametrize(
    "graph_sorted",
    [
        ("A", "B", "C"),
        (["A", ["B", "C"]], ["A", ["B", "D"]]),
        (["A", ["B", ["C", ["D", "E"]]]], ["A", ["B", ["C", ["D", "F"]]]]),
        (
            ["A", [["B", ["F", "F"]], ["C", ["D", "E"]]]],
            ["A", [["B", ["F", "G"]], ["C", ["D", "E"]]]],
        ),
        (
            [["A", ["F", "F"]], ["B", ["C", ["D", "E"]]]],
            [["A", ["F", "G"]], ["B", ["C", ["D", "E"]]]],
        ),
    ],
)
def test_compare_recursive(graph_sorted):
    graph = make_graph(graph_sorted)
    partition = sorted_partition(reversed(graph), key=cmp_to_key(compare_recursive))
    assert [next(iter(p)) for p in partition] == graph


@pytest.mark.parametrize(
    ["smiles1", "smiles2", "expected", "rule"],
    [
        ("H", "C", -1, compare_atomic_numbers),
        ("O", "C", 1, compare_atomic_numbers),
        ("C", "C", 0, compare_atomic_numbers),
        ("C", "C", 0, compare_atomic_mass),
        ("C", "[14C]", -1, compare_atomic_mass),
        ("[16O]", "O", -1, compare_atomic_mass),
    ],
)
def test_cip_rules(smiles1, smiles2, expected, rule):
    smiles = f"C({smiles1}){smiles2}"
    oemol = oemol_from_smiles(smiles)
    graph = HierarchicalDigraph(root=next(iter(oemol.GetAtoms())))
    node1 = next(
        c for c in graph.root_node.children if isinstance(c, AtomNode) and c.atom.GetIdx() == 1
    )
    node2 = next(
        c for c in graph.root_node.children if isinstance(c, AtomNode) and c.atom.GetIdx() == 2
    )
    assert rule(node1, node2) == expected
