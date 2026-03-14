import networkx as nx
from openeye import oechem

from .atom import AtomNode
from .duplicated import DuplicateNode
from .implicit import ImplicitHydrogenNode
from .node import DigraphNode
from .phantom import PhantomNode


class HierarchicalDigraph:
    """
    Implementation of a hierarchical digraph for evaluating CIP rules around a stereocenter.

    See section [P-92.1.4](https://iupac.qmul.ac.uk/BlueBook/P9.html#920104) of the IUPAC Blue Book for
    more information.

    This is a wrapper around a networkx directed graph, with utility methods for creating
    atoms and other nodes relevant to a hierarchical digraph.
    """

    graph: nx.DiGraph
    """Underlying directed graph."""
    root_node: DigraphNode
    """Root node of the graph."""

    def __init__(self, *, root: oechem.OEAtomBase):
        self.graph = nx.DiGraph()
        self.root_node = AtomNode(self, root)
        self.graph.add_node(self.root_node)

    def add_atom(
        self,
        *,
        atom: oechem.OEAtomBase,
        parent: DigraphNode,
    ) -> AtomNode:
        """
        Add a node representing an atom not yet encounted in a branch.

        Args:
            atom: OpenEye atom encounted.
            parent: Existing node that precedes this one.

        Returns:
            Newly added `AtomNode`.
        """
        node = AtomNode(digraph=self, atom=atom)
        self.graph.add_edge(parent, node)
        return node

    def add_duplicated(
        self,
        *,
        duplicated: AtomNode,
        parent: DigraphNode,
    ) -> DuplicateNode:
        """
        Add a node representing a duplicate of an `AtomNode`, either from multiple bonds or cyclic molecules.

        Args:
            duplicated: Previously encounted node.
            parent: Existing node that precedes this one.

        Returns:
            Newly added `DuplicateNode`.
        """
        node = DuplicateNode(digraph=self, duplicated=duplicated)
        self.graph.add_edge(parent, node)
        return node

    def add_implicit_h(self, *, parent: DigraphNode) -> ImplicitHydrogenNode:
        """
        Add a node representing an implicit hydrogen without an underlying `AtomBase`.

        Args:
            parent: Existing node that precedes this one.

        Returns:
            Newly added `ImplicitHydrogenNode`.
        """
        node = ImplicitHydrogenNode(digraph=self)
        self.graph.add_edge(parent, node)
        for _ in range(3):
            self.add_phantom(parent=node)
        return node

    def add_phantom(self, *, parent: DigraphNode) -> PhantomNode:
        """
        Add a phantom node used to bulk out the number of children to be 3.

        Args:
            parent: Existing node that precedes this one.

        Returns:
            Newly added `PhantomNode`.
        """
        node = PhantomNode(digraph=self)
        self.graph.add_edge(parent, node)
        return node
