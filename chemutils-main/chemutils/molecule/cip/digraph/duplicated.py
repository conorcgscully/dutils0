from __future__ import annotations

from collections.abc import Collection
from typing import TYPE_CHECKING

from .atom import AtomNode
from .node import DigraphNode
from .phantom import PhantomNode

if TYPE_CHECKING:
    from .digraph import HierarchicalDigraph


class DuplicateNode(DigraphNode):
    """
    Hierarchical digraph node that represents a duplicate node.

    See section [P-92.1.4.2](https://iupac.qmul.ac.uk/BlueBook/P9.html#92010402) of the IUPAC Blue Book
    for more details on the use of duplicate nodes.
    """

    duplicated: AtomNode
    """Existing atom node that has been duplicated."""

    def __init__(self, *, digraph: HierarchicalDigraph, duplicated: AtomNode):
        super().__init__(digraph=digraph)
        self.duplicated = duplicated

    @property
    def children(self) -> Collection[PhantomNode]:
        if len(super().children) == 0:
            for _ in range(3):
                self.digraph.add_phantom(parent=self)
        return super().children  # type: ignore

    @property
    def name(self) -> str:
        return f"({self.duplicated.name})"

    @property
    def root_distance(self) -> int:
        return self.duplicated.root_distance

    @property
    def atomic_number(self) -> int:
        return self.duplicated.atomic_number

    @property
    def atomic_mass(self) -> float:
        return self.duplicated.atomic_mass
