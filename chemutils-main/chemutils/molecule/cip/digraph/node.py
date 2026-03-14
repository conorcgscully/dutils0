from __future__ import annotations

from abc import ABC, abstractproperty
from collections.abc import Collection
from typing import TYPE_CHECKING

import networkx as nx

if TYPE_CHECKING:
    from .digraph import HierarchicalDigraph


class DigraphNode(ABC):
    """Base class for a node of a hierarchical digraph."""

    digraph: HierarchicalDigraph
    """Hierarchical digraph which this node belongs to."""

    def __init__(self, *, digraph: HierarchicalDigraph):
        self.digraph = digraph

    @property
    def children(self) -> Collection[DigraphNode]:
        """List of nodes which are direct children of this node."""
        return list(self.digraph.graph.successors(self))

    @property
    def parent(self) -> DigraphNode | None:
        """Node which is the parent of this node."""
        parents = list(self.digraph.graph.predecessors(self))
        if len(parents) == 0:
            return None
        return parents[0]  # type: ignore

    @property
    def ancestors(self) -> Collection[DigraphNode]:
        """Collection of nodes which precede this node in the directed graph."""
        return list(nx.algorithms.dag.ancestors(self.digraph.graph, self))

    @abstractproperty
    @property
    def name(self) -> str:
        """Small user-friendly name of the node."""
        raise NotImplementedError

    @abstractproperty
    @property
    def atomic_number(self) -> int:
        """Atomic number of the node, used to evaluate Rule 1a."""
        raise NotImplementedError

    @abstractproperty
    @property
    def atomic_mass(self) -> float:
        """Isotopic mass of the node, used to evaluate Rule 2."""
        raise NotImplementedError

    @property
    def root_distance(self) -> int:
        """Distance of the node from the root, used to evaluate Rule 1b."""
        return len(self.ancestors)

    def __str__(self) -> str:
        return self.name
