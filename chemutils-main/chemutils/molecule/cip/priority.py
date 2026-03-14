from __future__ import annotations

from functools import cmp_to_key, partial

from openeye import oechem

from .compare import (
    compare_atomic_mass,
    compare_atomic_numbers,
    compare_bond_stereo,
    compare_duplicate_nodes,
    compare_recursive,
)
from .digraph import AtomNode, HierarchicalDigraph
from .set import PartitionFunction, refine_partitions, sorted_partition


def neighbours_ordered_by_cip_priority(atom: oechem.OEAtomBase) -> list[oechem.OEAtomBase]:
    """
    Get the list of neighbouring atoms ranked by their CIP priority, from lowest to highest.

    Args:
        atom: Atom which is a tetrahedral stereocenter.

    Returns:
        List of neighbouring atoms ordered by their CIP priority.

    Raises:
        ValueError: Failed to uniquely order neighbours.
    """
    digraph = HierarchicalDigraph(root=atom)
    # Get all nodes connected to the root which represent actual OpenEye atoms.
    neighbour_nodes = [node for node in digraph.root_node.children if isinstance(node, AtomNode)]

    groups: list[list[AtomNode]] = [list(neighbour_nodes)]

    for comparison_rule in [
        compare_atomic_numbers,
        compare_duplicate_nodes,
        compare_atomic_mass,
        compare_bond_stereo,
    ]:
        # Key function that compares nodes recursively by rule
        key = cmp_to_key(partial(compare_recursive, comparison_function=comparison_rule))
        partition: PartitionFunction[AtomNode] = partial(sorted_partition, key=key)

        # Split nodes by rule
        groups = refine_partitions(groups, partition_function=partition)

        # Rule has successfully separated all branches - each group will contain one and exactly one node
        if len(groups) == len(neighbour_nodes):
            return [next(iter(group)).atom for group in groups]

    raise ValueError("Failed to uniquely order neighbours using CIP")
