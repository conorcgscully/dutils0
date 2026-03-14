from collections.abc import Callable, Collection
from enum import Enum
from functools import cmp_to_key, lru_cache, partial
from typing import Any

from openeye import oechem

from .digraph import AtomNode, DigraphNode, PhantomNode
from .set import refine_partitions, sorted_partition


class Comparison(int, Enum):
    ALessThanB = -1
    AMoreThanB = 1
    AEqualToB = 0


NodeComparisonFunction = Callable[[DigraphNode, DigraphNode], Comparison]


def default_comparison(a: Any, b: Any) -> Comparison:
    """
    Default comparison function that uses rich comparison.

    Args:
        a: The first object to compare.
        b: The second object to compare.

    Returns:
        -1 if `a` is less than `b`, 1 if `a` is more than `b` and 0 otherwise.
    """
    if a < b:
        return Comparison.ALessThanB
    if a > b:
        return Comparison.AMoreThanB
    return Comparison.AEqualToB


@lru_cache(128)
def compare_atomic_numbers(a: DigraphNode, b: DigraphNode, /) -> Comparison:
    """
    Compare two nodes by atomic number, ordering from low to high.

    Args:
        a: The first node to compare.
        b: The second node to compare.

    Returns:
        -1 if `a` is less than `b`, 1 if `a` is more than `b` and 0 otherwise.
    """
    if a.atomic_number < b.atomic_number:
        return Comparison.ALessThanB
    if a.atomic_number > b.atomic_number:
        return Comparison.AMoreThanB
    return Comparison.AEqualToB


@lru_cache(128)
def compare_duplicate_nodes(a: DigraphNode, b: DigraphNode, /) -> Comparison:
    """
    Compare two nodes using Subrule 1b of the Cahn-Ingold-Prelog rules.

    This rule is expressed as (IUPAC Blue Book, P-92.2.2 Sequence Subrule 1b):
        'a nearer duplicate atom node takes precedence over a further away duplicate atom node'

    This comparison functions only compares two nodes which both represent duplicated atom nodes.
    In this case, it ranks them (low to high) from the node with the most ancestors to the fewest.

    Args:
        a: The first node to compare.
        b: The second node to compare.

    Returns:
        -1 if `a` is less than `b`, 1 if `a` is more than `b` and 0 otherwise.
    """
    if a.root_distance > b.root_distance:
        return Comparison.ALessThanB
    if a.root_distance < b.root_distance:
        return Comparison.AMoreThanB
    return Comparison.AEqualToB


@lru_cache(128)
def compare_atomic_mass(a: DigraphNode, b: DigraphNode, /) -> Comparison:
    """
    Compare two nodes by atomic mass (accounting for isotopes).

    Args:
        a: The first node to compare.
        b: The second node to compare.

    Returns:
        -1 if `a` is less than `b`, 1 if `a` is more than `b` and 0 otherwise.
    """
    if a.atomic_mass < b.atomic_mass:
        return Comparison.ALessThanB
    if a.atomic_mass > b.atomic_mass:
        return Comparison.AMoreThanB
    return Comparison.AEqualToB


@lru_cache(128)
def compare_bond_stereo(a: DigraphNode, b: DigraphNode) -> Comparison:
    """
    Compare two nodes by double bond stereochemistry, with Z > E > non-stereo.

    Args:
        a: The first node to compare.
        b: The second node to compare.

    Returns:
        -1 if `a` is less than `b`, 1 if `a` is more than `b` and 0 otherwise.
    """
    if isinstance(a, AtomNode):
        bond_a = a.bond_to_parent
        a_stereo = oechem.OEPerceiveCIPStereo(bond_a.GetParent(), bond_a)
    else:
        a_stereo = oechem.OECIPBondStereo_NotStereo

    if isinstance(b, AtomNode):
        bond_b = b.bond_to_parent
        b_stereo = oechem.OEPerceiveCIPStereo(bond_b.GetParent(), bond_b)
    else:
        b_stereo = oechem.OECIPBondStereo_NotStereo

    if a_stereo == oechem.OECIPBondStereo_E and b_stereo == oechem.OECIPBondStereo_Z:
        return Comparison.ALessThanB
    if a_stereo == oechem.OECIPBondStereo_Z and b_stereo == oechem.OECIPBondStereo_E:
        return Comparison.AMoreThanB
    if (
        a_stereo in (oechem.OECIPBondStereo_E, oechem.OECIPBondStereo_Z)
        and b_stereo == oechem.OECIPBondStereo_NotStereo
    ):
        return Comparison.AMoreThanB
    if (
        b_stereo in (oechem.OECIPBondStereo_E, oechem.OECIPBondStereo_Z)
        and a_stereo == oechem.OECIPBondStereo_NotStereo
    ):
        return Comparison.ALessThanB

    return Comparison.AEqualToB


def compare_recursive(
    a: DigraphNode,
    b: DigraphNode,
    /,
    *,
    current_sphere: Collection[Collection[DigraphNode]] | None = None,
    comparison_function: NodeComparisonFunction = default_comparison,
) -> Comparison:
    """
    Compare two nodes, given a current set of descendant nodes in a given sphere.

    This function recursively compares two nodes by examining their children, and their children, in a
    breadth-first manner as described in the IUPAC Blue Book.

    The initial call (with `current_sphere = None`) compares the two nodes, first directly and then by
    comparing their children pair by pair. If this fails to differentiate the two nodes, the children
    are passed into a second call to `compare_recursive`. The algorithm procedes in this manner, comparing
    equivalent groups of descendents to find a difference.

    Args:
        a: The first node to compare.
        b: The second node to compare.
        current_sphere: Descendants of `a` and `b` in current sphere under consideration, partitioned by
            their ancestors.
        comparison_function: Comparison function used to compare nodes.

    Returns:
        -1 if `a` is less than `b`, 1 if `a` is more than `b` and 0 otherwise.
    """
    if current_sphere is None:
        current_sphere = [[a, b]]

    next_sphere = []

    # `current_sphere` contains lists of nodes. Each list represents all nodes which have equal ancestor
    # trees under the given rule.
    # For each of these groups, the comparison rule is used to split these nodes up and find the first one
    # where a group is either all descendants of a, or all descendants of b
    for nodes in current_sphere:
        # Group nodes by comparison function, ordered from high priority to low priority
        grouped_nodes = sorted_partition(nodes, key=cmp_to_key(comparison_function), reverse=True)
        # Break ties by comparing the nodes' children
        children_key = cmp_to_key(
            partial(compare_children, comparison_function=comparison_function)
        )
        grouped_nodes = refine_partitions(
            grouped_nodes,
            partition_function=partial(sorted_partition, key=children_key, reverse=True),
        )

        # Go through the groups of equivalent nodes from high priority to low priority
        # If any of the groups has either only descendants of a or descendants of b (not both), then
        # the comparison is solved. If not, add the children of the nodes to the next sphere to be examined.
        for group in grouped_nodes:
            descendant_a = a in group or any(a in g.ancestors for g in group)
            descendant_b = b in group or any(b in g.ancestors for g in group)
            if descendant_a and not descendant_b:
                return Comparison.AMoreThanB
            if descendant_b and not descendant_a:
                return Comparison.ALessThanB
            child_nodes = [
                child
                for node in group
                for child in node.children
                if not isinstance(child, PhantomNode)
            ]
            next_sphere.append(child_nodes)

    # Go through each group from high to low and recursively compare
    if (
        len(next_sphere) > 0
        and (
            compare := compare_recursive(
                a, b, current_sphere=next_sphere, comparison_function=comparison_function
            )
        )
        != Comparison.AEqualToB
    ):
        return compare

    return Comparison.AEqualToB


def compare_children(
    a: DigraphNode,
    b: DigraphNode,
    /,
    *,
    comparison_function: NodeComparisonFunction = default_comparison,
) -> Comparison:
    """
    Compare two objects by ordering their children based on a given comparison function and comparing them high to low.

    Args:
        a: The first node to compare.
        b: The second node to compare.
        comparison_function: Comparison function used to compare nodes.

    Returns:
        -1 if `a` is less than `b`, 1 if `a` is more than `b` and 0 otherwise.
    """
    a_children_high_to_low = sorted(a.children, key=cmp_to_key(comparison_function), reverse=True)
    b_children_high_to_low = sorted(b.children, key=cmp_to_key(comparison_function), reverse=True)
    for a_child, b_child in zip(a_children_high_to_low, b_children_high_to_low, strict=True):
        compare = comparison_function(a_child, b_child)
        if compare in [Comparison.ALessThanB, Comparison.AMoreThanB]:
            return compare
    return Comparison.AEqualToB
