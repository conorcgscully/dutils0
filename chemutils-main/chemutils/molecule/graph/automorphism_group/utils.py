from collections.abc import Callable
from typing import TypeVar

_T = TypeVar("_T")


def group_comparator(*, items: list[_T], comparator: Callable[[_T, _T], bool]) -> list[list[_T]]:
    """
    Group a set of items based on a function that determines if two items must be in the same group.

    This is effectively the same as grouping the vertices of a graph into components: the items are the vertices,
    and the comparator function determines if two vertices are connected by an edge.

    The relation does not have to be transitive - if `func(A, B)` and `func(B, C)` are both `True`, then A,
    B and C will be in the same group. This is true even if `func(A, C)` is `False`.

    The comparator function should be symmetric however, namely `func(x, y) = func(y, x)`.

    In group theory terms, the comparator should be a *partial equivalence relation*.:
    * reflexive, i.e. `func(x, x)` should always be `True`.
    * symmetric, i.e. `func(x, y) = func(y, x)`.
    """
    colours = list(range(len(items)))

    for i in range(len(items)):
        similar_colours = set()
        for j in range(i):
            if colours[j] in similar_colours:
                continue
            if comparator(items[i], items[j]):
                similar_colours.add(colours[j])
        if similar_colours:
            new_colour = min(similar_colours)
            colours[i] = new_colour
            for k in range(i):
                if colours[k] in similar_colours:
                    colours[k] = new_colour

    unique_colours = set(colours)
    return [
        [items[i] for i in range(len(items)) if colours[i] == colour] for colour in unique_colours
    ]
