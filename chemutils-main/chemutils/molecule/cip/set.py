from collections.abc import Callable, Collection
from itertools import groupby
from typing import TYPE_CHECKING, TypeVar, overload

if TYPE_CHECKING:
    from _typeshed import SupportsRichComparison

    _TComparableItem = TypeVar("_TComparableItem", bound=SupportsRichComparison)
else:
    _TComparableItem = TypeVar("_TComparableItem")

_TItem = TypeVar("_TItem")

PartitionFunction = Callable[[Collection[_TItem]], Collection[list[_TItem]]]
"""
Callable that partitions a collection into disjoint collections, maintaining input order.
"""


# Without a key function, items must be comparable
@overload
def sorted_partition(
    items: Collection[_TComparableItem], /, *, key: None = ..., reverse: bool = ...
) -> list[list[_TComparableItem]]: ...


# With a key function, items can be any type
@overload
def sorted_partition(
    items: Collection[_TItem],
    /,
    *,
    key: Callable[[_TItem], _TComparableItem] = ...,
    reverse: bool = ...,
) -> list[list[_TItem]]: ...


def sorted_partition(
    items: Collection[_TItem],
    /,
    *,
    key: Callable[[_TItem], _TComparableItem] | None = None,
    reverse: bool = False,
) -> list[list[_TItem]]:
    """
    Sort a set of items and group them together by the same key.

    This is the concept of a 'weak ordering' in set theory - the key allows each item to be ranked
    as equal than, less than or more than any other item, allowing for multiple items to have the
    same key.

    This is best paired with the `functools.cmp_to_key` function, which converts a comparison function
    into a key function.

    For example:
    ```
    grouped_sort([4, 1, 2, 4, 3, 1])
    # [[1, 1], [2], [3], [4, 4]]
    ```

    The order of values in each partition is the same order as the input values.

    Args:
        items: Collection of items to sort and group.
        key: Custom key function for sorting and grouping items. If not provided, the items themselves
            are used as a key.
        reverse: Should the key be sorted in reverse order.

    Returns:
        List of lists, where each list contains items which are equal under the
        given key and the lists are arranged in increasing value of the key.
    """
    items = sorted(items, key=key, reverse=reverse)  # type: ignore
    return [list(group) for _, group in groupby(items, key=key)]


def refine_partitions(
    groups: Collection[Collection[_TItem]], /, *, partition_function: PartitionFunction[_TItem]
) -> list[list[_TItem]]:
    """
    Refine an (ordered) grouping (or partitioning) of a collection of items by passing each group to a partitioning function.

    This takes each partition and breaks it into smaller parititions in place - the original order of the partitions is maintained.

    For example we could first partition a list of characters by letter using `sorted_partition`:

    ```
    [B, c, b, A, b, A, a, c] -> [[A, a, A], [B, b, b], [c, c]]
    ```

    To further partition these partitions to separate upper from lower case, we would use `refine_partition`:

    ```
    [[A, A, a], [B, b, b], [c, c]] -> [[A, A], [a], [B], [b, b], [c, c]]
    ```

    Args:
        groups: Collection of collections of items.
        partition_function: Function to partition each of the collections in `groups`.

    Returns:
        List of lists, where each group in `groups` may have been split in-place into multiple lists using `partition_function`.
    """
    new_groups = []
    for group in groups:
        if len(group) == 1:
            new_groups.append(list(group))
        else:
            for new_group in partition_function(group):
                new_groups.append(new_group)
    return new_groups
