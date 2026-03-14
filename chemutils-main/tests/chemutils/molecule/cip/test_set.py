from functools import partial

import pytest

from chemutils.molecule.cip.set import refine_partitions, sorted_partition


@pytest.mark.parametrize(
    ["original_data", "key", "partitioned"],
    [
        (
            ["ben", "fred", "david", "betty", "anna", "carrie", "david", "alex", "betty"],
            lambda s: s[0],
            [
                ["anna", "alex"],
                ["ben", "betty", "betty"],
                ["carrie"],
                ["david", "david"],
                ["fred"],
            ],
        ),
        (
            [2, 0, 3, 4, 2, 1, 1, 3],
            None,
            [
                [0],
                [1, 1],
                [2, 2],
                [3, 3],
                [4],
            ],
        ),
    ],
)
@pytest.mark.parametrize("reverse", [True, False])
def test_sorted_partition(original_data, key, partitioned, reverse):
    if reverse:
        assert sorted_partition(original_data, key=key, reverse=reverse) == list(
            reversed(partitioned)
        )
    else:
        assert sorted_partition(original_data, key=key, reverse=reverse) == partitioned


@pytest.mark.parametrize(
    ["original_groups", "key", "partitioned"],
    [
        (
            [["Ben", "fred", "glen"], ["james"], ["Anna", "emma"]],
            lambda s: s[-1],
            [
                ["fred"],
                ["Ben", "glen"],
                ["james"],
                ["Anna", "emma"],
            ],
        )
    ],
)
def test_refine_partitions(original_groups, key, partitioned):
    partition_function = partial(sorted_partition, key=key)
    assert refine_partitions(original_groups, partition_function=partition_function) == partitioned
