"""
"""

import pytest

import pandas as pd

from BioDendro.cluster import Tree


# Test Tree methods


@pytest.mark.parametrize("sample,expected", [
    ([3, 3, 3], "3.0000_3.0000_3.0000"),
    ([1.234, 2.345678, 3.45678910], "2.3455_1.2340_3.4568"),
    ])
def test_Tree__bin_name(sample, expected):
    actual = Tree._bin_name(sample)
    assert actual == expected
    return


@pytest.mark.parametrize("sample,expected", [
    ([1, 3, 5, 7], [0, 1, 2, 3]),  # Case where no clusters
    ([1, 3, 4, 7], [0, 1, 3]),  # Single cluster
    ([1, 2, 5, 6], [0, 2]),  # match at start or end
    ([1, 3, 4, 5], [0, 1]),  # Multiple matches
    ([1, 3, 4.5, 5], [0, 1, 2]),  # Check threshold equality case
    ([1, 3, 4.4999, 5], [0, 1]),  # Check threshold near equality case
    ])
def test_Tree__bin_starts(sample, expected):
    actual = Tree._bin_starts(pd.Series(sample), threshold=1.5)

    assert len(actual) == len(expected)
    for act, exp in zip(actual, expected):
        assert act == exp
    return


@pytest.mark.parametrize("sample,starts,expected", [
    (
        [1, 3, 5],
        [0, 1, 2],
        ["1.0000_1.0000_1.0000",
         "3.0000_3.0000_3.0000",
         "5.0000_5.0000_5.0000"]
    ),
    (
        [1, 3, 4, 7],
        [0, 1, 3],
        ["1.0000_1.0000_1.0000",
         "3.5000_3.0000_4.0000",
         "3.5000_3.0000_4.0000",
         "7.0000_7.0000_7.0000"]
    ),
    (
        [1, 2, 5, 6],
        [0, 2],
        ["1.5000_1.0000_2.0000",
         "1.5000_1.0000_2.0000",
         "5.5000_5.0000_6.0000",
         "5.5000_5.0000_6.0000"]
    ),  # match at start or end
    ])
def test_Tree__bin_names(sample, starts, expected):
    actual = Tree._bin_names(sample, starts)

    assert len(actual) == len(expected)

    for act, exp in zip(actual, expected):
        assert act == exp
    return


@pytest.mark.parametrize("sample,bins,expected", [
    (
        {"sample": [1, 2, 3, 4]},
        ['a', 'b', 'c', 'd'],
        [[True, False, False, False],
         [False, True, False, False],
         [False, False, True, False],
         [False, False, False, True]],
    ),
    (
        {"sample": [1, 2, 3, 4]},
        ['a', 'a', 'b', 'c'],
        [[True, False, False],
         [True, False, False],
         [False, True, False],
         [False, False, True]],
    ),
    ])
def test_Tree__pivot(sample, bins, expected):
    actual = Tree._pivot(pd.DataFrame(sample), bins, "sample")

    assert all([e in actual.columns for e in bins])
    assert all([e in actual.index for e in sample["sample"]])

    for i, (_, row) in enumerate(actual.iterrows()):
        print("i", i, "act", row)
        print("i", i, "exp", expected[i])
        for j, act in enumerate(row):
            print(i, j)
            assert act == expected[i][j]
    return


@pytest.mark.parametrize("columns,values,expected", [
    (
        ['a', 'b', 'c', 'd'],
        [[True, False, False, False],
         [False, True, False, False],
         [False, False, True, False],
         [False, False, False, True]],
        ['a', 'b', 'c', 'd'],
    ),
    (
        ['a', 'b', 'c'],
        [[True, False, False],
         [True, False, False],
         [False, False, False],
         [False, False, False]],
        ['a'],
    ),
    (
        ['a', 'a', 'b'],
        [[False, False, False],
         [False, False, False],
         [False, False, False],
         [False, False, False]],
        [],
    ),
    ])
def test_Tree__exclude_false_columns(columns, values, expected):
    df = pd.DataFrame([dict(zip(columns, row)) for row in values])

    actual = Tree._exclude_false_columns(df)
    assert len(actual.columns) == len(expected)

    for col in expected:
        assert col in actual.columns
    return
