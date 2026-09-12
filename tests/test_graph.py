from itertools import combinations, permutations
from random import Random

import pytest

from babai.graph import Graph, graph_isomorphisms, graph_string_problem


def explicit_isomorphisms(a, b):
    # Check the adjacency predicate, independently of the implemented reduction.
    return {p for p in permutations(range(a.n)) if all(((i, j) in a.edges) == (tuple(sorted((p[i], p[j]))) in b.edges) for i, j in combinations(range(a.n), 2))}


def test_graphs_including_tiny_nonfaithful_pair_actions():
    rng = Random(9986)
    for n in range(6):
        for _ in range(5):
            a = Graph(n, [pair for pair in combinations(range(n), 2) if rng.randrange(2)])
            p = list(range(n))
            rng.shuffle(p)
            b = a.relabel(tuple(p))
            result = graph_isomorphisms(a, b)
            expected = explicit_isomorphisms(a, b)
            assert result is not None
            assert result.order == len(expected)
            assert {p for p in permutations(range(n)) if result.contains(p)} == expected
    assert graph_string_problem(Graph(2), Graph(2)).action.kernel.order == 2


def test_regular_nonisomorphic_graphs():
    cycle = Graph(6, [(i, (i + 1) % 6) for i in range(6)])
    triangles = Graph(6, [(0, 1), (1, 2), (0, 2), (3, 4), (4, 5), (3, 5)])
    assert graph_isomorphisms(cycle, triangles) is None


def test_input_validation():
    with pytest.raises(ValueError):
        Graph(3, [(0, 3)])
    with pytest.raises(ValueError):
        Graph(3, [(1, 1)])
    assert graph_isomorphisms(Graph(1), Graph(2)) is None
    assert graph_isomorphisms(Graph(2), Graph(2, [(0, 1)])) is None
