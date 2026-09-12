from fractions import Fraction
from itertools import combinations, product
from random import Random

import pytest

from babai.coherent_split import BipartiteGraph
from babai.configuration import Configuration, refine_joint
from babai.johnson import johnson_configuration
from babai.split_johnson import Choice, IncompatibleChoices, bipartite_split, split_or_johnson


def normalized_parts(parts):
    return tuple(frozenset(frozenset(block) for block in blocks) for blocks in parts)


def verify_split(result):
    assert {v for blocks in result.parts for block in blocks for v in block} == set(result.domain)
    assert sum(len(block) for blocks in result.parts for block in blocks) == len(result.domain)
    special = set() if result.johnson is None else set(result.johnson.vertices)
    for blocks in result.parts:
        vertices = {v for block in blocks for v in block}
        if vertices == special:
            continue
        assert all(len(block) <= result.threshold for block in blocks)
        if len(vertices) > result.threshold:
            assert all(len(block) >= 2 for block in blocks)


def rook_configuration(m):
    return Configuration(m * m, 2, tuple(0 if a == b else 1 if a // m == b // m or a % m == b % m else 2 for a, b in product(range(m * m), repeat=2)))


def test_johnson_and_imprimitive_outcomes():
    result = split_or_johnson(johnson_configuration(7, 2))
    verify_split(result)
    assert result.johnson is not None and result.johnson.model.atoms == 7
    assert not result.choices
    # The octahedral scheme has a canonical partition into antipodal pairs.
    subsets = tuple(map(frozenset, combinations(range(4), 2)))
    x = Configuration(6, 2, tuple(2 - len(a.intersection(b)) for a in subsets for b in subsets))
    result = split_or_johnson(x)
    verify_split(result)
    assert result.johnson is None and not result.choices
    assert all(len(block) == 2 for blocks in result.parts for block in blocks)


@pytest.mark.parametrize("m", [3, 7])
def test_complete_split_is_equivariant_when_choices_are_replayed(m):
    x = rook_configuration(m)
    result = split_or_johnson(x)
    verify_split(result)
    assert result.johnson is None
    assert result.choices[0].kind == "point"
    if m == 7:
        assert any(choice.kind == "all" for choice in result.choices)
    p = list(range(x.degree))
    Random(673 + m).shuffle(p)
    mapped_plan = tuple(Choice(choice.kind, tuple(p[v] for v in choice.vertices)) for choice in result.choices)
    other = split_or_johnson(x.relabel(tuple(p)), plan=mapped_plan)
    assert normalized_parts(tuple(tuple(tuple(p[v] for v in block) for block in blocks) for blocks in result.parts)) == normalized_parts(other.parts)
    with pytest.raises(IncompatibleChoices):
        split_or_johnson(x, plan=())


def incidence_graph(m, subsets):
    subsets = tuple(map(frozenset, subsets))
    n = len(subsets)
    return BipartiteGraph(tuple(range(n)), tuple(frozenset((n + i,)) for i in range(m)), subsets)


def test_hypergraph_complete_and_balanced_cases(monkeypatch):
    # Bypass only the early small-right optimization to exercise these
    # exact algebraic branches on tractable examples. No timing claim is
    # inferred from tests using the altered cutoff.
    monkeypatch.setattr("babai.split_johnson._small_right_limit", lambda _: 0)
    result = bipartite_split(incidence_graph(5, combinations(range(5), 2)))
    verify_split(result)
    assert result.johnson is not None
    assert result.johnson.model.subset_size == 2
    result = bipartite_split(incidence_graph(6, combinations(range(6), 3)))
    verify_split(result)
    assert result.johnson is None
    assert all(len(block) == 2 for blocks in result.parts for block in blocks)


@pytest.mark.parametrize("m,missing", [(6, [(0, 1), (2, 3), (4, 5)]), (7, [(0, 1)])])
def test_design_and_twin_reductions(monkeypatch, m, missing):
    monkeypatch.setattr("babai.split_johnson._small_right_limit", lambda _: 0)
    graph = incidence_graph(m, [pair for pair in combinations(range(m), 2) if pair not in missing])
    result = bipartite_split(graph)
    verify_split(result)
    assert result.johnson is None
    if m == 6:
        assert any(choice.kind == "design" for choice in result.choices)
    else:
        assert not result.choices


def test_rejects_a_clique_and_invalid_bipartite_hypotheses():
    clique = Configuration(5, 2, tuple(int(a != b) for a, b in product(range(5), repeat=2)))
    with pytest.raises(ValueError, match="nontrivial"):
        split_or_johnson(clique)
    graph = BipartiteGraph((0, 1, 2), (frozenset((3,)),), (frozenset(),) * 3)
    with pytest.raises(ValueError, match="hypotheses"):
        bipartite_split(graph)
