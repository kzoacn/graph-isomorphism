from collections import Counter
from itertools import combinations, permutations
from math import factorial
from random import Random

import pytest

from babai.action import Homomorphism, block_action
from babai.group import PermutationGroup, symmetric_group
from babai.permutation import from_cycles, inverse, pullback
from babai.solver import BabaiSolver
from babai.standard_blocks import subset_action


def pair_model(m, vertex_group=None):
    vertex_group = vertex_group or symmetric_group(m)
    pairs = tuple(combinations(range(m), 2))
    index = {pair: i for i, pair in enumerate(pairs)}
    images = {tuple(index[tuple(sorted((g[a], g[b])))] for a, b in pairs): g for g in vertex_group.generators}
    group = PermutationGroup(len(pairs), images)
    action = Homomorphism(group, m, (images[g] for g in group.generators))
    return pairs, subset_action(action, tuple(map(frozenset, pairs)))


def verify_result(model, source, target, result, expected_order):
    if not expected_order:
        assert result is None
        return
    assert result is not None and result.order == expected_order
    assert model.action.domain.contains(result.representative)
    assert all(source[v] == target[result.representative[v]] for v in range(len(source)))
    assert all(model.action.domain.contains(g) and all(source[v] == source[g[v]] for v in range(len(source))) for g in result.subgroup.generators)


@pytest.mark.parametrize("m", [3, 4, 5, 6, 7])
def test_modelled_graph_solver_against_all_vertex_permutations(m):
    rng = Random(733 + m)
    pairs, model = pair_model(m)
    pair_index = {pair: i for i, pair in enumerate(pairs)}
    for _ in range(4):
        source = tuple(rng.randrange(2) for _ in pairs)
        p = list(range(m))
        rng.shuffle(p)
        actual = tuple(pair_index[tuple(sorted((p[a], p[b])))] for a, b in pairs)
        target = pullback(source, inverse(actual))
        expected = sum(all(source[i] == target[pair_index[tuple(sorted((p[a], p[b])))]] for i, (a, b) in enumerate(pairs)) for p in permutations(range(m)))
        solver = BabaiSolver()
        result = solver.solve_model(model, source, target)
        verify_result(model, source, target, result, expected)
        assert solver.current_tower is None


def test_regular_nonisomorphism_through_the_giant_recursion():
    pairs, model = pair_model(6)
    cycle_edges = {tuple(sorted((i, (i + 1) % 6))) for i in range(6)}
    triangles = {(0, 1), (0, 2), (1, 2), (3, 4), (3, 5), (4, 5)}
    source = tuple(int(p in cycle_edges) for p in pairs)
    target = tuple(int(p in triangles) for p in pairs)
    assert BabaiSolver().solve_model(model, source, target) is None


@pytest.mark.parametrize("m", [3, 4, 5])
def test_nonconstant_fibers_top_action_and_small_ideal(m):
    n = 2 * m
    group = PermutationGroup(n, [from_cycles(n, (0, 1)), from_cycles(n, range(0, n, 2), range(1, n, 2)), from_cycles(n, (0, 2), (1, 3))])
    action = block_action(group, tuple((2 * i, 2 * i + 1) for i in range(m)))
    model = subset_action(action, tuple(frozenset((i // 2,)) for i in range(n)))
    for source in [tuple("01" * m), tuple("00" + "01" * (m - 1))]:
        p = from_cycles(n, range(0, n, 2), range(1, n, 2))
        target = pullback(source, inverse(p))
        types = Counter(tuple(sorted(source[2 * i:2 * i + 2])) for i in range(m))
        expected = 1
        for block_type, count in types.items():
            expected *= factorial(count) * ((2 if block_type[0] == block_type[1] else 1) ** count)
        solver = BabaiSolver()
        result = solver.solve_model(model, source, target)
        verify_result(model, source, target, result, expected)


@pytest.mark.parametrize("kind", ["cyclic", "affine", "imprimitive", "johnson"])
def test_non_giant_group_structure_reductions(kind):
    if kind == "cyclic":
        group = PermutationGroup(7, [from_cycles(7, range(7))])
    elif kind == "affine":
        group = PermutationGroup(7, [tuple((i + 1) % 7 for i in range(7)), tuple(3 * i % 7 for i in range(7))])
    elif kind == "imprimitive":
        group = PermutationGroup(6, [from_cycles(6, (0, 1)), from_cycles(6, (0, 2, 4), (1, 3, 5)), from_cycles(6, (0, 2), (1, 3))])
    else:
        _, pair_action = pair_model(5)
        group = pair_action.action.domain
    action = Homomorphism(group, group.degree, group.generators)
    model = subset_action(action, tuple(frozenset((v,)) for v in range(group.degree)))
    source = tuple(v % 3 for v in range(group.degree))
    target = pullback(source, inverse(group.generators[0]))
    expected = sum(all(source[v] == target[p[v]] for v in range(group.degree)) for p in group.elements(max_order=1000))
    solver = BabaiSolver(prefer_structural=True)
    result = solver.solve_model(model, source, target)
    verify_result(model, source, target, result, expected)
    expected_event = {"cyclic": "group_individualization", "affine": "ideal_orbit", "imprimitive": "block_descent", "johnson": "johnson_descent"}[kind]
    assert any(event.kind == expected_event for event in solver.events)


def test_the_large_quotient_continuation_is_wired_into_luks(monkeypatch):
    monkeypatch.setattr("babai.string_iso.small_quotient_bound", lambda _: 1)
    cyclic = PermutationGroup(7, [from_cycles(7, range(7))])
    result = BabaiSolver(prefer_structural=True).solve(cyclic, "abcdefg", "gabcdef")
    assert result is not None and result.order == 1


@pytest.mark.parametrize("kind", ["blocks", "johnson"])
def test_union_models_are_used_by_the_recursive_solver(kind):
    if kind == "blocks":
        vertices = PermutationGroup(8, [from_cycles(8, (0, 1)), from_cycles(8, (0, 2, 4, 6), (1, 3, 5, 7)), from_cycles(8, (0, 2), (1, 3))])
        _, model = pair_model(8, vertices)
        expected_ideal, child_degrees = 4, {4, 24}
    else:
        _, first = pair_model(5)
        _, model = pair_model(10, first.action.domain)
        expected_ideal, child_degrees = 5, {15, 30}
    rng = Random(10073)
    source = tuple(rng.randrange(2) for _ in model.supports)
    p = model.action.domain.generators[0]
    target = pullback(source, inverse(p))
    expected = sum(all(source[v] == target[g[v]] for v in range(len(source))) for g in model.action.domain.elements(max_order=1000))
    solver = BabaiSolver(prefer_structural=True)
    result = solver.solve_model(model, source, target)
    verify_result(model, source, target, result, expected)
    assert any(e.ideal_degree == expected_ideal and e.actual_degree in child_degrees for e in solver.events)
