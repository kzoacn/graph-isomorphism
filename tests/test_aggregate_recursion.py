from itertools import combinations
from math import factorial

import pytest

from babai.action import Homomorphism, block_action
from babai.aggregation import certificate_table, fullness_structure
from babai.group import PermutationGroup, symmetric_group
from babai.local_certificates import LocalCertificates
from babai.permutation import compose, from_cycles, inverse, pullback
from babai.solver import BabaiSolver
from babai.standard_blocks import subset_action
from babai.tower import BlockTower


def pair_model(m):
    vertices = symmetric_group(m)
    pairs = tuple(combinations(range(m), 2))
    index = {pair: i for i, pair in enumerate(pairs)}
    permutations = [tuple(index[tuple(sorted((g[a], g[b])))] for a, b in pairs) for g in vertices.generators]
    group = PermutationGroup(len(pairs), permutations)
    action = Homomorphism(group, m, vertices.generators)
    return pairs, subset_action(action, tuple(map(frozenset, pairs)))


@pytest.mark.parametrize("edges,expected_kind,expected_order", [
    ([(0, 1), (1, 2), (2, 3), (3, 4), (4, 5)], "small_support", 2),
    ([(0, 1), (0, 2), (1, 2), (3, 4), (3, 5), (4, 5)], "orbit_partition", 72),
])
def test_certificate_tables_feed_the_main_recursion(monkeypatch, edges, expected_kind, expected_order):
    # Exercise aggregation on small actual computations. Only the theorem
    # threshold is lowered for this test. The returned full certificates
    # still undergo their constructive global-automorphism and image
    # checks; no asymptotic conclusion is inferred from this setting.
    monkeypatch.setattr("babai.local_certificates.required_test_size", lambda _: 3)
    pairs, model = pair_model(6)
    source = tuple(int(pair in edges) for pair in pairs)
    vertex_permutation = (3, 1, 5, 0, 4, 2)
    actual = model.action.lift(vertex_permutation)
    target = pullback(source, inverse(actual))
    solver = BabaiSolver()
    engine = LocalCertificates(model.action, solver)
    left = certificate_table(engine, source, 3)
    right = certificate_table(engine, target, 3)
    assert left.structure.kind == right.structure.kind == expected_kind
    solver._towers.append(BlockTower(len(source)))
    try:
        result = solver._aggregate_tables(model, left, right)
    finally:
        solver._towers.pop()
    assert result is not None and result.order == expected_order
    assert result.contains(actual)
    assert all(source[v] == target[result.representative[v]] for v in range(len(source)))


def test_large_nongiant_certificate_orbit_on_the_petersen_graph():
    pairs, model = pair_model(10)
    vertices = tuple(map(frozenset, combinations(range(5), 2)))
    source = tuple(int(not vertices[a].intersection(vertices[b])) for a, b in pairs)
    indices = {subset: i for i, subset in enumerate(vertices)}
    full_group = PermutationGroup(10, [tuple(indices[frozenset(g[v] for v in subset)] for subset in vertices) for g in symmetric_group(5).generators])
    # This is the entire, canonical vertex automorphism group of Petersen.
    # We test the aggregate outcome with actual group generators rather
    # than constructing all high-locality certificates for this fixture.
    left = fullness_structure(full_group)
    assert left.kind == "large_nongiant_orbit"
    p = (7, 2, 8, 0, 9, 5, 1, 6, 4, 3)
    actual = model.action.lift(p)
    target = pullback(source, inverse(actual))
    conjugate = PermutationGroup(10, (compose(compose(inverse(p), g), p) for g in full_group.generators))
    right = fullness_structure(conjugate)
    solver = BabaiSolver()
    solver._towers.append(BlockTower(len(source)))
    try:
        result = solver._nongiant_certificate_orbit(model, source, target, left, right)
    finally:
        solver._towers.pop()
    assert result is not None and result.order == factorial(5)
    assert result.contains(actual)
    assert any(event.kind == "johnson_descent" for event in solver.events)


def test_large_symmetry_rejects_equal_edge_counts_with_different_internal_colors():
    pairs, model = pair_model(10)
    clique = tuple(int(a < 7 and b < 7) for a, b in pairs)
    bipartite = tuple(int((a < 7) != (b < 7)) for a, b in pairs)
    assert sum(clique) == sum(bipartite) == 21
    solver = BabaiSolver()
    assert solver.solve_model(model, clique, bipartite) is None


def test_nonfull_local_relations_reject_a_nonisomorphic_tree(monkeypatch):
    monkeypatch.setattr("babai.local_certificates.required_test_size", lambda _: 3)
    pairs, model = pair_model(6)
    path = {(0, 1), (1, 2), (2, 3), (3, 4), (4, 5)}
    fork = {(0, 1), (1, 2), (2, 3), (2, 4), (4, 5)}
    source = tuple(int(pair in path) for pair in pairs)
    target = tuple(int(pair in fork) for pair in pairs)
    solver = BabaiSolver()
    engine = LocalCertificates(model.action, solver)
    left = certificate_table(engine, source, 3)
    right = certificate_table(engine, target, 3)
    assert left.structure.kind == right.structure.kind == "small_support"
    solver._towers.append(BlockTower(len(source)))
    try:
        assert solver._aggregate_tables(model, left, right) is None
    finally:
        solver._towers.pop()


def test_large_symmetry_pullback_preserves_nontrivial_external_kernel():
    group = PermutationGroup(20, [from_cycles(20, (0, 1)), from_cycles(20, range(0, 20, 2), range(1, 20, 2)), from_cycles(20, (0, 2), (1, 3))])
    action = block_action(group, tuple((2 * i, 2 * i + 1) for i in range(10)))
    model = subset_action(action, tuple(frozenset((v // 2,)) for v in range(20)))
    source = tuple("01" * 7 + "00" * 3)
    p = from_cycles(20, (0, 14), (1, 15), (4, 5))
    target = pullback(source, inverse(p))
    left_domain = tuple(range(7))
    right_domain = tuple(action(p)[v] for v in left_domain)
    solver = BabaiSolver()
    solver._towers.append(BlockTower(20))
    try:
        result = solver._large_symmetry(model, source, target, left_domain, right_domain)
    finally:
        solver._towers.pop()
    assert result is not None and result.order == factorial(7) * factorial(3) * 2**3
    assert result.contains(p)
