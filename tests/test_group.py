"""Small explicit groups are an independent oracle for the implicit algorithms."""

from itertools import permutations
from math import factorial
from random import Random

import pytest

from babai.group import Coset, EnumerationLimitExceeded, PermutationGroup, alternating_group, symmetric_group
from babai.permutation import compose, from_cycles, identity, inverse, pullback, sign


def closure(n, generators):
    found = {tuple(range(n))}
    pending = list(found)
    while pending:
        p = pending.pop()
        for q in generators:
            # Independent definition, deliberately not compose().
            candidate = tuple(q[p[i]] for i in range(n))
            if candidate not in found:
                found.add(candidate)
                pending.append(candidate)
    return found


def random_group(rng, n):
    generators = []
    for _ in range(rng.randrange(4)):
        g = list(range(n))
        rng.shuffle(g)
        generators.append(tuple(g))
    return PermutationGroup(n, generators)


def test_right_action_convention():
    p = from_cycles(4, (0, 1, 2))
    q = from_cycles(4, (1, 3))
    assert compose(p, q) == (3, 2, 0, 1)
    assert compose(p, inverse(p)) == identity(4)
    assert pullback("abcd", p) == ("b", "c", "a", "d")
    assert sign(p) == 1
    assert sign(q) == -1


@pytest.mark.parametrize("n", range(8))
def test_chain_against_explicit_group(n):
    rng = Random(81413 + n)
    for _ in range(12):
        group = random_group(rng, n)
        explicit = closure(n, group.generators)
        assert group.order == len(explicit)
        assert set(group.elements(max_order=factorial(n))) == explicit
        assert group.compact().same_group(group)
        assert len(group.compact().generators) <= n * (n - 1) // 2
        for _ in range(8):
            candidate = list(range(n))
            rng.shuffle(candidate)
            assert group.contains(candidate) == (tuple(candidate) in explicit)
        points = list(range(n))
        rng.shuffle(points)
        points = points[:min(n, 2)]
        stabilizer = group.pointwise_stabilizer(points)
        expected = {g for g in explicit if all(g[p] == p for p in points)}
        assert set(stabilizer.elements(max_order=factorial(n))) == expected
        expected_orbits = {frozenset(g[p] for g in explicit) for p in range(n)}
        assert {frozenset(orbit) for orbit in group.orbits()} == expected_orbits


def test_transporters_and_non_normal_cosets():
    rng = Random(427)
    for _ in range(8):
        group = random_group(rng, 5)
        explicit = closure(5, group.generators)
        for source in [(0,), (2, 0)]:
            for target in permutations(range(5), len(source)):
                coset = group.transporter(source, target)
                expected = {g for g in explicit if tuple(g[p] for p in source) == target}
                assert (coset is None) == (not expected)
                if coset:
                    actual = {compose(g, coset.representative) for g in coset.subgroup.elements(max_order=120)}
                    assert actual == expected
        subgroup = group.pointwise_stabilizer([2])
        representatives = group.coset_representatives(subgroup, max_index=5)
        cosets = [{g for g in explicit if Coset(subgroup, r).contains(g)} for r in representatives]
        assert set().union(*cosets) == explicit
        assert sum(map(len, cosets)) == len(explicit)


def test_large_group_is_implicit():
    group = symmetric_group(16)
    assert group.order == factorial(16)
    assert group.pointwise_stabilizer([7, 2, 14]).order == factorial(13)
    assert alternating_group(12).order == factorial(12) // 2
    with pytest.raises(EnumerationLimitExceeded):
        list(group.elements(max_order=1000))
    with pytest.raises(EnumerationLimitExceeded):
        group.coset_representatives(group.pointwise_stabilizer([0, 1]), max_index=16)


def test_validation():
    with pytest.raises(ValueError):
        PermutationGroup(3, [(1, 1, 0)])
    with pytest.raises(ValueError):
        symmetric_group(3).orbits([0, 1])
    with pytest.raises(ValueError):
        symmetric_group(3).pointwise_stabilizer([0, 0])
    with pytest.raises(ValueError):
        symmetric_group(3).transporter([0, 1], [0])
