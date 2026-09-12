from itertools import permutations
from random import Random

import pytest

from babai.group import Coset, PermutationGroup, alternating_group, symmetric_group
from babai.permutation import from_cycles, identity, inverse, pullback
from babai.string_iso import GiantActionRequired, LuksSolver


def assert_coset_is_exact(group, source, target, result, window=None, representative=None):
    # Exhaust every ambient permutation, independent of the group enumerator.
    window = range(group.degree) if window is None else window
    representative = identity(group.degree) if representative is None else representative
    ambient = Coset(group, representative)
    expected = {p for p in permutations(range(group.degree)) if ambient.contains(p) and all(source[i] == target[p[i]] for i in window)}
    actual = set() if result is None else {p for p in permutations(range(group.degree)) if result.contains(p)}
    assert actual == expected


def test_random_string_isomorphisms_and_full_cosets():
    rng = Random(47116)
    solver = LuksSolver()
    for n in range(1, 7):
        for _ in range(10):
            generators = []
            for _ in range(2):
                p = list(range(n))
                rng.shuffle(p)
                generators.append(tuple(p))
            group = PermutationGroup(n, generators)
            source = tuple(rng.randrange(3) for _ in range(n))
            target = list(source)
            rng.shuffle(target)
            result = solver.solve(group, source, target)
            assert_coset_is_exact(group, source, target, result)


def test_alternating_parity_obstruction_and_repeated_colors():
    group = alternating_group(6)
    solver = LuksSolver()
    assert solver.solve(group, "abcdef", "bacdef") is None
    result = solver.solve(group, "aabcde", "aacbde")
    assert_coset_is_exact(group, "aabcde", "aacbde", result)
    result = solver.solve(group, "aabbcc", "ccbbaa")
    assert result.order == 4
    assert_coset_is_exact(group, "aabbcc", "ccbbaa", result)


def test_imprimitive_weak_reduction():
    group = PermutationGroup(6, [from_cycles(6, (0, 1)), from_cycles(6, (0, 2, 4), (1, 3, 5)), from_cycles(6, (0, 2), (1, 3))])
    solver = LuksSolver()
    for source, target in [("001122", "221100"), ("001122", "012012"), ("011011", "110110")]:
        result = solver.solve(group, source, target)
        assert_coset_is_exact(group, source, target, result)


def test_window_outside_coset_and_nonfaithful_lifting():
    group = PermutationGroup(6, [from_cycles(6, (0, 1, 2)), from_cycles(6, (3, 4)), from_cycles(6, (3, 4, 5))])
    r = from_cycles(6, (0, 3), (1, 4), (2, 5))
    source = "abcxyz"
    target = "uvwabc"
    result = LuksSolver().solve_window(Coset(group, r), source, target, [0, 1, 2])
    assert result.order == 6
    assert_coset_is_exact(group, source, target, result, [0, 1, 2], r)


def test_missing_giant_continuation_is_not_a_negative_answer(monkeypatch):
    # Force the large-quotient branch on a small instance so its contract
    # can be tested without constructing an enormous permutation action.
    monkeypatch.setattr("babai.string_iso.small_quotient_bound", lambda _: 1)
    cyclic = PermutationGroup(5, [from_cycles(5, range(5))])
    with pytest.raises(GiantActionRequired) as exception:
        LuksSolver().solve(cyclic, "ababa", "aabab")
    assert exception.value.context.quotient.image.order == 5

