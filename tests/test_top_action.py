from itertools import permutations

from babai.action import block_action
from babai.group import PermutationGroup
from babai.permutation import from_cycles
from babai.string_iso import LuksSolver
from babai.top_action import giant_top_action, top_action


def action_on_three_blocks():
    group = PermutationGroup(6, [from_cycles(6, (0, 1)), from_cycles(6, (0, 2, 4), (1, 3, 5)), from_cycles(6, (0, 2), (1, 3))])
    return block_action(group, [(0, 1), (2, 3), (4, 5)])


def test_top_action_generators_lift_to_the_entire_isomorphism_coset():
    action = action_on_three_blocks()
    result = top_action(action, tuple("010101"), tuple("101001"), LuksSolver(), max_orbit_size=2)
    assert result.surjective
    expected = {p for p in permutations(range(6)) if action.domain.contains(p) and all("010101"[v] == "101001"[p[v]] for v in range(6))}
    assert result.isomorphisms.order == len(expected) == 6
    assert all(result.isomorphisms.contains(p) for p in expected)


def test_a_failed_surjectivity_test_does_not_claim_nonisomorphism():
    action = action_on_three_blocks()
    result = top_action(action, tuple("001122"), tuple("001122"), LuksSolver(), max_orbit_size=2)
    assert not result.surjective
    result = giant_top_action(action, tuple("001122"), tuple("221100"), LuksSolver(), max_orbit_size=2)
    assert not result.has_alternating_automorphisms
    result = giant_top_action(action, tuple("010101"), tuple("001011"), LuksSolver(), max_orbit_size=2)
    assert result.has_alternating_automorphisms
    assert result.isomorphisms is None
