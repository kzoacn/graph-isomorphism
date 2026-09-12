from itertools import combinations, permutations
from math import comb

import pytest

from babai.action import Homomorphism, block_action
from babai.group import PermutationGroup, symmetric_group
from babai.permutation import from_cycles, identity
from babai.standard_blocks import descend_action, standard_blocks


def tuple_representation(m, ordered=False):
    points = tuple(permutations(range(m), 2) if ordered else combinations(range(m), 2))
    index = {p: i for i, p in enumerate(points)}
    group = symmetric_group(m)
    generators = [tuple(index[tuple(g[v] for v in p) if ordered else tuple(sorted(g[v] for v in p))] for p in points) for g in group.generators]
    actual = PermutationGroup(len(points), generators)
    return points, Homomorphism(actual, m, group.generators)


@pytest.mark.parametrize("ordered", [False, True])
def test_standard_supports_and_ordered_pair_fibers(ordered):
    points, action = tuple_representation(7, ordered)
    model = standard_blocks(action)
    assert model.subset_size == 2
    assert model.supports == tuple(map(frozenset, points))
    assert len(model.blocks) == comb(7, 2)
    assert all(len(block) == (2 if ordered else 1) for block in model.blocks)
    colors = (0, 0, 0, 1, 1, 1, 1)
    windows = model.color_windows(colors)
    assert sum(len(w) for _, w in windows) == len(points)
    assert all(3 * len(w) <= 2 * len(points) for _, w in windows)
    assert model.fiber_colors(tuple(0 for _ in points)) == (0,) * comb(7, 2)
    if ordered:
        assert model.fiber_colors(tuple(range(len(points)))) is None


def test_descending_an_action_requires_its_kernel_to_factor():
    group = PermutationGroup(6, [from_cycles(6, (0, 1, 2), (3, 4, 5)), from_cycles(6, (0, 1), (3, 4))])
    action = Homomorphism(group, 3, [(1, 2, 0), (1, 0, 2)])
    projection, descended = descend_action(action, (3, 4, 5))
    assert descended.image.same_group(symmetric_group(3))
    assert all(descended(projection(g)) == action(g) for g in group.generators)
    group = PermutationGroup(6, [from_cycles(6, (0, 1, 2)), from_cycles(6, (3, 4, 5))])
    action = Homomorphism(group, 3, [(1, 2, 0), identity(3)])
    with pytest.raises(ValueError, match="does not factor"):
        descend_action(action, (3, 4, 5))

