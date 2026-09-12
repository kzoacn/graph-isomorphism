from itertools import combinations
from math import factorial

import pytest

from babai.action import Homomorphism, block_action, restriction_action
from babai.blocks import maximal_block_system, nontrivial_block_system
from babai.group import PermutationGroup, alternating_group, symmetric_group
from babai.permutation import compose, from_cycles, identity, sign


def test_sign_homomorphism_kernel_and_lifts():
    group = symmetric_group(7)
    phi = Homomorphism(group, 2, (identity(2) if sign(g) == 1 else (1, 0) for g in group.generators))
    assert phi.image.order == 2
    assert phi.kernel.order == factorial(7) // 2
    assert phi.kernel.same_group(alternating_group(7))
    for p in [(0, 1), (1, 0)]:
        lift = phi.lift(p)
        assert lift is not None and group.contains(lift)
        assert phi(lift) == p
    assert phi.preimage(PermutationGroup(2)).same_group(phi.kernel)
    assert phi.preimage(phi.image).same_group(group)


def test_pair_action_and_constructive_inverse():
    group = symmetric_group(6)
    pairs = list(combinations(range(6), 2))
    index = {pair: i for i, pair in enumerate(pairs)}
    images = [tuple(index[tuple(sorted((g[a], g[b])))] for a, b in pairs) for g in group.generators]
    phi = Homomorphism(group, len(pairs), images)
    assert phi.kernel.order == 1
    for g in list(group.elements(max_order=720))[::29]:
        p = phi(g)
        assert phi.lift(p) == g
        assert p == tuple(index[tuple(sorted((g[a], g[b])))] for a, b in pairs)
    assert phi.lift(from_cycles(len(pairs), (0, 1))) is None
    restricted = phi.restrict(group.pointwise_stabilizer([2]))
    assert restricted.image.order == 120


def test_relation_consistency_is_checked():
    cyclic3 = PermutationGroup(3, [from_cycles(3, (0, 1, 2))])
    with pytest.raises(ValueError, match="do not define a homomorphism"):
        Homomorphism(cyclic3, 2, [(1, 0)])
    with pytest.raises(ValueError):
        Homomorphism(cyclic3, 2, [])


def test_nonfaithful_restriction_lifts_whole_coset():
    group = PermutationGroup(6, [from_cycles(6, (0, 1, 2)), from_cycles(6, (3, 4)), from_cycles(6, (3, 4, 5))])
    phi = restriction_action(group, [2, 0, 1])
    assert phi.image.order == 3
    assert phi.kernel.order == 6
    target = phi.image.transporter([0], [2])
    result = phi.preimage_coset(target)
    explicit = list(group.elements(max_order=18))
    assert {g for g in explicit if result.contains(g)} == {g for g in explicit if phi(g)[0] == 2}
    assert phi.then(Homomorphism(phi.image, 1, [identity(1) for _ in phi.image.generators])).kernel.same_group(group)


def test_blocks():
    for n in [2, 5, 7]:
        assert maximal_block_system(symmetric_group(n)) == tuple((i,) for i in range(n))
    # The cyclic degree-12 action has nontrivial blocks, then a primitive quotient.
    group = PermutationGroup(12, [from_cycles(12, range(12))])
    blocks = maximal_block_system(group)
    assert 1 < len(blocks) < 12
    phi = block_action(group, blocks)
    assert nontrivial_block_system(phi.image) is None
    assert phi.kernel.order * phi.image.order == group.order
    for g in group.elements(max_order=12):
        assert {frozenset(g[p] for p in b) for b in blocks} == {frozenset(b) for b in blocks}


def test_wreath_product_blocks():
    group = PermutationGroup(8, [from_cycles(8, (0, 1)), from_cycles(8, (0, 2, 4, 6), (1, 3, 5, 7)), from_cycles(8, (0, 2), (1, 3))])
    assert group.order == 2**4 * factorial(4)
    assert maximal_block_system(group) == ((0, 1), (2, 3), (4, 5), (6, 7))
    phi = block_action(group, maximal_block_system(group))
    assert phi.image.same_group(symmetric_group(4))
    assert phi.kernel.order == 16
    with pytest.raises(ValueError):
        block_action(group, [(0, 2), (1, 3), (4, 5), (6, 7)])
    with pytest.raises(ValueError):
        restriction_action(group, [0, 1])
    with pytest.raises(ValueError):
        maximal_block_system(PermutationGroup(4))
