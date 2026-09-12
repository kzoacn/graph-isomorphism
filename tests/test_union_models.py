from itertools import combinations
from math import comb

from babai.action import Homomorphism, block_action
from babai.group import PermutationGroup, symmetric_group
from babai.permutation import from_cycles
from babai.standard_blocks import descend_action, subset_action


def test_block_union_descent_retains_complete_uniform_subset_models():
    vertices = PermutationGroup(8, [from_cycles(8, (0, 1)), from_cycles(8, (0, 2, 4, 6), (1, 3, 5, 7)), from_cycles(8, (0, 2), (1, 3))])
    pairs = tuple(combinations(range(8), 2))
    index = {pair: i for i, pair in enumerate(pairs)}
    images = [tuple(index[tuple(sorted((g[a], g[b])))] for a, b in pairs) for g in vertices.generators]
    actual = PermutationGroup(len(pairs), images)
    vertex_action = Homomorphism(actual, 8, vertices.generators)
    block_projection = block_action(vertices, [(0, 1), (2, 3), (4, 5), (6, 7)])
    action = vertex_action.then(block_projection)
    supports = tuple(frozenset((a // 2, b // 2)) for a, b in pairs)
    sizes = {}
    for v, support in enumerate(supports):
        sizes.setdefault(len(support), []).append(v)
    assert sorted(map(len, sizes.values())) == [4, 24]
    for size, window in sizes.items():
        _, descended = descend_action(action, tuple(window))
        model = subset_action(descended, tuple(supports[v] for v in window))
        assert model.subset_size == size
        assert len(model.blocks) == comb(4, size)
        assert len({len(block) for block in model.blocks}) == 1


def test_johnson_union_descent_and_complement_normalization():
    atoms = symmetric_group(5)
    johnson_vertices = tuple(map(frozenset, combinations(range(5), 2)))
    index = {subset: i for i, subset in enumerate(johnson_vertices)}
    pair_vertices = tuple(combinations(range(10), 2))
    pair_index = {pair: i for i, pair in enumerate(pair_vertices)}
    actual_generators = []
    for g in atoms.generators:
        first = tuple(index[frozenset(g[v] for v in subset)] for subset in johnson_vertices)
        actual_generators.append(tuple(pair_index[tuple(sorted((first[a], first[b])))] for a, b in pair_vertices))
    actual = PermutationGroup(45, actual_generators)
    action = Homomorphism(actual, 5, atoms.generators)
    supports = tuple(johnson_vertices[a] | johnson_vertices[b] for a, b in pair_vertices)
    for union_size, expected_size in [(3, 30), (4, 15)]:
        window = tuple(v for v, support in enumerate(supports) if len(support) == union_size)
        assert len(window) == expected_size
        _, descended = descend_action(action, window)
        model = subset_action(descended, tuple(supports[v] for v in window))
        assert model.subset_size == 5 - union_size
        assert all(len(block) == 3 for block in model.blocks)
