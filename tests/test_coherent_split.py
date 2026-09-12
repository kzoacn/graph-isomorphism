from itertools import combinations
from random import Random

from babai.coherent_split import coherent_split
from babai.configuration import Configuration


def subset_fibers(m, k, l):
    left = tuple(map(frozenset, combinations(range(m), k)))
    right = tuple(map(frozenset, combinations(range(m), l)))
    domain = left + right
    signatures = [(int(a >= len(left)), int(b >= len(left)), len(x.intersection(y))) for a, x in enumerate(domain) for b, y in enumerate(domain)]
    palette = {s: i for i, s in enumerate(sorted(set(signatures)))}
    return Configuration(len(domain), 2, tuple(palette[s] for s in signatures)), tuple(range(len(left))), tuple(range(len(left), len(domain)))


def test_corrected_primitive_case_and_choice_equivariance():
    x, left, right = subset_fibers(7, 3, 2)
    result = coherent_split(x, left, right)
    assert result.case == "primitive_right_fix"
    assert len(result.fixed) == 1
    graph = result.bipartite
    assert 2 * len(graph.left) > len(left)
    assert 2 * len(graph.right) < len(right)
    assert graph.is_nontrivial_semiregular()
    assert 2 * max(map(len, graph.twin_classes())) <= len(graph.left)
    # Fix the corresponding point on an arbitrarily relabeled structure.
    p = list(range(x.degree))
    Random(75).shuffle(p)
    other = coherent_split(x.relabel(tuple(p)), [p[v] for v in left], [p[v] for v in right], choose=lambda _: p[result.fixed[0]])
    assert {p[v] for v in graph.left} == set(other.bipartite.left)
    edges = {(p[v], frozenset(p[w] for w in graph.right[j])) for v, neighbors in zip(graph.left, graph.neighborhoods) for j in neighbors}
    other_edges = {(v, other.bipartite.right[j]) for v, neighbors in zip(other.bipartite.left, other.bipartite.neighborhoods) for j in neighbors}
    assert edges == other_edges


def test_imprimitive_left_partition():
    x, left, right = subset_fibers(6, 3, 2)
    result = coherent_split(x, left, right)
    assert result.case == "left_imprimitive"
    assert result.bipartite is None and not result.fixed
    assert {v for blocks in result.parts for block in blocks for v in block} == set(left)
    assert all(len(block) == 2 for blocks in result.parts for block in blocks)


def test_right_block_contraction():
    m = 6
    left_subsets = tuple(map(frozenset, combinations(range(m), 2)))
    right_pairs = tuple((i, bit) for i in range(m) for bit in range(2))
    nleft = len(left_subsets)
    n = nleft + len(right_pairs)
    signatures = []
    for a in range(n):
        for b in range(n):
            if a < nleft and b < nleft:
                s = (0, 0, len(left_subsets[a].intersection(left_subsets[b])), 0)
            elif a < nleft:
                s = (0, 1, int(right_pairs[b - nleft][0] in left_subsets[a]), 0)
            elif b < nleft:
                s = (1, 0, int(right_pairs[a - nleft][0] in left_subsets[b]), 0)
            else:
                i, bit = right_pairs[a - nleft]
                j, other_bit = right_pairs[b - nleft]
                s = (1, 1, int(i == j), int(bit == other_bit))
            signatures.append(s)
    # Put the relation joining copies of one point first, so it is the
    # canonical disconnected relation selected by this fixture's palette.
    palette_values = sorted(set(signatures), key=lambda s: (s != (1, 1, 1, 0), s))
    palette = {s: i for i, s in enumerate(palette_values)}
    x = Configuration(n, 2, tuple(palette[s] for s in signatures))
    result = coherent_split(x, range(nleft), range(nleft, n))
    assert result.case == "contracted_right_blocks"
    assert not result.fixed
    assert len(result.bipartite.right) == m
    assert all(len(fiber) == 2 for fiber in result.bipartite.right)
    assert result.bipartite.is_nontrivial_semiregular()


def test_balanced_cross_coloring():
    x, left, right = subset_fibers(9, 3, 2)
    result = coherent_split(x, left, right)
    assert result.case == "balanced_cross_colors"
    assert result.bipartite is None
    assert len(result.fixed) == 1
    assert all(2 * len(block) <= len(left) for blocks in result.parts for block in blocks)
