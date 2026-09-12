from collections import Counter
from fractions import Fraction
from itertools import product
from random import Random

import pytest

from babai.configuration import Configuration, is_clique, is_coherent, refine_joint
from babai.design import design_at, design_lemma


def graph_configuration(n, edges):
    edges = {frozenset(e) for e in edges}
    return Configuration(n, 2, tuple(0 if a == b else 1 if frozenset((a, b)) in edges else 2 for a, b in product(range(n), repeat=2)))


def test_joint_refinement_is_equivariant():
    rng = Random(5181)
    for n in range(1, 9):
        x = graph_configuration(n, [(a, b) for a in range(n) for b in range(a + 1, n) if rng.randrange(2)])
        p = list(range(n))
        rng.shuffle(p)
        y = x.relabel(tuple(p))
        result = refine_joint(x, y)
        a, b = result.configurations
        assert a.relabel(tuple(p)) == b
        assert is_coherent(a) and is_coherent(b)
        assert result.rounds <= 2 * n**2
        assert refine_joint(a, b).rounds == 1


def test_distinct_original_colors_are_not_accidentally_identified():
    a, b = refine_joint(Configuration(3, 2, (7,) * 9), Configuration(3, 2, (12,) * 9)).configurations
    assert not set(a.colors).intersection(b.colors)


def test_two_wl_separates_cycle_from_two_triangles():
    cycle = graph_configuration(6, [(i, (i + 1) % 6) for i in range(6)])
    triangles = graph_configuration(6, [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3)])
    a, b = refine_joint(cycle, triangles).configurations
    assert Counter(a.colors) != Counter(b.colors)


def test_coherence_all_coordinate_maps_and_conditioning():
    rng = Random(34)
    x = Configuration(4, 3, tuple(rng.randrange(3) for _ in range(4**3)))
    assert not is_coherent(x)
    x = refine_joint(x).configurations[0]
    assert is_coherent(x)
    assert is_coherent(x.skeleton(2))
    assert is_coherent(x.condition((2,)))
    for _, vertices in x.vertex_classes():
        assert is_coherent(x.induced(vertices))


def test_design_lemma_and_twins():
    cycle = refine_joint(graph_configuration(7, [(i, (i + 1) % 7) for i in range(7)])).configurations[0]
    witness = design_lemma(cycle)
    assert witness.fixed == ()
    assert witness.dominant == tuple(range(7))
    assert is_coherent(witness.configuration)
    assert not is_clique(witness.configuration)
    assert design_at(cycle, (0,)).dominant is None
    complete = refine_joint(graph_configuration(7, [(a, b) for a in range(7) for b in range(a + 1, 7)])).configurations[0]
    assert complete.twin_classes() == (tuple(range(7)),)
    with pytest.raises(ValueError, match="twin class"):
        design_lemma(complete)
    x = Configuration(6, 2, tuple(6 * a + b for a, b in product(range(6), repeat=2)))
    assert design_lemma(x, Fraction(2, 3)).dominant is None


def test_empty_structure_and_validation():
    assert refine_joint(Configuration(0, 2, ())).rounds == 0
    with pytest.raises(ValueError):
        Configuration(2, 2, (0, 0))
    with pytest.raises(ValueError):
        Configuration(2, 2, (0, 1, 1, 0)).color((0, 2))
