from itertools import combinations, product
from math import comb
from random import Random

import pytest

from babai.configuration import Configuration, refine_joint
from babai.johnson import johnson_configuration, recognize_johnson
from babai.permutation import from_cycles


@pytest.mark.parametrize("m,k", [(5, 2), (6, 2), (7, 2), (7, 3), (8, 3), (9, 4)])
def test_recovers_relabelled_johnson_schemes(m, k):
    x = johnson_configuration(m, k)
    rng = Random(100 * m + k)
    p = list(range(x.degree))
    rng.shuffle(p)
    x = x.relabel(tuple(p))
    # Relation labels need not agree with the numerical Johnson distances.
    color_permutation = list(range(k + 1))
    rng.shuffle(color_permutation)
    x = Configuration(x.degree, 2, tuple(color_permutation[c] for c in x.colors))
    model = recognize_johnson(x)
    assert model is not None
    assert (model.atoms, model.subset_size) == (m, k)
    assert len(set(model.subsets)) == comb(m, k)
    atom_permutation = from_cycles(m, (0, 2, 1))
    vertex_permutation = model.induce(atom_permutation)
    assert x.relabel(vertex_permutation) == x
    assert model.recover(vertex_permutation) == atom_permutation
    assert model.recover(from_cycles(x.degree, (0, 1))) is None


def test_wl_of_a_johnson_graph_recovers_the_scheme():
    full = johnson_configuration(7, 3)
    adjacency = Configuration(full.degree, 2, tuple(min(c, 2) for c in full.colors))
    coherent = refine_joint(adjacency).configurations[0]
    assert recognize_johnson(coherent) is not None


def test_rejects_false_parameter_matches_and_corruption():
    x = johnson_configuration(7, 2)
    colors = list(x.colors)
    colors[1] = colors[x.degree] = 2 if colors[1] == 1 else 1
    assert recognize_johnson(Configuration(x.degree, 2, tuple(colors))) is None
    # A 21-cycle has the right vertex count for J(7,2), but is not Johnson.
    cycle = Configuration(21, 2, tuple(0 if a == b else 1 if (a - b) % 21 in (1, 20) else 2 for a, b in product(range(21), repeat=2)))
    assert recognize_johnson(cycle) is None
    assert recognize_johnson(Configuration(5, 2, tuple(int(a != b) for a, b in product(range(5), repeat=2)))) is None
