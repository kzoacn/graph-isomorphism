from collections import Counter
from itertools import combinations, product

from babai.configuration import Configuration, individualization_compatible, refine_joint
from babai.graph import Graph, graph_isomorphisms


def examples():
    rook = Graph(16, [(a, b) for a, b in combinations(range(16), 2) if a // 4 == b // 4 or a % 4 == b % 4])
    steps = ((1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (-1, -1))
    cayley = Graph(16, [(4 * i + j, 4 * ((i + di) % 4) + (j + dj) % 4) for i in range(4) for j in range(4) for di, dj in steps])
    return rook, cayley


def test_equal_wl_profiles_do_not_produce_a_false_isomorphism():
    rook, shrikhande = examples()
    structures = [Configuration(16, 2, tuple(0 if a == b else 1 if tuple(sorted((a, b))) in graph.edges else 2 for a, b in product(range(16), repeat=2))) for graph in (rook, shrikhande)]
    a, b = refine_joint(*structures).configurations
    assert Counter(a.colors) == Counter(b.colors)
    # An independent obstruction: one graph has 4-cliques, the other does not.
    has_k4 = lambda graph: any(all(tuple(sorted(edge)) in graph.edges for edge in combinations(vertices, 2)) for vertices in combinations(range(16), 4))
    assert has_k4(rook) and not has_k4(shrikhande)
    assert not individualization_compatible(a, b, (0,), (0,))
    assert graph_isomorphisms(rook, shrikhande) is None
