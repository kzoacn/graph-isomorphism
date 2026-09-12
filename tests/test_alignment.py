from math import factorial
from random import Random

from babai.alignment import align_splits, split_stabilizer
from babai.configuration import Configuration
from babai.johnson import johnson_configuration
from babai.permutation import sign
from babai.split_johnson import Choice, split_or_johnson


def test_aligns_johnson_schemes_and_generates_their_full_group():
    x = johnson_configuration(7, 2)
    p = list(range(x.degree))
    Random(463).shuffle(p)
    a = split_or_johnson(x)
    b = split_or_johnson(x.relabel(tuple(p)))
    result = align_splits(a, b)
    assert result is not None and result.order == factorial(7)
    assert result.contains(tuple(p))
    assert all(x.relabel(g) == x for g in split_stabilizer(a).generators)
    even = align_splits(a, b, even=True)
    assert even is not None and even.order == factorial(7) // 2
    assert sign(even.representative) == 1


def test_alignment_respects_individualized_points():
    x = Configuration(9, 2, tuple(0 if a == b else 1 if a // 3 == b // 3 or a % 3 == b % 3 else 2 for a in range(9) for b in range(9)))
    p = (3, 1, 8, 2, 5, 0, 6, 4, 7)
    a = split_or_johnson(x)
    plan = tuple(Choice(choice.kind, tuple(p[v] for v in choice.vertices)) for choice in a.choices)
    b = split_or_johnson(x.relabel(p), plan=plan)
    result = align_splits(a, b)
    assert result is not None and result.contains(p)
    for g in result.subgroup.generators:
        assert all(g[v] == v for v in a.fixed_points)
    assert all(result.representative[v] == p[v] for v in a.fixed_points)
