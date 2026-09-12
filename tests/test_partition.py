from itertools import permutations

from babai.partition import ColoredPartition, young_subgroup
from babai.permutation import sign


def test_colored_wreath_partition_and_even_alignment():
    a = ColoredPartition(6, ((0, ((0, 1), (2, 3))), (1, ((4, 5),))))
    b = ColoredPartition(6, ((0, ((1, 4), (0, 5))), (1, ((2, 3),))))
    assert a.stabilizer().order == 16
    result = a.align(b, even=True)
    assert result.order == 8
    assert sign(result.representative) == 1
    for p in permutations(range(6)):
        expected = sign(p) == 1 and all({frozenset(p[v] for v in block) for block in blocks} == {frozenset(block) for block in b.classes[i][1]} for i, (_, blocks) in enumerate(a.classes))
        assert result.contains(p) == expected


def test_young_subgroup_parity_and_discrete_partition():
    assert young_subgroup(6, [[0, 1, 2], [3, 4, 5]], even=True).order == 18
    a = ColoredPartition(2, ((0, ((0,),)), (1, ((1,),))))
    b = ColoredPartition(2, ((0, ((1,),)), (1, ((0,),))))
    assert a.align(b, even=True) is None
