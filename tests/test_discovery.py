from itertools import product
from random import Random

from babai.alignment import align_splits
from babai.configuration import Configuration, is_coherent, refine_joint
from babai.discovery import discover_relation
from babai.group import PermutationGroup
from babai.orbitals import break_double_transitivity, is_doubly_transitive, orbital_configuration
from babai.permutation import from_cycles
from babai.split_johnson import Choice


def test_relation_discovery_replays_all_individualization_choices():
    n = 9
    raw = Configuration(n, 2, tuple(0 if a == b else 1 if a // 3 == b // 3 or a % 3 == b % 3 else 2 for a, b in product(range(n), repeat=2)))
    p = list(range(n))
    Random(885).shuffle(p)
    a, b = refine_joint(raw, raw.relabel(tuple(p))).configurations
    left = discover_relation(a)
    right = discover_relation(b, fixed=tuple(p[v] for v in left.choices[0].vertices), split_plan=tuple(Choice(c.kind, tuple(p[v] for v in c.vertices)) for c in left.choices[1:]))
    alignment = align_splits(left, right)
    assert alignment is not None and alignment.contains(tuple(p))


def test_orbitals_and_transitivity_reduction_on_affine_group():
    # AGL(1,7) is 2-transitive, but its point stabilizer is not 2-transitive.
    group = PermutationGroup(7, [tuple((i + 1) % 7 for i in range(7)), tuple(3 * i % 7 for i in range(7))])
    assert group.order == 42
    assert is_doubly_transitive(group)
    assert break_double_transitivity(group) == (0,)
    assert is_coherent(orbital_configuration(group))
    cyclic = PermutationGroup(7, [from_cycles(7, range(7))])
    assert break_double_transitivity(cyclic) == ()
