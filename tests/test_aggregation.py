from math import factorial

from babai.action import Homomorphism
from babai.aggregation import certificate_table, fullness_structure, local_guide_relations
from babai.configuration import refine_joint
from babai.group import PermutationGroup, symmetric_group
from babai.local_certificates import LocalCertificates
from babai.permutation import from_cycles
from babai.string_iso import LuksSolver


def test_fullness_cases():
    assert fullness_structure(PermutationGroup(10)).kind == "small_support"
    group = PermutationGroup(10, [from_cycles(10, (0, 1), (2, 3), (4, 5))])
    structure = fullness_structure(group)
    assert structure.kind == "orbit_partition"
    assert structure.partition.stabilizer().order == 2**3 * factorial(3) * factorial(4)
    group = PermutationGroup(10, [from_cycles(10, range(7))])
    structure = fullness_structure(group)
    assert structure.kind == "large_nongiant_orbit"
    assert structure.large_orbit == tuple(range(7))


def test_full_certificate_table_generates_global_automorphisms():
    group = symmetric_group(10)
    action = Homomorphism(group, 10, group.generators)
    engine = LocalCertificates(action, LuksSolver())
    table = certificate_table(engine, "a" * 10, 9)
    assert len(table.certificates) == 10
    assert table.full_group.same_group(group)
    assert table.structure.kind == "large_symmetry"
    table = certificate_table(engine, "abcdefghij", 9)
    assert all(not c.full for c in table.certificates.values())
    assert table.full_group.order == 1
    assert table.structure.kind == "small_support"


def test_local_guides_share_equivalence_classes_across_both_inputs():
    labels = (("a", "b", "c", "d"), ("c", "a", "d", "b"))

    def equivalent(i, a, j, b):
        return tuple(labels[i][v] for v in a) == tuple(labels[j][v] for v in b)

    a, b = local_guide_relations((tuple(range(4)), tuple(range(4))), 2, equivalent)
    p = (1, 3, 0, 2)
    assert a.relabel(p) == b
    a, b = refine_joint(a, b).configurations
    assert a.relabel(p) == b
    assert all(len(twins) == 1 for twins in a.twin_classes())
