from math import factorial

import pytest

from babai.action import Homomorphism, block_action
from babai.group import Coset, PermutationGroup, symmetric_group
from babai.local_certificates import LocalCertificates, _kernel_window_isomorphisms, affected_points, contains_alternating
from babai.permutation import from_cycles, identity, inverse, pullback
from babai.string_iso import LuksSolver


def identity_action(n):
    group = symmetric_group(n)
    return Homomorphism(group, n, group.generators)


def test_full_and_nonfull_certificates_with_theorem_threshold():
    action = identity_action(10)
    certificates = LocalCertificates(action, LuksSolver())
    full = certificates.build("aaaaaaaaaz", range(9))
    assert full.full
    assert full.window == tuple(range(9))
    assert full.certificate_group.order == factorial(9)
    assert all(g[9] == 9 for g in full.certificate_group.generators)
    nonfull = certificates.build("aaaaabbbbz", range(9))
    assert not nonfull.full
    assert nonfull.window == tuple(range(9))
    assert nonfull.automorphisms.order == factorial(5) * factorial(4)
    assert not contains_alternating(nonfull.image)
    with pytest.raises(ValueError, match="Theorem requires"):
        certificates.build("aaaaaaaaaz", range(8))


def test_affected_and_unaffected_points():
    certificates = LocalCertificates(identity_action(11), LuksSolver())
    test_action = certificates.test_action(range(9))
    assert affected_points(test_action) == tuple(range(9))
    assert test_action.kernel.order == 2
    full = certificates.build("aaaaaaaaaxy", range(9))
    # The certificate need not preserve x,y through its window, but it
    # must fix these unaffected points before claiming global symmetry.
    assert full.full and full.certificate_group.order == factorial(9)


def test_certificate_comparison_and_ordered_tuple_constraint():
    action = identity_action(10)
    certificates = LocalCertificates(action, LuksSolver())
    source = tuple("aaaaabbbbz")
    permutation = from_cycles(10, (0, 9, 6, 2), (3, 8))
    target = pullback(source, inverse(permutation))
    left = certificates.build(source, range(9))
    right = certificates.build(target, [permutation[p] for p in range(9)])
    result = certificates.compare(source, target, left, right)
    assert result is not None
    assert result.contains(permutation)
    assert result.order == factorial(5) * factorial(4)
    assert all(source[p] == target[result.representative[p]] for p in left.window)
    ordered = certificates.compare(source, target, left, right, ordered=(tuple(range(9)), tuple(permutation[p] for p in range(9))))
    assert ordered is not None and ordered.order == 1 and ordered.representative == permutation
    different = tuple("aaaaaabbbz")
    right = certificates.build(different, range(9))
    assert certificates.compare(source, different, left, right) is None


def test_kernel_reduction_only_sends_small_windows_to_recursion():
    group = PermutationGroup(6, [from_cycles(6, (0, 1)), from_cycles(6, (0, 2, 4), (1, 3, 5)), from_cycles(6, (0, 2), (1, 3))])
    action = block_action(group, [(0, 1), (2, 3), (4, 5)])

    class Recorder(LuksSolver):
        windows = []

        def solve_window(self, coset, source, target, window):
            self.windows.append(tuple(window))
            return super().solve_window(coset, source, target, window)

    solver = Recorder()
    result = _kernel_window_isomorphisms(action, Coset(group, identity(6)), "001122", "221100", tuple(range(6)), solver, max_orbit_size=2, max_image_order=6)
    assert result is not None and result.order == 8
    assert all(len(window) <= 2 for window in solver.windows)
    assert all("001122"[i] == "221100"[result.representative[i]] for i in range(6))

