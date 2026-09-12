"""Babai's local certificates (Helfgott section 6.1).

The recursive work is explicitly restricted to kernel orbits of size at
most n/k. A certificate is not a global string-isomorphism answer. Its
aggregation into the main giant-action recursion is a separate stage.
"""

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from math import factorial

from .action import Homomorphism, restriction_action
from .group import Coset, PermutationGroup, _points, alternating_group, merge_isomorphism_cosets
from .partition import young_subgroup
from .permutation import compose, identity
from .string_iso import LuksSolver


def required_test_size(degree: int) -> int:
    """The strict Unaffected Stabilizer Theorem threshold, in integers."""
    return max(9, degree.bit_length() + 2)


def contains_alternating(group: PermutationGroup) -> bool:
    return alternating_group(group.degree).is_subgroup_of(group)


def affected_points(action: Homomorphism) -> tuple[int, ...]:
    """Points whose stabilizer does not map onto a group containing A_k."""
    if not contains_alternating(action.image):
        raise ValueError("affected-point test requires a giant image")
    affected = []
    for orbit in action.domain.orbits():
        stabilizer = action.domain.pointwise_stabilizer((orbit[0],))
        image = action.restrict(stabilizer).image
        if not contains_alternating(image):
            affected.extend(orbit)
    return tuple(sorted(affected))


def _kernel_window_isomorphisms(action: Homomorphism, coset: Coset, source: Sequence, target: Sequence, window: tuple[int, ...], solver: LuksSolver, *, max_orbit_size: int, max_image_order: int) -> Coset | None:
    """Weak reduction under an explicitly bounded image, then small orbits."""
    if not action.domain.same_group(coset.subgroup):
        raise ValueError("action must be defined on the candidate subgroup")
    kernel = action.kernel
    orbits = kernel.orbits(window)
    if any(len(orbit) > max_orbit_size for orbit in orbits):
        raise AssertionError("affected-kernel orbit bound violated")

    def branches():
        for image in action.image.elements(max_order=max_image_order):
            lift = action.lift(image)
            if lift is None:
                raise AssertionError("image element has no lift")
            result = Coset(kernel, compose(lift, coset.representative))
            for orbit in orbits:
                result = solver.solve_window(result, source, target, orbit)
                if result is None:
                    break
            yield result

    return merge_isomorphism_cosets(branches())


def _update_window(action: Homomorphism, coset: Coset, source: tuple, target: tuple, window: tuple[int, ...], solver: LuksSolver, k: int) -> Coset | None:
    group = coset.subgroup
    # Polynomial shortcuts. These do not invoke an unrestricted SI oracle.
    if all(source[p] == target[coset.representative[p]] for p in window) and all(source[g[p]] == source[p] for g in group.generators for p in window):
        return coset
    restriction = restriction_action(group, window)
    image_order = restriction.image.order
    if image_order == factorial(len(window)) or contains_alternating(restriction.image):
        # solve_window dispatches directly to the natural S_w/A_w case.
        return solver.solve_window(coset, source, target, window)
    return _kernel_window_isomorphisms(action, coset, source, target, window, solver, max_orbit_size=group.degree // k, max_image_order=factorial(k))


@dataclass(frozen=True)
class LocalCertificate:
    test_set: tuple[int, ...]
    full: bool
    window: tuple[int, ...]
    # Exact Aut_(G_T)^window(source); image is on the ordered test_set.
    automorphisms: PermutationGroup
    image: PermutationGroup
    # Present precisely when full. It fixes the complement of the window
    # pointwise and consists of global source-string automorphisms.
    certificate_group: PermutationGroup | None
    windows: tuple[tuple[int, ...], ...]


class LocalCertificates:
    """Local certificate construction and comparison for a given giant map."""

    def __init__(self, action: Homomorphism, solver: LuksSolver):
        if not contains_alternating(action.image):
            raise ValueError("local certificates require a giant representation")
        self.action = action
        self.solver = solver

    def test_action(self, test_set: Iterable[int]) -> Homomorphism:
        test_set = tuple(sorted(_points(test_set, self.action.image_degree)))
        complement = tuple(p for p in range(self.action.image_degree) if p not in test_set)
        even = self.action.image.order != factorial(self.action.image_degree)
        stabilizer = young_subgroup(self.action.image_degree, (test_set, complement), even=even)
        restricted = self.action.restrict(self.action.preimage(stabilizer))
        return restricted.then(restriction_action(restricted.image, test_set))

    def build(self, source: Sequence, test_set: Iterable[int]) -> LocalCertificate:
        source = tuple(source)
        n = self.action.domain.degree
        if len(source) != n:
            raise ValueError("string length differs from the action degree")
        test_set = tuple(sorted(_points(test_set, self.action.image_degree)))
        k = len(test_set)
        # Strict inequality k > max(8, 2+log2(n)), using integer arithmetic.
        if k < required_test_size(n):
            raise ValueError("the Unaffected Stabilizer Theorem requires k > max(8, 2+log2(n))")
        test_action = self.test_action(test_set)
        current = test_action.domain
        window: tuple[int, ...] = ()
        windows = []
        while True:
            restricted = test_action.restrict(current)
            if not contains_alternating(restricted.image):
                return LocalCertificate(test_set, False, window, current, restricted.image, None, tuple(windows))
            affected = affected_points(restricted)
            new_window = tuple(sorted(set(window).union(affected)))
            if new_window == window:
                certificate = current.pointwise_stabilizer(p for p in range(n) if p not in window)
                if not contains_alternating(test_action.restrict(certificate).image):
                    raise AssertionError("Unaffected Stabilizer Theorem conclusion violated")
                if not all(all(source[p] == source[g[p]] for p in range(n)) for g in certificate.generators):
                    raise AssertionError("full certificate is not a global automorphism group")
                return LocalCertificate(test_set, True, window, current, restricted.image, certificate, tuple(windows))
            result = _update_window(restricted, Coset(current, identity(n)), source, source, new_window, self.solver, k)
            if result is None:
                raise AssertionError("an automorphism problem must contain the identity")
            current = result.subgroup
            window = new_window
            windows.append(window)

    def compare(self, source: Sequence, target: Sequence, left: LocalCertificate, right: LocalCertificate, *, ordered: tuple[tuple[int, ...], tuple[int, ...]] | None = None) -> Coset | None:
        """Isomorphisms of the canonically windowed strings mapping T to T'.

        Certificates must have been built by this LocalCertificates instance
        from the corresponding strings. The windows are replayed together;
        all recursive calls obey the same affected-kernel size bound.
        """
        source, target = tuple(source), tuple(target)
        if len(source) != self.action.domain.degree or len(target) != len(source):
            raise ValueError("invalid string lengths")
        if len(left.test_set) != len(right.test_set):
            raise ValueError("test sets must have the same size")
        if len(left.windows) != len(right.windows) or left.full != right.full:
            return None
        image_coset = self.action.image.transporter(left.test_set, right.test_set)
        if image_coset is None:
            return None
        representative = self.action.lift(image_coset.representative)
        if representative is None:
            raise AssertionError("test-set alignment has no lift")
        test_action = self.test_action(left.test_set)
        result = Coset(test_action.domain, representative)
        for window, target_window in zip(left.windows, right.windows):
            if {result.representative[p] for p in window} != set(target_window):
                return None
            action = test_action.restrict(result.subgroup)
            result = _update_window(action, result, source, target, window, self.solver, len(left.test_set))
            if result is None:
                return None
        if ordered is not None:
            source_tuple, target_tuple = ordered
            if len(source_tuple) != len(left.test_set) or set(source_tuple) != set(left.test_set) or len(target_tuple) != len(right.test_set) or set(target_tuple) != set(right.test_set):
                raise ValueError("ordered tuples must order the certificate test sets")
            action = self.action.restrict(result.subgroup)
            representative_image = self.action(result.representative)
            inverse_image = [0] * len(representative_image)
            for a, b in enumerate(representative_image):
                inverse_image[b] = a
            matching = action.image.transporter(source_tuple, tuple(inverse_image[b] for b in target_tuple))
            matching = action.preimage_coset(matching)
            return None if matching is None else matching.multiply_right(result.representative)
        return result
