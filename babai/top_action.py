"""Lift image generators to string isomorphisms (Babai section 9)."""

from dataclasses import dataclass
from math import factorial

from .action import Homomorphism
from .group import Coset, PermutationGroup, alternating_group, coset_hull
from .local_certificates import contains_alternating
from .permutation import identity, pullback
from .string_iso import LuksSolver


@dataclass(frozen=True)
class TopActionResult:
    surjective: bool
    isomorphisms: Coset | None


def lift_image_isomorphisms(action: Homomorphism, source: tuple, target: tuple, image: tuple[int, ...], solver: LuksSolver, *, max_orbit_size: int) -> Coset | None:
    if len(source) != action.domain.degree or len(target) != len(source):
        raise ValueError("invalid string lengths")
    lift = action.lift(image)
    if lift is None:
        return None
    orbits = action.kernel.orbits()
    if any(len(orbit) > max_orbit_size for orbit in orbits):
        raise ValueError("kernel orbit exceeds the permitted recursive window size")
    result = Coset(action.kernel, lift)
    for orbit in orbits:
        result = solver.solve_window(result, source, target, orbit)
        if result is None:
            return None
    return result


def top_action(action: Homomorphism, source: tuple, target: tuple, solver: LuksSolver, *, max_orbit_size: int) -> TopActionResult:
    """Decide whether Iso_G(x,y) projects onto the *entire* image group.

    Testing the identity and image generators suffices. A failed test is
    not a nonisomorphism certificate unless surjectivity of Aut_G(x) is
    already established by the caller.
    """
    fibers = []
    for image in (identity(action.image_degree),) + action.image.generators:
        result = lift_image_isomorphisms(action, source, target, image, solver, max_orbit_size=max_orbit_size)
        if result is None:
            return TopActionResult(False, None)
        fibers.append(result)
    # These fibers generate the full isomorphism coset. Their union need
    # not itself be a coset, so use the hull rather than a union routine.
    return TopActionResult(True, coset_hull(fibers))


@dataclass(frozen=True)
class GiantTopActionResult:
    has_alternating_automorphisms: bool
    isomorphisms: Coset | None


def giant_top_action(action: Homomorphism, source: tuple, target: tuple, solver: LuksSolver, *, max_orbit_size: int) -> GiantTopActionResult:
    """Verify alternating top symmetry, then solve SI if it is present."""
    if not contains_alternating(action.image):
        raise ValueError("a giant image is required")
    even_group = action.preimage(alternating_group(action.image_degree))
    even_action = action.restrict(even_group)
    automorphisms = top_action(even_action, source, source, solver, max_orbit_size=max_orbit_size)
    if not automorphisms.surjective:
        return GiantTopActionResult(False, None)
    solutions = []
    for representative in action.domain.coset_representatives(even_group, max_index=2):
        result = top_action(even_action, source, pullback(target, representative), solver, max_orbit_size=max_orbit_size)
        if result.surjective:
            solutions.append(result.isomorphisms.multiply_right(representative))
    return GiantTopActionResult(True, coset_hull(solutions))

