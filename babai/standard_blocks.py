"""Construct the standard Johnson blocks of a giant representation.

Babai Theorem 8.5.1 and Proposition 8.5.4. The support of a point is the
complement of the unique large orbit of the image of its stabilizer.
The Jordan--Liebeck subgroup inclusions are checked constructively, so
small instances satisfying the same conclusion can also be used safely.
"""

from collections import Counter
from dataclasses import dataclass
from math import comb

from .action import Homomorphism, restriction_action
from .group import PermutationGroup, alternating_group
from .local_certificates import contains_alternating
from .permutation import identity


class StandardBlocksUnavailable(ValueError):
    """A small representation need not admit the large-degree conclusion."""


@dataclass(frozen=True)
class StandardBlocks:
    action: Homomorphism
    supports: tuple[frozenset[int], ...]
    subsets: tuple[frozenset[int], ...]
    blocks: tuple[tuple[int, ...], ...]
    subset_size: int

    def color_windows(self, ideal_colors: tuple[int, ...]) -> tuple[tuple[tuple[tuple[int, int], ...], tuple[int, ...]], ...]:
        if len(ideal_colors) != self.action.image_degree:
            raise ValueError("ideal coloring has the wrong degree")
        windows = {}
        for v, support in enumerate(self.supports):
            signature = tuple(sorted(Counter(ideal_colors[p] for p in support).items()))
            windows.setdefault(signature, []).append(v)
        return tuple((signature, tuple(vertices)) for signature, vertices in sorted(windows.items()))

    def fiber_colors(self, string: tuple) -> tuple | None:
        """Return colors on standard blocks if the string is constant on each."""
        if len(string) != self.action.domain.degree:
            raise ValueError("string has the wrong degree")
        if any(any(string[v] != string[block[0]] for v in block) for block in self.blocks):
            return None
        return tuple(string[block[0]] for block in self.blocks)

    def restrict_group(self, subgroup: PermutationGroup) -> "StandardBlocks":
        return StandardBlocks(self.action.restrict(subgroup), self.supports, self.subsets, self.blocks, self.subset_size)


def subset_action(action: Homomorphism, supports: tuple[frozenset[int], ...]) -> StandardBlocks:
    """Verify an equivariant complete-subset model, for an arbitrary image group.

    The actual group need not act transitively. The *domain* must consist
    of equally many points over every t-subset of the ideal domain. This
    combinatorial invariant survives the structural reductions even when
    the image ceases to be giant. Complements normalize t to at most m/2.
    """
    n, m = action.domain.degree, action.image_degree
    supports = tuple(frozenset(support) for support in supports)
    if len(supports) != n or not supports or len({len(s) for s in supports}) != 1:
        raise ValueError("provide one equal-size subset per actual point")
    t = len(supports[0])
    if not 0 < t < m or any(any(type(p) is not int or not 0 <= p < m for p in support) for support in supports):
        raise ValueError("supports must be nonempty proper ideal subsets")
    if 2 * t > m:
        universe = frozenset(range(m))
        supports = tuple(universe - support for support in supports)
        t = m - t
    classes = {}
    for v, support in enumerate(supports):
        classes.setdefault(support, []).append(v)
    if len(classes) != comb(m, t) or len({len(block) for block in classes.values()}) != 1:
        raise ValueError("the subset model must be complete with uniform fibers")
    for g in action.domain.generators:
        image = action(g)
        if any(supports[g[v]] != frozenset(image[p] for p in supports[v]) for v in range(n)):
            raise ValueError("the proposed subset model is not equivariant")
    subsets = tuple(sorted(classes, key=lambda support: tuple(sorted(support))))
    return StandardBlocks(action, supports, subsets, tuple(tuple(classes[s]) for s in subsets), t)


def standard_blocks(action: Homomorphism) -> StandardBlocks:
    group = action.domain
    m = action.image_degree
    if group.degree == 0 or len(group.orbits()) != 1 or not contains_alternating(action.image):
        raise ValueError("standard blocks require a transitive giant representation")
    stabilizer_image = action.restrict(group.pointwise_stabilizer((0,))).image
    large = next((orbit for orbit in stabilizer_image.orbits() if 2 * len(orbit) > m), None)
    if large is None:
        raise StandardBlocksUnavailable("the point-stabilizer image has no dominant orbit")
    support = frozenset(range(m)) - frozenset(large)
    if not support or 2 * len(support) >= m:
        raise StandardBlocksUnavailable("the standard support is empty or too large")
    generators = []
    for g in alternating_group(len(large)).generators:
        p = list(range(m))
        for i, v in enumerate(large):
            p[v] = large[g[i]]
        generators.append(tuple(p))
    if not PermutationGroup(m, generators).is_subgroup_of(stabilizer_image):
        raise StandardBlocksUnavailable("Jordan--Liebeck stabilizer inclusion does not hold")
    transversal = group.orbit_transversal(0)
    supports = tuple(frozenset(action(transversal[v])[p] for p in support) for v in range(group.degree))
    for g in group.generators:
        image = action(g)
        if any(supports[g[v]] != frozenset(image[p] for p in supports[v]) for v in range(group.degree)):
            raise AssertionError("standard supports are not equivariant")
    classes = {}
    for v, s in enumerate(supports):
        classes.setdefault(s, []).append(v)
    subsets = tuple(sorted(classes, key=lambda s: tuple(sorted(s))))
    blocks = tuple(tuple(classes[s]) for s in subsets)
    if len(subsets) != comb(m, len(support)) or len({len(block) for block in blocks}) != 1:
        raise AssertionError("standard blocks do not form a full uniform Johnson action")
    return StandardBlocks(action, supports, subsets, blocks, len(support))


def descend_action(action: Homomorphism, window: tuple[int, ...]) -> tuple[Homomorphism, Homomorphism]:
    """Factor a giant/auxiliary action through an invariant-window restriction.

    Returns the restriction pi:G->G|window and psi:G|window->Sym(m).
    Rejects an action that does not factor; dropping such a kernel would
    otherwise silently give an incorrect homomorphism.
    """
    restriction = restriction_action(action.domain, window)
    e = identity(action.image_degree)
    if any(action(g) != e for g in restriction.kernel.generators):
        raise ValueError("ideal action does not factor through this window")
    images = []
    for g in restriction.image.generators:
        lift = restriction.lift(g)
        if lift is None:
            raise AssertionError("restriction generator has no lift")
        images.append(action(lift))
    descended = Homomorphism(restriction.image, action.image_degree, images)
    return restriction, descended
