"""Luks reductions used by the Babai string-isomorphism algorithm.

This module implements exact orbit-by-orbit and small-quotient reductions.
The large primitive quotient is a separate, explicit continuation. Until
that continuation is integrated, such inputs raise GiantActionRequired;
there is no exponential backtracking fallback and no guessed answer.
"""

from collections import Counter
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass
from math import factorial

from .action import Homomorphism, block_action, restriction_action
from .blocks import maximal_block_system
from .group import Coset, PermutationGroup, merge_isomorphism_cosets
from .permutation import compose, from_cycles, identity, pullback, sign
from .tower import BlockTower


@dataclass(frozen=True)
class GiantActionContext:
    group: PermutationGroup
    source: tuple
    target: tuple
    blocks: tuple[tuple[int, ...], ...]
    quotient: Homomorphism
    tower: BlockTower | None = None


class GiantActionRequired(NotImplementedError):
    """The remaining Babai stages are needed for a large primitive quotient."""

    def __init__(self, context: GiantActionContext):
        self.context = context
        super().__init__(f"Babai giant-action continuation is not integrated: degree={context.group.degree}, blocks={len(context.blocks)}, quotient_order={context.quotient.image.order}")


def small_quotient_bound(degree: int) -> int:
    """n^(1+ceil(log2(n))); an integer upper bound for the small-group case."""
    if degree < 1:
        raise ValueError("quotient degree must be positive")
    return degree ** (1 + (degree - 1).bit_length())


def _natural_isomorphisms(group: PermutationGroup, source: tuple, target: tuple, even: bool) -> Coset | None:
    classes: dict[object, list[int]] = {}
    target_classes: dict[object, list[int]] = {}
    for i, c in enumerate(source):
        classes.setdefault(c, []).append(i)
    for i, c in enumerate(target):
        target_classes.setdefault(c, []).append(i)
    if {c: len(v) for c, v in classes.items()} != {c: len(v) for c, v in target_classes.items()}:
        return None
    representative = list(range(group.degree))
    generators = []
    for color, vertices in classes.items():
        for i, j in zip(vertices, target_classes[color]):
            representative[i] = j
        for v in vertices[1:]:
            generators.append(from_cycles(group.degree, (vertices[0], v)))
    representative = tuple(representative)
    automorphisms = PermutationGroup(group.degree, generators)
    if even:
        if sign(representative) == -1:
            if not generators:
                return None
            representative = compose(generators[0], representative)
        parity = Homomorphism(automorphisms, 2, ((1, 0) for _ in automorphisms.generators))
        automorphisms = parity.kernel
    return Coset(automorphisms, representative)


class LuksSolver:
    """A recursive SI engine; the giant continuation must satisfy the same API.

    It returns all isomorphisms as Aut_G(source) * representative. An
    installed continuation is responsible for both correctness and its
    recurrence bound. None means proven nonisomorphism, never unsupported.
    """

    def __init__(self, giant_handler: Callable[[GiantActionContext], Coset | None] | None = None):
        self.giant_handler = giant_handler
        self.calls = 0
        self._towers: list[BlockTower] = []

    @property
    def current_tower(self) -> BlockTower | None:
        return self._towers[-1] if self._towers else None

    def solve(self, group: PermutationGroup, source: Sequence, target: Sequence, *, tower: BlockTower | None = None) -> Coset | None:
        source, target = tuple(source), tuple(target)
        if len(source) != group.degree or len(target) != group.degree:
            raise ValueError("strings must have the permutation group's degree")
        if tower is None:
            tower = self.current_tower
            if tower is None or tower.degree != group.degree:
                tower = BlockTower(group.degree)
        tower.validate_for(group)
        self._towers.append(tower)
        try:
            return self._solve_instance(group, source, target)
        finally:
            self._towers.pop()

    def _solve_instance(self, group: PermutationGroup, source: tuple, target: tuple) -> Coset | None:
        self.calls += 1
        if Counter(source) != Counter(target):
            return None
        if source == target and all(all(source[g[p]] == source[p] for p in range(group.degree)) for g in group.generators):
            return Coset(group, identity(group.degree))
        if not group.generators:
            return None

        # Natural symmetric and alternating actions have a direct solution.
        order = group.order
        if order == factorial(group.degree):
            return _natural_isomorphisms(group, source, target, even=False)
        if order * 2 == factorial(group.degree) and all(sign(g) == 1 for g in group.generators):
            return _natural_isomorphisms(group, source, target, even=True)

        orbits = group.orbits()
        if len(orbits) > 1:
            result = Coset(group, identity(group.degree))
            for orbit in orbits:
                result = self.solve_window(result, source, target, orbit)
                if result is None:
                    return None
            return result

        blocks = maximal_block_system(group, self.current_tower.top)
        self._towers[-1] = self.current_tower.extend(blocks)
        quotient = block_action(group, blocks)
        bound = small_quotient_bound(len(blocks))
        if quotient.image.order > bound:
            context = GiantActionContext(group, source, target, blocks, quotient, self.current_tower)
            if self.giant_handler is None:
                raise GiantActionRequired(context)
            return self.giant_handler(context)

        # Weak Luks reduction: enumerate the *bounded quotient*, not G.
        kernel = quotient.kernel

        def branches():
            for image in quotient.image.elements(max_order=bound):
                representative = quotient.lift(image)
                if representative is None:
                    raise AssertionError("quotient element has no lift")
                result = self.solve(kernel, source, pullback(target, representative))
                yield None if result is None else result.multiply_right(representative)

        return merge_isomorphism_cosets(branches())

    def solve_window(self, coset: Coset, source: Sequence, target: Sequence, window: Iterable[int]) -> Coset | None:
        """Impose a new invariant window on candidates H*r (Luks's chain rule)."""
        group = coset.subgroup
        if len(source) != group.degree or len(target) != group.degree:
            raise ValueError("strings must have the permutation group's degree")
        window = tuple(window)
        action = restriction_action(group, window)  # Also checks invariance.
        source_window = tuple(source[p] for p in window)
        target_window = tuple(target[coset.representative[p]] for p in window)
        tower = self.current_tower
        if tower is None or tower.degree != group.degree:
            tower = BlockTower(group.degree)
        result = self.solve(action.image, source_window, target_window, tower=tower.restrict(window))
        result = action.preimage_coset(result)
        return None if result is None else result.multiply_right(coset.representative)
