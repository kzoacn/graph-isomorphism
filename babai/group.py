"""Deterministic permutation group algorithms, without group enumeration.

The stabilizer chain uses Schreier generators and deterministic sifting.
It deliberately favors an inspectable polynomial-time construction over
randomized, faster variants of Schreier--Sims. Enumeration is a separate
operation and always requires an explicit upper bound.
"""

from collections import deque
from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass
from functools import cached_property
from itertools import product
from math import prod

from .permutation import Permutation, compose, from_cycles, identity, inverse, validate


class EnumerationLimitExceeded(ValueError):
    """An explicitly bounded enumeration would exceed its bound."""


def _points(points: Iterable[int], n: int) -> tuple[int, ...]:
    result = tuple(points)
    if len(set(result)) != len(result) or any(type(p) is not int or not 0 <= p < n for p in result):
        raise ValueError("points must be distinct elements of the permutation domain")
    return result


def _orbit_transversal(n: int, generators: Sequence[Permutation], point: int) -> dict[int, Permutation]:
    representatives = {point: identity(n)}
    pending = deque([point])
    while pending:
        current = pending.popleft()
        representative = representatives[current]
        for g in generators:
            target = g[current]
            if target not in representatives:
                representatives[target] = compose(representative, g)
                pending.append(target)
    return representatives


@dataclass(frozen=True)
class StabilizerChain:
    base: tuple[int, ...]
    transversals: tuple[dict[int, Permutation], ...]
    inverse_transversals: tuple[dict[int, Permutation], ...]

    def sift(self, g: Permutation, start: int = 0) -> Permutation:
        for i in range(start, len(self.base)):
            representative_inverse = self.inverse_transversals[i].get(g[self.base[i]])
            if representative_inverse is None:
                return g
            g = compose(g, representative_inverse)
        return g

    @property
    def order(self) -> int:
        return prod(map(len, self.transversals))

    def generators_from(self, level: int = 0) -> tuple[Permutation, ...]:
        e = identity(len(self.base))
        return tuple(dict.fromkeys(g for table in self.transversals[level:] for g in table.values() if g != e))


def _schreier_sims(n: int, generators: Sequence[Permutation], base: tuple[int, ...]) -> StabilizerChain:
    """Complete a strong generating set by repeatedly closing Schreier tests.

    Each new residual expands at least one basic orbit. Basic orbits grow
    monotonically and there are at most n(n-1)/2 such expansions. All loops
    between expansions range over points and a polynomial number of generators.
    In particular this routine never loops over the elements of G.
    """
    e = identity(n)
    strong = list(dict.fromkeys(g for g in generators if g != e))
    while True:
        levels: list[tuple[Permutation, ...]] = []
        active = tuple(strong)
        tables = []
        for point in base:
            levels.append(active)
            tables.append(_orbit_transversal(n, active, point))
            active = tuple(g for g in active if g[point] == point)
        chain = StabilizerChain(base, tuple(tables), tuple({p: inverse(t) for p, t in table.items()} for table in tables))
        residual = e
        for i in range(n - 1, -1, -1):
            for point, representative in tables[i].items():
                for g in levels[i]:
                    # t(point) * g * t(g(point))^-1 fixes this base point.
                    schreier = compose(compose(representative, g), chain.inverse_transversals[i][g[point]])
                    residual = chain.sift(schreier, i + 1)
                    if residual != e:
                        break
                if residual != e:
                    break
            if residual != e:
                break
        if residual == e:
            return chain
        strong.append(residual)


@dataclass(frozen=True, eq=False)
class PermutationGroup:
    """A finite permutation group specified by generators, including its degree."""

    degree: int
    generators: tuple[Permutation, ...]

    def __init__(self, degree: int, generators: Iterable[Sequence[int]] = ()):
        e = identity(degree)
        object.__setattr__(self, "degree", degree)
        object.__setattr__(self, "generators", tuple(dict.fromkeys(p for g in generators if (p := validate(g, degree)) != e)))

    def __repr__(self) -> str:
        return f"PermutationGroup(degree={self.degree}, generators={len(self.generators)})"

    @cached_property
    def chain(self) -> StabilizerChain:
        return _schreier_sims(self.degree, self.generators, tuple(range(self.degree)))

    def chain_with_prefix(self, points: Iterable[int]) -> StabilizerChain:
        prefix = _points(points, self.degree)
        selected = set(prefix)
        base = prefix + tuple(p for p in range(self.degree) if p not in selected)
        if base == tuple(range(self.degree)):
            return self.chain
        return _schreier_sims(self.degree, self.generators, base)

    @property
    def order(self) -> int:
        return self.chain.order

    def contains(self, g: Sequence[int]) -> bool:
        g = validate(g, self.degree)
        return self.chain.sift(g) == identity(self.degree)

    def is_subgroup_of(self, other: "PermutationGroup") -> bool:
        return self.degree == other.degree and all(other.contains(g) for g in self.generators)

    def same_group(self, other: "PermutationGroup") -> bool:
        return self.is_subgroup_of(other) and self.order == other.order

    def compact(self) -> "PermutationGroup":
        """An equivalent group with at most n(n-1)/2 generators."""
        return PermutationGroup(self.degree, self.chain.generators_from())

    def orbit_transversal(self, point: int) -> dict[int, Permutation]:
        _points((point,), self.degree)
        return _orbit_transversal(self.degree, self.generators, point)

    def orbits(self, domain: Iterable[int] | None = None) -> tuple[tuple[int, ...], ...]:
        points = set(range(self.degree) if domain is None else _points(domain, self.degree))
        if any(g[p] not in points for g in self.generators for p in points):
            raise ValueError("domain is not invariant under the group")
        result = []
        while points:
            first = min(points)
            orbit = {first}
            pending = [first]
            while pending:
                p = pending.pop()
                for g in self.generators:
                    if g[p] not in orbit:
                        orbit.add(g[p])
                        pending.append(g[p])
            points.difference_update(orbit)
            result.append(tuple(sorted(orbit)))
        return tuple(result)

    def pointwise_stabilizer(self, points: Iterable[int]) -> "PermutationGroup":
        points = _points(points, self.degree)
        if not points:
            return self
        chain = self.chain_with_prefix(points)
        return PermutationGroup(self.degree, chain.generators_from(len(points)))

    def transporter(self, source: Iterable[int], target: Iterable[int]) -> "Coset | None":
        """All g in G mapping an ordered tuple source to target, as H * r."""
        source, target = _points(source, self.degree), _points(target, self.degree)
        if len(source) != len(target):
            raise ValueError("source and target have different lengths")
        chain = self.chain_with_prefix(source)
        # After matching i points, candidates are G_(source[:i]) * representative.
        representative = identity(self.degree)
        for i, (a, b) in enumerate(zip(source, target)):
            preimage = inverse(representative)[b]
            t = chain.transversals[i].get(preimage)
            if t is None:
                return None
            representative = compose(t, representative)
        return Coset(PermutationGroup(self.degree, chain.generators_from(len(source))), representative)

    def elements(self, *, max_order: int) -> Iterator[Permutation]:
        """Enumerate only if the caller's explicit bound permits the entire group."""
        if not isinstance(max_order, int) or max_order < 0:
            raise ValueError("max_order must be a nonnegative integer")
        if self.order > max_order:
            raise EnumerationLimitExceeded(f"group order {self.order} exceeds enumeration bound {max_order}")
        e = identity(self.degree)
        for representatives in product(*(tuple(table.values()) for table in self.chain.transversals)):
            result = e
            for r in reversed(representatives):
                result = compose(result, r)
            yield result

    def coset_representatives(self, subgroup: "PermutationGroup", *, max_index: int) -> tuple[Permutation, ...]:
        """Representatives of right cosets H*r, with a mandatory index bound."""
        if not isinstance(max_index, int) or max_index < 0:
            raise ValueError("max_index must be a nonnegative integer")
        if not subgroup.is_subgroup_of(self):
            raise ValueError("argument is not a subgroup")
        index = self.order // subgroup.order
        if index > max_index:
            raise EnumerationLimitExceeded(f"subgroup index {index} exceeds enumeration bound {max_index}")
        representatives = [identity(self.degree)]
        inverse_representatives = [identity(self.degree)]
        pos = 0
        while pos < len(representatives):
            r = representatives[pos]
            pos += 1
            for g in self.generators:
                candidate = compose(r, g)
                if not any(subgroup.contains(compose(candidate, ri)) for ri in inverse_representatives):
                    representatives.append(candidate)
                    inverse_representatives.append(inverse(candidate))
        if len(representatives) != index:
            raise AssertionError("incomplete coset traversal")
        return tuple(representatives)


@dataclass(frozen=True)
class Coset:
    """The right coset ``subgroup * representative``; an empty coset is None."""

    subgroup: PermutationGroup
    representative: Permutation

    def __post_init__(self):
        object.__setattr__(self, "representative", validate(self.representative, self.subgroup.degree))

    @property
    def order(self) -> int:
        return self.subgroup.order

    def contains(self, g: Sequence[int]) -> bool:
        return self.subgroup.contains(compose(validate(g, self.subgroup.degree), inverse(self.representative)))

    def multiply_right(self, g: Permutation) -> "Coset":
        return Coset(self.subgroup, compose(self.representative, g))


def symmetric_group(n: int) -> PermutationGroup:
    identity(n)
    if n < 2:
        return PermutationGroup(n)
    return PermutationGroup(n, (from_cycles(n, (0, 1)), from_cycles(n, range(n))))


def alternating_group(n: int) -> PermutationGroup:
    identity(n)
    return PermutationGroup(n, (from_cycles(n, (0, 1, i)) for i in range(2, n)))


def merge_isomorphism_cosets(cosets: Iterable[Coset | None]) -> Coset | None:
    """Join branches whose union is known to be one isomorphism coset.

    This is *not* a union operation for arbitrary cosets: the caller must
    establish that their union is a coset, e.g. by Luks's weak reduction.
    """
    return coset_hull(cosets)


def coset_hull(cosets: Iterable[Coset | None]) -> Coset | None:
    """The least right coset containing the supplied nonempty cosets.

    Unlike set union this operation generally adds elements. For image
    generator fibers in TopAction the algebraic argument proves that
    their hull is precisely the desired isomorphism coset.
    """
    first = None
    generators = []
    for coset in cosets:
        if coset is None:
            continue
        if first is None:
            first = coset
        if coset.subgroup.degree != first.subgroup.degree:
            raise ValueError("cosets have different degrees")
        generators.extend(coset.subgroup.generators)
        generators.append(compose(coset.representative, inverse(first.representative)))
    if first is None:
        return None
    group = PermutationGroup(first.subgroup.degree, generators).compact()
    return Coset(group, first.representative)
