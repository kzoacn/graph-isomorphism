"""Relational configurations and correlated k-ary Weisfeiler--Leman.

These are Babai/Helfgott's F2 and F3 functors, not a GI decision procedure.
Comparisons must refine both structures together to share color meanings.
"""

from collections import Counter
from dataclasses import dataclass
from itertools import product
from collections.abc import Iterable

from .group import _points
from .permutation import Permutation, validate


def tuple_index(points: tuple[int, ...], n: int) -> int:
    result = 0
    for p in points:
        result = result * n + p
    return result


def _equality_pattern(points: tuple[int, ...]) -> tuple[int, ...]:
    first = {}
    return tuple(first.setdefault(p, i) for i, p in enumerate(points))


@dataclass(frozen=True)
class Configuration:
    degree: int
    arity: int
    colors: tuple[int, ...]

    def __post_init__(self):
        if type(self.degree) is not int or self.degree < 0 or type(self.arity) is not int or self.arity < 1:
            raise ValueError("require degree >= 0 and arity >= 1")
        colors = tuple(self.colors)
        if len(colors) != self.degree**self.arity or any(type(c) is not int or c < 0 for c in colors):
            raise ValueError("provide one nonnegative integer color per tuple")
        object.__setattr__(self, "colors", colors)

    def color(self, points: Iterable[int]) -> int:
        points = tuple(points)
        if len(points) != self.arity or any(type(p) is not int or not 0 <= p < self.degree for p in points):
            raise ValueError("invalid tuple")
        return self.colors[tuple_index(points, self.degree)]

    def skeleton(self, arity: int) -> "Configuration":
        if not 1 <= arity <= self.arity:
            raise ValueError("invalid skeleton arity")
        return Configuration(self.degree, arity, tuple(self.color(t + (t[-1],) * (self.arity - arity)) for t in product(range(self.degree), repeat=arity)))

    def condition(self, prefix: Iterable[int]) -> "Configuration":
        """Fix coordinates, reducing arity; no new refinement is necessary."""
        prefix = tuple(prefix)
        if len(prefix) >= self.arity or any(type(p) is not int or not 0 <= p < self.degree for p in prefix):
            raise ValueError("invalid conditioning tuple")
        k = self.arity - len(prefix)
        return Configuration(self.degree, k, tuple(self.color(prefix + t) for t in product(range(self.degree), repeat=k)))

    def induced(self, vertices: Iterable[int]) -> "Configuration":
        vertices = _points(vertices, self.degree)
        return Configuration(len(vertices), self.arity, tuple(self.color(t) for t in product(vertices, repeat=self.arity)))

    def relabel(self, permutation: Permutation) -> "Configuration":
        permutation = validate(permutation, self.degree)
        colors = [0] * len(self.colors)
        for points, color in zip(product(range(self.degree), repeat=self.arity), self.colors):
            colors[tuple_index(tuple(permutation[p] for p in points), self.degree)] = color
        return Configuration(self.degree, self.arity, tuple(colors))

    def vertex_classes(self) -> tuple[tuple[int, tuple[int, ...]], ...]:
        classes: dict[int, list[int]] = {}
        for v in range(self.degree):
            classes.setdefault(self.color((v,) * self.arity), []).append(v)
        return tuple((c, tuple(vertices)) for c, vertices in sorted(classes.items()))

    def twin_classes(self) -> tuple[tuple[int, ...], ...]:
        """Classes of vertices whose transposition is an automorphism.

        This exact O(n^(k+2) k) test is used to verify Design Lemma inputs.
        It is unrelated to equality of WL vertex colors, a weaker condition.
        """
        remaining = set(range(self.degree))
        classes = []
        tuples = tuple(product(range(self.degree), repeat=self.arity))
        while remaining:
            a = min(remaining)
            twins = [a]
            for b in sorted(remaining - {a}):
                if all(c == self.colors[tuple_index(tuple(b if v == a else a if v == b else v for v in t), self.degree)] for t, c in zip(tuples, self.colors)):
                    twins.append(b)
            remaining.difference_update(twins)
            classes.append(tuple(twins))
        return tuple(classes)


def _index_joint(signature_lists: tuple[tuple, ...]) -> tuple[tuple[tuple[int, ...], ...], tuple]:
    # Sorted full signatures, never hash digests: no hash collision can equate colors.
    palette = tuple(sorted({signature for signatures in signature_lists for signature in signatures}))
    index = {signature: i for i, signature in enumerate(palette)}
    return tuple(tuple(index[s] for s in signatures) for signatures in signature_lists), palette


def _configuration_signatures(configuration: Configuration) -> tuple:
    n, k = configuration.degree, configuration.arity
    maps = tuple(product(range(k), repeat=k))
    return tuple((_equality_pattern(t), tuple(configuration.color(tuple(t[i] for i in mapping)) for mapping in maps)) for t in product(range(n), repeat=k))


def _wl_signatures(configuration: Configuration) -> tuple:
    n, k = configuration.degree, configuration.arity
    signatures = []
    powers = tuple(n**i for i in reversed(range(k)))
    for index, points in enumerate(product(range(n), repeat=k)):
        counts = Counter(tuple(configuration.colors[index + (z - points[j]) * powers[j]] for j in range(k)) for z in range(n))
        signatures.append((configuration.colors[index], tuple(sorted(counts.items()))))
    return tuple(signatures)


@dataclass(frozen=True)
class Refinement:
    configurations: tuple[Configuration, ...]
    rounds: int
    # Each palette records the meaning of new colors in terms of old colors.
    palettes: tuple[tuple, ...]


def refine_joint(*structures: Configuration) -> Refinement:
    """F2 followed by correlated k-WL, using one palette for all arguments.

    F2 records *all* coordinate maps, including noninjective maps. F3
    records the joint distribution over coordinate replacements at the
    *same* z. Independent marginal counts would be a weaker algorithm.

    The bound is n^O(k), including F2 when k <= n. For k=O(log n) this
    subroutine is quasipolynomial; it does not solve GI on its own.
    """
    if not structures:
        raise ValueError("provide at least one structure")
    if len({s.arity for s in structures}) != 1:
        raise ValueError("structures must have the same arity")
    colors, palette = _index_joint(tuple(_configuration_signatures(s) for s in structures))
    current = tuple(Configuration(s.degree, s.arity, c) for s, c in zip(structures, colors))
    palettes = [palette]
    color_count = len(palette)
    rounds = 0
    if not color_count:
        return Refinement(current, rounds, tuple(palettes))
    while True:
        colors, palette = _index_joint(tuple(_wl_signatures(s) for s in current))
        current = tuple(Configuration(s.degree, s.arity, c) for s, c in zip(current, colors))
        palettes.append(palette)
        rounds += 1
        if len(palette) == color_count:
            return Refinement(current, rounds, tuple(palettes))
        color_count = len(palette)


def is_coherent(configuration: Configuration) -> bool:
    """Verify equality-pattern, coordinate-map, and intersection-number axioms."""
    for signatures in (_configuration_signatures(configuration), _wl_signatures(configuration)):
        seen = {}
        for color, signature in zip(configuration.colors, signatures):
            if color in seen and seen[color] != signature:
                return False
            seen[color] = signature
    return True


def is_clique(configuration: Configuration) -> bool:
    if configuration.arity != 2:
        raise ValueError("clique test requires a binary configuration")
    n = configuration.degree
    diagonal = {configuration.color((i, i)) for i in range(n)}
    off_diagonal = {configuration.color((i, j)) for i in range(n) for j in range(n) if i != j}
    return len(diagonal) <= 1 and len(off_diagonal) <= 1 and not diagonal.intersection(off_diagonal)


def individualization_compatible(left: Configuration, right: Configuration, left_fixed: tuple[int, ...], right_fixed: tuple[int, ...]) -> bool:
    """Necessary isomorphism condition after a prescribed point mapping.

    Input color labels must already have shared meanings. This polynomial
    binary refinement is pruning *inside* the bounded Babai choice family;
    it does not introduce further individualizations or backtracking.
    """
    if left.degree != right.degree or left.arity < 2 or right.arity < 2:
        raise ValueError("require equal-degree configurations of arity at least two")
    n = left.degree
    left_fixed, right_fixed = _points(left_fixed, n), _points(right_fixed, n)
    if len(left_fixed) != len(right_fixed):
        raise ValueError("fixed tuples have different lengths")
    signatures = []
    for configuration, fixed in ((left, left_fixed), (right, right_fixed)):
        binary = configuration.skeleton(2)
        tags = {v: i + 1 for i, v in enumerate(fixed)}
        signatures.append(tuple((binary.color((a, b)), tags.get(a, 0), tags.get(b, 0)) for a in range(n) for b in range(n)))
    colors, _ = _index_joint(tuple(signatures))
    a, b = refine_joint(*(Configuration(n, 2, c) for c in colors)).configurations
    return Counter(a.colors) == Counter(b.colors)
