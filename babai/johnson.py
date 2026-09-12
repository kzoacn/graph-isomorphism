"""Constructive recognition of Johnson schemes, without clique enumeration.

For m > 2k, a Johnson graph edge's common neighbors split into two
cliques, of sizes m-k-1 and k-1. The larger clique, together with the edge,
is a star corresponding to a (k-1)-subset. Repeating this construction
recovers the underlying m points. Every recovered representation is
checked against every pair color before it is returned.
"""

from dataclasses import dataclass
from itertools import combinations
from math import comb

from .configuration import Configuration
from .permutation import Permutation, validate


@dataclass(frozen=True)
class JohnsonModel:
    atoms: int
    subset_size: int
    subsets: tuple[frozenset[int], ...]
    distance_colors: tuple[int, ...]

    def induce(self, permutation: Permutation) -> Permutation:
        permutation = validate(permutation, self.atoms)
        index = {s: i for i, s in enumerate(self.subsets)}
        return tuple(index[frozenset(permutation[p] for p in subset)] for subset in self.subsets)

    def recover(self, permutation: Permutation) -> Permutation | None:
        """Recover the atom permutation of a scheme automorphism, if it exists."""
        permutation = validate(permutation, len(self.subsets))
        fibers = tuple(frozenset(i for i, subset in enumerate(self.subsets) if p in subset) for p in range(self.atoms))
        index = {fiber: p for p, fiber in enumerate(fibers)}
        images = []
        for fiber in fibers:
            image = index.get(frozenset(permutation[i] for i in fiber))
            if image is None:
                return None
            images.append(image)
        images = tuple(images)
        return images if self.induce(images) == permutation else None


def johnson_configuration(atoms: int, subset_size: int) -> Configuration:
    if type(atoms) is not int or type(subset_size) is not int or subset_size < 2 or atoms <= 2 * subset_size:
        raise ValueError("a Johnson scheme requires m > 2k and k >= 2")
    subsets = tuple(map(frozenset, combinations(range(atoms), subset_size)))
    return Configuration(len(subsets), 2, tuple(subset_size - len(a.intersection(b)) for a in subsets for b in subsets))


def _components(vertices: frozenset[int], adjacency: tuple[frozenset[int], ...]) -> tuple[frozenset[int], ...]:
    remaining = set(vertices)
    parts = []
    while remaining:
        seed = min(remaining)
        pending = [seed]
        part = {seed}
        remaining.remove(seed)
        while pending:
            v = pending.pop()
            neighbors = adjacency[v].intersection(remaining)
            pending.extend(neighbors)
            part.update(neighbors)
            remaining.difference_update(neighbors)
        parts.append(frozenset(part))
    return tuple(parts)


def _descend(adjacency: tuple[frozenset[int], ...], atoms: int, rank: int):
    n = len(adjacency)
    if n != comb(atoms, rank) or any(len(neighbors) != rank * (atoms - rank) for neighbors in adjacency):
        return None
    stars = set()
    for a in range(n):
        for b in adjacency[a]:
            if b <= a:
                continue
            common = adjacency[a].intersection(adjacency[b])
            if len(common) != atoms - 2:
                return None
            components = _components(common, adjacency)
            if len(components) != 2 or sorted(map(len, components)) != [rank - 1, atoms - rank - 1]:
                return None
            if any(any(not (part - {v}).issubset(adjacency[v]) for v in part) for part in components):
                return None
            stars.add(max(components, key=len).union((a, b)))
            if len(stars) > comb(atoms, rank - 1):
                return None
    if len(stars) != comb(atoms, rank - 1):
        return None
    stars = tuple(sorted(stars, key=lambda s: tuple(sorted(s))))
    incidence = tuple(tuple(i for i, star in enumerate(stars) if v in star) for v in range(n))
    if any(len(indices) != rank for indices in incidence):
        return None
    new_adjacency = [set() for _ in stars]
    for indices in incidence:
        for a, b in combinations(indices, 2):
            new_adjacency[a].add(b)
            new_adjacency[b].add(a)
    return incidence, tuple(map(frozenset, new_adjacency))


def _from_relation(configuration: Configuration, atoms: int, rank: int, color: int) -> JohnsonModel | None:
    n = configuration.degree
    adjacency = tuple(frozenset(b for b in range(n) if a != b and configuration.color((a, b)) == color) for a in range(n))
    levels = []
    for r in range(rank, 1, -1):
        result = _descend(adjacency, atoms, r)
        if result is None:
            return None
        incidence, adjacency = result
        levels.append(incidence)
    if len(adjacency) != atoms or any(len(neighbors) != atoms - 1 for neighbors in adjacency):
        return None
    subsets = tuple(frozenset((p,)) for p in range(atoms))
    for incidence in reversed(levels):
        subsets = tuple(frozenset().union(*(subsets[i] for i in indices)) for indices in incidence)
    if len(set(subsets)) != n or any(len(s) != rank for s in subsets):
        return None
    distance_colors = {}
    for a in range(n):
        for b in range(n):
            distance = rank - len(subsets[a].intersection(subsets[b]))
            c = configuration.color((a, b))
            if distance in distance_colors and distance_colors[distance] != c:
                return None
            distance_colors[distance] = c
    if len(distance_colors) != rank + 1 or len(set(distance_colors.values())) != rank + 1:
        return None
    return JohnsonModel(atoms, rank, subsets, tuple(distance_colors[d] for d in range(rank + 1)))


def recognize_johnson(configuration: Configuration) -> JohnsonModel | None:
    """Recognize a full binary Johnson scheme J(m,k), with m > 2k >= 4.

    The number of relation colors gives k. Numerical intersection data
    alone is not accepted as a certificate: the vertex-to-subset map is
    constructed and its entire coloring is verified. Time is polynomial
    in the number of configuration vertices.
    """
    if configuration.arity != 2:
        raise ValueError("Johnson recognition requires a binary configuration")
    n = configuration.degree
    rank = len(set(configuration.colors)) - 1
    if rank < 2 or n < comb(2 * rank + 1, rank):
        return None
    diagonal = {configuration.color((i, i)) for i in range(n)}
    if len(diagonal) != 1:
        return None
    if any(configuration.color((a, b)) != configuration.color((b, a)) for a in range(n) for b in range(a + 1, n)):
        return None
    low, high = 2 * rank + 1, n
    while low < high:
        middle = (low + high) // 2
        if comb(middle, rank) < n:
            low = middle + 1
        else:
            high = middle
    atoms = low
    if comb(atoms, rank) != n:
        return None
    for color in sorted(set(configuration.colors) - diagonal):
        if sum(configuration.color((0, v)) == color for v in range(1, n)) != rank * (atoms - rank):
            continue
        model = _from_relation(configuration, atoms, rank, color)
        if model is not None:
            return model
    return None

