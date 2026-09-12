"""Permutations acting on the right, on the domain ``range(n)``.

``compose(p, q)[i] == q[p[i]]``: apply p first, then q. All group,
coset, and string operations in this package use this convention.
"""

from collections.abc import Iterable, Sequence

Permutation = tuple[int, ...]


def identity(n: int) -> Permutation:
    if not isinstance(n, int) or isinstance(n, bool) or n < 0:
        raise ValueError("degree must be a nonnegative integer")
    return tuple(range(n))


def validate(p: Sequence[int], n: int | None = None) -> Permutation:
    p = tuple(p)
    if n is None:
        n = len(p)
    if len(p) != n or any(type(i) is not int for i in p) or set(p) != set(range(n)):
        raise ValueError(f"not a permutation of range({n})")
    return p


def compose(p: Permutation, q: Permutation) -> Permutation:
    if len(p) != len(q):
        raise ValueError("permutations have different degrees")
    return tuple(q[j] for j in p)


def inverse(p: Permutation) -> Permutation:
    result = [0] * len(p)
    for i, j in enumerate(p):
        result[j] = i
    return tuple(result)


def from_cycles(n: int, *cycles: Iterable[int]) -> Permutation:
    """Compose the listed cycles from left to right; cycles may overlap."""
    result = identity(n)
    for cycle in cycles:
        cycle = tuple(cycle)
        if len(set(cycle)) != len(cycle) or any(type(i) is not int or not 0 <= i < n for i in cycle):
            raise ValueError("invalid cycle")
        p = list(range(n))
        for i, j in zip(cycle, cycle[1:] + cycle[:1]):
            p[i] = j
        result = compose(result, tuple(p))
    return result


def sign(p: Permutation) -> int:
    visited: set[int] = set()
    result = 1
    for i in range(len(p)):
        if i in visited:
            continue
        length = 0
        while i not in visited:
            visited.add(i)
            length += 1
            i = p[i]
        if length % 2 == 0:
            result = -result
    return result


def pullback(values: Sequence, p: Permutation) -> tuple:
    """Return ``values[p[i]]``; this is the string ``values^(p^-1)``."""
    if len(values) != len(p):
        raise ValueError("string and permutation have different degrees")
    return tuple(values[j] for j in p)

