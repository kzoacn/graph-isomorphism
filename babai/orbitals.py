"""Orbital configurations and the bounded transitivity reduction.

Orbital color numbers are invariant under the given group, but are not
canonical under arbitrary conjugation of groups. For two different
certificate groups the main algorithm compares individual orbital graphs
against every possible matching orbital, instead of comparing these IDs.
"""

from collections import deque
from math import ceil, log

from .action import restriction_action
from .configuration import Configuration
from .group import PermutationGroup
from .local_certificates import contains_alternating


def orbital_configuration(group: PermutationGroup) -> Configuration:
    n = group.degree
    colors = [-1] * (n * n)
    color = 0
    for start in range(n * n):
        if colors[start] != -1:
            continue
        colors[start] = color
        pending = deque([divmod(start, n)])
        while pending:
            a, b = pending.popleft()
            for g in group.generators:
                i, j = g[a], g[b]
                index = i * n + j
                if colors[index] == -1:
                    colors[index] = color
                    pending.append((i, j))
        color += 1
    return Configuration(n, 2, tuple(colors))


def orbital_graphs(group: PermutationGroup) -> tuple[Configuration, ...]:
    configuration = orbital_configuration(group)
    n = group.degree
    off_diagonal = sorted({configuration.color((a, b)) for a in range(n) for b in range(n) if a != b})
    return tuple(Configuration(n, 2, tuple(0 if a == b else 1 if configuration.color((a, b)) == color else 2 for a in range(n) for b in range(n))) for color in off_diagonal)


def is_doubly_transitive(group: PermutationGroup) -> bool:
    if group.degree < 2 or len(group.orbits()) != 1:
        return False
    stabilizer = group.pointwise_stabilizer((0,))
    return len(stabilizer.orbits(range(1, group.degree))) == 1


def break_double_transitivity(group: PermutationGroup) -> tuple[int, ...]:
    """Find a prefix after which the remaining transitive action is not 2-transitive."""
    if group.degree < 3 or len(group.orbits()) != 1 or contains_alternating(group):
        raise ValueError("require a transitive non-giant group of degree at least three")
    fixed = []
    current = group
    while True:
        remaining = tuple(v for v in range(group.degree) if v not in fixed)
        image = restriction_action(current, remaining).image
        if not is_doubly_transitive(image):
            return tuple(fixed)
        fixed.append(remaining[0])
        if len(fixed) > ceil(3 * log(group.degree)):
            raise AssertionError("Wielandt's transitivity bound was violated")
        current = current.pointwise_stabilizer((fixed[-1],))

