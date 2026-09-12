"""Polynomial-time block systems using orbits of unordered pairs."""

from collections import deque

from .action import block_action
from .group import PermutationGroup
from collections.abc import Iterable


def _pair_components(group: PermutationGroup, a: int, b: int) -> tuple[tuple[int, ...], ...]:
    parent = list(range(group.degree))

    def root(p):
        while p != parent[p]:
            parent[p] = parent[parent[p]]
            p = parent[p]
        return p

    first = tuple(sorted((a, b)))
    pending = deque([first])
    seen = {first}
    while pending:
        a, b = pending.popleft()
        parent[root(a)] = root(b)
        for g in group.generators:
            pair = tuple(sorted((g[a], g[b])))
            if pair not in seen:
                seen.add(pair)
                pending.append(pair)
    parts: dict[int, list[int]] = {}
    for p in range(group.degree):
        parts.setdefault(root(p), []).append(p)
    return tuple(sorted(map(tuple, parts.values())))


def nontrivial_block_system(group: PermutationGroup) -> tuple[tuple[int, ...], ...] | None:
    """Find proper nonsingleton blocks of a transitive group, or certify primitivity."""
    if group.degree == 0 or len(group.orbits()) != 1:
        raise ValueError("block search requires a nonempty transitive action")
    for b in range(1, group.degree):
        blocks = _pair_components(group, 0, b)
        if len(blocks) > 1:
            return blocks
    return None


def maximal_block_system(group: PermutationGroup, initial: Iterable[Iterable[int]] | None = None) -> tuple[tuple[int, ...], ...]:
    """Return blocks with primitive quotient action (singletons when primitive).

    Maximal refers to block size, so the quotient has a minimal number of
    blocks greater than one. Choices depend on the ambient group, not strings.
    """
    if group.degree == 0 or len(group.orbits()) != 1:
        raise ValueError("block search requires a nonempty transitive action")
    blocks = tuple((p,) for p in range(group.degree)) if initial is None else tuple(tuple(block) for block in initial)
    quotient = group if initial is None else block_action(group, blocks).image
    if group.degree > 1 and len(blocks) < 2:
        raise ValueError("an initial proper block system is required")
    while True:
        coarser = nontrivial_block_system(quotient)
        if coarser is None:
            return blocks
        blocks = tuple(tuple(sorted(p for i in part for p in blocks[i])) for part in coarser)
        quotient = block_action(group, blocks).image
