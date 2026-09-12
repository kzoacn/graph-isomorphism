"""Invariant block towers carried through Luks and Babai recursion."""

from dataclasses import dataclass
from collections.abc import Iterable

from .group import PermutationGroup, _points
from .permutation import identity

Partition = tuple[tuple[int, ...], ...]


def _partition(degree: int, blocks: Iterable[Iterable[int]]) -> Partition:
    blocks = tuple(sorted(tuple(sorted(block)) for block in blocks))
    points = tuple(v for block in blocks for v in block)
    _points(points, degree)
    if len(points) != degree or any(not block for block in blocks):
        raise ValueError("blocks must partition the domain")
    return blocks


def refines(fine: Partition, coarse: Partition) -> bool:
    location = {v: i for i, block in enumerate(coarse) for v in block}
    return all(len({location.get(v, -1) for v in block}) == 1 and all(v in location for v in block) for block in fine)


@dataclass(frozen=True)
class BlockTower:
    degree: int
    levels: tuple[Partition, ...] = ()

    def __post_init__(self):
        identity(self.degree)
        singleton = tuple((v,) for v in range(self.degree))
        levels = [singleton]
        previous = singleton
        for raw in self.levels:
            level = _partition(self.degree, raw)
            if not refines(previous, level):
                raise ValueError("block tower levels must be successive coarsenings")
            previous = level
            # The universal partition conveys no quotient information and
            # must not prevent extending the tower after an orbit descent.
            if level != levels[-1] and len(level) > 1:
                levels.append(level)
        object.__setattr__(self, "levels", tuple(levels))

    @property
    def top(self) -> Partition:
        return self.levels[-1]

    def extend(self, blocks: Iterable[Iterable[int]]) -> "BlockTower":
        return BlockTower(self.degree, self.levels + (_partition(self.degree, blocks),))

    def restrict(self, window: Iterable[int]) -> "BlockTower":
        window = _points(window, self.degree)
        positions = {v: i for i, v in enumerate(window)}
        levels = []
        for level in self.levels:
            blocks = [tuple(positions[v] for v in block if v in positions) for block in level]
            levels.append(tuple(block for block in blocks if block))
        return BlockTower(len(window), tuple(levels))

    def validate_for(self, group: PermutationGroup) -> None:
        if group.degree != self.degree:
            raise ValueError("block tower has the wrong degree")
        for level in self.levels:
            expected = {frozenset(block) for block in level}
            for g in group.generators:
                if {frozenset(g[v] for v in block) for block in level} != expected:
                    raise ValueError("block tower is not invariant under this group")

