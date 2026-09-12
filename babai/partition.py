"""Explicit stabilizers and alignments of colored partitions.

These have known product/wreath-product structure. They must not be
confused with finding an arbitrary set stabilizer in an arbitrary group.
"""

from collections.abc import Iterable
from dataclasses import dataclass

from .action import Homomorphism
from .group import Coset, PermutationGroup, _points
from .permutation import compose, from_cycles, identity, sign


def even_subgroup(group: PermutationGroup) -> PermutationGroup:
    return Homomorphism(group, 2, (identity(2) if sign(g) == 1 else (1, 0) for g in group.generators)).kernel


def young_subgroup(degree: int, cells: Iterable[Iterable[int]], *, even: bool = False) -> PermutationGroup:
    """The product of symmetric groups on a partition's cells, intersected with A_n if requested."""
    cells = tuple(tuple(cell) for cell in cells)
    points = tuple(p for cell in cells for p in cell)
    _points(points, degree)
    if len(points) != degree:
        raise ValueError("cells must partition the domain")
    group = PermutationGroup(degree, (from_cycles(degree, (cell[0], p)) for cell in cells for p in cell[1:]))
    return even_subgroup(group) if even else group


@dataclass(frozen=True)
class ColoredPartition:
    degree: int
    # A color class consists of equally sized blocks; blocks of one color
    # may be permuted. Different colors may not be exchanged.
    classes: tuple[tuple[int, tuple[tuple[int, ...], ...]], ...]

    def __post_init__(self):
        identity(self.degree)
        classes = tuple(sorted((color, tuple(sorted(tuple(sorted(block)) for block in blocks))) for color, blocks in self.classes))
        if len({c for c, _ in classes}) != len(classes) or any(type(c) is not int or c < 0 for c, _ in classes):
            raise ValueError("colors must be distinct nonnegative integers")
        if any(not blocks or any(not block for block in blocks) or len({len(b) for b in blocks}) != 1 for _, blocks in classes):
            raise ValueError("a class must have nonempty blocks of equal size")
        points = tuple(p for _, blocks in classes for block in blocks for p in block)
        _points(points, self.degree)
        if len(points) != self.degree:
            raise ValueError("classes must partition the domain")
        object.__setattr__(self, "classes", classes)

    def stabilizer(self, *, even: bool = False) -> PermutationGroup:
        generators = []
        for _, blocks in self.classes:
            for block in blocks:
                generators.extend(from_cycles(self.degree, (block[0], p)) for p in block[1:])
            for block in blocks[1:]:
                generators.append(from_cycles(self.degree, *zip(blocks[0], block)))
        group = PermutationGroup(self.degree, generators)
        return even_subgroup(group) if even else group

    def align(self, target: "ColoredPartition", *, even: bool = False) -> Coset | None:
        if self.degree != target.degree:
            return None
        shape = lambda p: tuple((c, len(blocks), len(blocks[0])) for c, blocks in p.classes)
        if shape(self) != shape(target):
            return None
        representative = list(range(self.degree))
        for (_, source_blocks), (_, target_blocks) in zip(self.classes, target.classes):
            for source_block, target_block in zip(source_blocks, target_blocks):
                for a, b in zip(source_block, target_block):
                    representative[a] = b
        representative = tuple(representative)
        stabilizer = self.stabilizer()
        if even:
            if sign(representative) == -1:
                odd = next((g for g in stabilizer.generators if sign(g) == -1), None)
                if odd is None:
                    return None
                representative = compose(odd, representative)
            stabilizer = even_subgroup(stabilizer)
        return Coset(stabilizer, representative)

