"""The corrected coherent Split-or-Johnson reduction.

This is Helfgott Proposition 5.8, including Babai's 2017 repair of the
primitive-right-side case. That case makes *one* reduction with right
side less than half as large, rather than a forked recursive call.
"""

from collections import Counter
from collections.abc import Callable, Iterable
from dataclasses import dataclass

from .configuration import Configuration, is_clique, is_coherent
from .group import _points


@dataclass(frozen=True)
class BipartiteGraph:
    # Left vertices retain the input configuration's labels. Right vertices
    # can be contracted blocks, hence are represented by disjoint fibers.
    left: tuple[int, ...]
    right: tuple[frozenset[int], ...]
    neighborhoods: tuple[frozenset[int], ...]  # indices into right, one per left

    def __post_init__(self):
        if len(self.left) != len(self.neighborhoods) or len(set(self.left)) != len(self.left):
            raise ValueError("invalid bipartite left domain")
        right_vertices = [v for fiber in self.right for v in fiber]
        if any(not fiber for fiber in self.right) or len(set(right_vertices)) != len(right_vertices) or set(right_vertices).intersection(self.left):
            raise ValueError("right fibers must be nonempty, disjoint, and disjoint from the left domain")
        if any(any(type(i) is not int or not 0 <= i < len(self.right) for i in neighbors) for neighbors in self.neighborhoods):
            raise ValueError("invalid neighbor index")

    def twin_classes(self) -> tuple[tuple[int, ...], ...]:
        classes = {}
        for vertex, neighbors in zip(self.left, self.neighborhoods):
            classes.setdefault(neighbors, []).append(vertex)
        return tuple(tuple(vertices) for vertices in classes.values())

    def is_nontrivial_semiregular(self) -> bool:
        left_degrees = set(map(len, self.neighborhoods))
        right_degrees = {sum(i in neighbors for neighbors in self.neighborhoods) for i in range(len(self.right))}
        return len(left_degrees) == len(right_degrees) == 1 and 0 < next(iter(left_degrees)) < len(self.right)


@dataclass(frozen=True)
class CoherentSplit:
    # Distinct entries are distinct color classes. Blocks within one entry
    # can be permuted. Labels always refer to the input left domain.
    parts: tuple[tuple[tuple[int, ...], ...], ...]
    bipartite: BipartiteGraph | None
    fixed: tuple[int, ...]
    case: str


def homogeneous_partition(configuration: Configuration, vertices: Iterable[int]) -> tuple[tuple[int, ...], ...] | None:
    """Components of the first disconnected off-diagonal relation, if any."""
    vertices = tuple(vertices)
    colors = sorted({configuration.color((a, b)) for a in vertices for b in vertices if a != b})
    for color in colors:
        remaining = set(vertices)
        parts = []
        while remaining:
            seed = min(remaining)
            component = {seed}
            pending = [seed]
            remaining.remove(seed)
            while pending:
                a = pending.pop()
                neighbors = {b for b in remaining if configuration.color((a, b)) == color or configuration.color((b, a)) == color}
                pending.extend(neighbors)
                component.update(neighbors)
                remaining.difference_update(neighbors)
            parts.append(tuple(sorted(component)))
        if len(parts) > 1:
            if len({len(part) for part in parts}) != 1 or len(parts[0]) < 2:
                raise ValueError("input is not a homogeneous coherent configuration")
            return tuple(parts)
    return None


def _bipartite(configuration: Configuration, left: tuple[int, ...], fibers: tuple[frozenset[int], ...], color: int) -> BipartiteGraph:
    return BipartiteGraph(left, fibers, tuple(frozenset(i for i, fiber in enumerate(fibers) if any(configuration.color((v, w)) == color for w in fiber)) for v in left))


def coherent_split(configuration: Configuration, left: Iterable[int], right: Iterable[int], *, choose: Callable[[tuple[int, ...]], int] = min) -> CoherentSplit:
    """Reduce a coherent two-fiber configuration, exposing any chosen vertex.

    Input: |left| > |right|, nonconstant cross colors, and a nontrivial
    coherent configuration on right. Output: a colored half-partition of
    left, or a bipartite graph on W1,W2 with |W1| > |left|/2,
    |W2| <= |right|/2 and no twin class larger than |W1|/2.

    ``choose`` allows a caller to replay the one noncanonical choice on
    the second structure. It is never silently assumed to be invariant.
    """
    if configuration.arity != 2 or not is_coherent(configuration):
        raise ValueError("a binary coherent configuration is required")
    left, right = tuple(sorted(_points(left, configuration.degree))), tuple(sorted(_points(right, configuration.degree)))
    if len(left) <= len(right) or not right:
        raise ValueError("require |left| > |right| > 0")
    if {frozenset(left), frozenset(right)} != {frozenset(vertices) for _, vertices in configuration.vertex_classes()}:
        raise ValueError("left and right must be exactly the two vertex color classes")
    cross_colors = sorted({configuration.color((a, b)) for a in left for b in right})
    if len(cross_colors) < 2 or is_clique(configuration.induced(right)):
        raise ValueError("cross colors and the right configuration must be nontrivial")
    if is_clique(configuration.induced(left)):
        raise AssertionError("Fisher's inequality forbids a left clique")
    partition = homogeneous_partition(configuration, left)
    if partition is not None:
        return CoherentSplit((partition,), None, (), "left_imprimitive")

    def select(candidates):
        candidates = tuple(sorted(candidates))
        selected = choose(candidates)
        if selected not in candidates:
            raise ValueError("chosen point is not among the admissible choices")
        return selected

    degrees = Counter(configuration.color((v, right[0])) for v in left)
    violet = next((c for c in cross_colors if degrees[c] * 2 > len(left)), None)
    if violet is None:
        pivot = select(right)
        cells = {}
        for v in left:
            cells.setdefault(configuration.color((v, pivot)), []).append(v)
        return CoherentSplit(tuple((tuple(vertices),) for _, vertices in sorted(cells.items())), None, (pivot,), "balanced_cross_colors")

    right_partition = homogeneous_partition(configuration, right)
    if right_partition is not None:
        good_blocks = []
        for block in right_partition:
            graph = _bipartite(configuration, left, tuple(frozenset((w,)) for w in block), violet)
            if max(map(len, graph.twin_classes())) * 2 <= len(left):
                good_blocks.append(block)
        if good_blocks:
            pivot = select(v for block in good_blocks for v in block)
            block = next(block for block in good_blocks if pivot in block)
            graph = _bipartite(configuration, left, tuple(frozenset((w,)) for w in block), violet)
            return CoherentSplit((), graph, (pivot,), "right_block")
        green = next(c for c in cross_colors if c != violet)
        graph = _bipartite(configuration, left, tuple(map(frozenset, right_partition)), green)
        if not graph.is_nontrivial_semiregular() or max(map(len, graph.twin_classes())) * 2 > len(left):
            raise AssertionError("contracted graph fails the coherent reduction conditions")
        return CoherentSplit((), graph, (), "contracted_right_blocks")

    # Babai's corrected primitive-right case. The two neighborhoods use
    # the same pivot. There is only one resulting bipartite instance.
    pivot = select(right)
    new_left = tuple(v for v in left if configuration.color((v, pivot)) == violet)
    blue_candidates = sorted({configuration.color((w, pivot)) for w in right if w != pivot})
    blue = next((c for c in blue_candidates if 0 < 2 * sum(configuration.color((w, pivot)) == c for w in right) < len(right)), None)
    if blue is None:
        raise AssertionError("nontrivial right configuration has no small relation")
    new_right = tuple(w for w in right if configuration.color((w, pivot)) == blue)
    graph = _bipartite(configuration, new_left, tuple(frozenset((w,)) for w in new_right), violet)
    if not graph.is_nontrivial_semiregular() or max(map(len, graph.twin_classes())) * 2 > len(new_left):
        raise AssertionError("corrected primitive reduction conditions violated")
    remainder = tuple(v for v in left if v not in new_left)
    parts = ((remainder,),) if remainder else ()
    return CoherentSplit(parts, graph, (pivot,), "primitive_right_fix")

