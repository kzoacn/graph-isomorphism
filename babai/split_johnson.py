"""Split-or-Johnson with explicit, replayable noncanonical choices.

The recursion follows Helfgott sections 5.2 and 5.3, using the corrected
coherent reduction. The small-right cutoff is enlarged to O(log^2 n),
which preserves a quasipolynomial index bound but is not an optimization
of Babai's original exponent.
"""

from collections import Counter
from collections.abc import Iterable
from dataclasses import dataclass
from fractions import Fraction
from itertools import product
from math import comb

from .coherent_split import BipartiteGraph, coherent_split, homogeneous_partition
from .configuration import Configuration, is_clique, is_coherent, refine_joint
from .design import design_at, design_lemma
from .johnson import JohnsonModel, recognize_johnson
from .partition import ColoredPartition


class IncompatibleChoices(ValueError):
    """The proposed images of individualized points cannot reproduce this branch."""


@dataclass(frozen=True)
class Choice:
    kind: str
    vertices: tuple[int, ...]


@dataclass(frozen=True)
class EmbeddedJohnson:
    vertices: tuple[int, ...]
    model: JohnsonModel


@dataclass(frozen=True)
class SplitResult:
    domain: tuple[int, ...]
    parts: tuple[tuple[tuple[int, ...], ...], ...]
    johnson: EmbeddedJohnson | None
    choices: tuple[Choice, ...]
    threshold: Fraction

    @property
    def fixed_points(self) -> tuple[int, ...]:
        return tuple(dict.fromkeys(v for choice in self.choices for v in choice.vertices))

    def colored_partition(self) -> ColoredPartition:
        positions = {v: i for i, v in enumerate(self.domain)}
        return ColoredPartition(len(self.domain), tuple((c, tuple(tuple(positions[v] for v in block) for block in blocks)) for c, blocks in enumerate(self.parts)))


class _Choices:
    def __init__(self, replay: Iterable[Choice] | None):
        self.replay = None if replay is None else tuple(replay)
        self.recorded = []

    def _take(self, kind, fibers, default_vertices):
        if self.replay is None:
            vertices = tuple(default_vertices)
        else:
            if len(self.recorded) >= len(self.replay) or self.replay[len(self.recorded)].kind != kind:
                raise IncompatibleChoices("choice types disagree")
            vertices = self.replay[len(self.recorded)].vertices
        positions = []
        for v in vertices:
            candidates = [i for i, fiber in enumerate(fibers) if v in fiber]
            if len(candidates) != 1:
                raise IncompatibleChoices("chosen vertex is outside the eligible fibers")
            positions.append(candidates[0])
        if len(set(positions)) != len(positions):
            raise IncompatibleChoices("distinct chosen points must specify distinct fibers")
        self.recorded.append(Choice(kind, tuple(vertices)))
        return tuple(positions)

    def point(self, fibers):
        indices = self._take("point", fibers, (min(min(fiber) for fiber in fibers),))
        if len(indices) != 1:
            raise IncompatibleChoices("one point is required")
        return indices[0]

    def order_all(self, fibers):
        indices = self._take("all", fibers, sorted(min(fiber) for fiber in fibers))
        if len(indices) != len(fibers):
            raise IncompatibleChoices("every right fiber must be individualized")
        return indices

    def design(self, configuration, fibers):
        if self.replay is None:
            witness = design_lemma(configuration, Fraction(2, 3))
            vertices = tuple(min(fibers[i]) for i in witness.fixed)
        else:
            vertices = ()
        fixed = self._take("design", fibers, vertices)
        witness = design_at(configuration, fixed, Fraction(2, 3))
        if witness is None:
            raise IncompatibleChoices("the proposed Design Lemma tuple is not a witness")
        return witness

    def finish(self):
        if self.replay is not None and len(self.recorded) != len(self.replay):
            raise IncompatibleChoices("unused noncanonical choices")


def _small_right_limit(left_size: int) -> int:
    return 64 * max(1, (left_size - 1).bit_length()) ** 2


def _restrict(graph: BipartiteGraph, left=None, right=None) -> BipartiteGraph:
    left = tuple(range(len(graph.left))) if left is None else tuple(left)
    right = tuple(range(len(graph.right))) if right is None else tuple(right)
    positions = {j: i for i, j in enumerate(right)}
    return BipartiteGraph(tuple(graph.left[i] for i in left), tuple(graph.right[j] for j in right), tuple(frozenset(positions[j] for j in graph.neighborhoods[i] if j in positions) for i in left))


def _no_large_twins(graph: BipartiteGraph) -> bool:
    return all(3 * len(part) <= 2 * len(graph.left) for part in graph.twin_classes())


def _choose_right_half(graph: BipartiteGraph, subset: Iterable[int]) -> BipartiteGraph:
    subset = tuple(sorted(subset))
    complement = tuple(i for i in range(len(graph.right)) if i not in subset)
    for vertices in (subset, complement):
        if not vertices or len(vertices) == len(graph.right):
            continue
        restricted = _restrict(graph, right=vertices)
        if _no_large_twins(restricted):
            return restricted
    raise AssertionError("twin separation lemma failed")


def _balanced_union(classes, size):
    classes = tuple(classes)
    for _, vertices in classes:
        if 3 * len(vertices) > size:
            if 3 * len(vertices) > 2 * size:
                raise AssertionError("the coloring has a dominant class")
            return vertices
    selected = []
    for _, vertices in classes:
        selected.extend(vertices)
        if 3 * len(selected) > size:
            return tuple(selected)
    raise AssertionError("invalid right coloring")


def _hybrid_configuration(graph: BipartiteGraph, right_configuration: Configuration) -> Configuration:
    a, b = len(graph.left), len(graph.right)
    signatures = []
    for i, j in product(range(a + b), repeat=2):
        if i < a and j < a:
            signature = (0, int(i == j))
        elif i < a:
            signature = (1, int(j - a in graph.neighborhoods[i]))
        elif j < a:
            signature = (2, int(i - a in graph.neighborhoods[j]))
        else:
            signature = (3, right_configuration.color((i - a, j - a)))
        signatures.append(signature)
    palette = {s: i for i, s in enumerate(sorted(set(signatures)))}
    return refine_joint(Configuration(a + b, 2, tuple(palette[s] for s in signatures))).configurations[0]


class _Split:
    def __init__(self, domain: tuple[int, ...], threshold: Fraction, plan):
        self.domain = domain
        self.threshold = threshold
        self.choices = _Choices(plan)
        self.parts = []
        self.johnson = None

    def add_cell(self, vertices):
        vertices = tuple(sorted(vertices))
        if vertices:
            self.parts.append((vertices,))

    def embed(self, graph):
        rank = len(graph.neighborhoods[0])
        self.johnson = EmbeddedJohnson(graph.left, JohnsonModel(len(graph.right), rank, graph.neighborhoods, tuple(range(rank + 1))))
        self.add_cell(graph.left)

    def bipartite(self, graph: BipartiteGraph):
        """Iterate; every charged reduction shrinks the right side by a constant."""
        while True:
            left_size, right_size = len(graph.left), len(graph.right)
            if left_size <= self.threshold:
                self.add_cell(graph.left)
                return
            if not right_size < self.threshold or not _no_large_twins(graph):
                raise ValueError("bipartite Split-or-Johnson hypotheses are not satisfied")
            if right_size <= _small_right_limit(left_size):
                order = self.choices.order_all(graph.right)
                cells = {}
                for v, neighbors in zip(graph.left, graph.neighborhoods):
                    cells.setdefault(tuple(int(i in neighbors) for i in order), []).append(v)
                for _, vertices in sorted(cells.items()):
                    self.add_cell(vertices)
                return

            # Twin classes and degrees are canonical. A dominant class can
            # require more work only when every vertex in it has no twin.
            twin_classes = graph.twin_classes()
            positions = {v: i for i, v in enumerate(graph.left)}
            classes = {}
            for twins in twin_classes:
                degree = len(graph.neighborhoods[positions[twins[0]]])
                classes.setdefault((len(twins), degree), []).append(tuple(sorted(twins)))
            dominant = next((key for key, blocks in classes.items() if key[0] == 1 and sum(map(len, blocks)) > self.threshold), None)
            if dominant is None:
                self.parts.extend(tuple(sorted(blocks)) for _, blocks in sorted(classes.items()))
                return
            keep = {v for block in classes[dominant] for v in block}
            self.add_cell(v for v in graph.left if v not in keep)
            graph = _restrict(graph, left=(i for i, v in enumerate(graph.left) if v in keep))
            left_size = len(graph.left)
            rank = dominant[1]
            if 2 * rank > right_size:
                universe = frozenset(range(right_size))
                graph = BipartiteGraph(graph.left, graph.right, tuple(universe - neighbors for neighbors in graph.neighborhoods))
                rank = right_size - rank
            if not 2 <= rank <= right_size // 2:
                raise AssertionError("a large twin-free uniform hypergraph has invalid rank")
            if left_size == comb(right_size, rank):
                if 2 * rank == right_size:
                    # Balanced subsets have a complement involution. Its
                    # pairs are an equivariant nontrivial block partition.
                    index = {neighbors: v for v, neighbors in zip(graph.left, graph.neighborhoods)}
                    universe = frozenset(range(right_size))
                    blocks = {tuple(sorted((v, index[universe - neighbors]))) for v, neighbors in zip(graph.left, graph.neighborhoods)}
                    self.parts.append(tuple(sorted(blocks)))
                else:
                    self.embed(graph)
                return

            exponent, power = 0, 1
            while power < left_size:
                exponent += 1
                power *= right_size
            arity = min(rank, 6 * exponent)
            colors = []
            for points in product(range(right_size), repeat=arity):
                subset = frozenset(points)
                colors.append(0 if len(subset) != arity else 1 + sum(subset.issubset(neighbors) for neighbors in graph.neighborhoods))
            relation = Configuration(right_size, arity, tuple(colors))
            twins = relation.twin_classes()
            large_twins = next((part for part in twins if 2 * len(part) > right_size), None)
            if large_twins is not None:
                if len(large_twins) == right_size:
                    raise AssertionError("the design counting bound was violated")
                # This step is canonical, so a decrease by just one costs
                # polynomial work but introduces no further group index.
                graph = _choose_right_half(graph, large_twins)
                continue

            refined = refine_joint(relation).configurations[0]
            witness = self.choices.design(refined, graph.right)
            conditioned = refined.condition(witness.fixed)
            if witness.dominant is None:
                subset = _balanced_union(conditioned.vertex_classes(), right_size)
                graph = _choose_right_half(graph, subset)
                continue

            hybrid = _hybrid_configuration(graph, conditioned.skeleton(2))
            left_classes, right_classes = [], []
            for c, vertices in hybrid.vertex_classes():
                if vertices[0] < left_size:
                    left_classes.append((c, vertices))
                else:
                    right_classes.append((c, tuple(v - left_size for v in vertices)))
            dominant_left = next((vertices for _, vertices in left_classes if len(vertices) > self.threshold), None)
            if dominant_left is None:
                for _, vertices in left_classes:
                    self.add_cell(graph.left[i] for i in vertices)
                return
            dominant_right = next((vertices for _, vertices in right_classes if 3 * len(vertices) > 2 * right_size), None)
            if dominant_right is None:
                graph = _choose_right_half(graph, _balanced_union(right_classes, right_size))
                continue
            keep = set(dominant_left)
            self.add_cell(v for i, v in enumerate(graph.left) if i not in keep)
            if len({hybrid.color((i, left_size + j)) for i in dominant_left for j in dominant_right}) == 1:
                graph = _restrict(graph, left=dominant_left, right=(j for j in range(right_size) if j not in dominant_right))
                continue

            # Feed the two homogeneous color classes to the corrected
            # coherent reduction and lift its labels/fibers back afterward.
            vertices = dominant_left + tuple(left_size + j for j in dominant_right)
            local = hybrid.induced(vertices)
            a = len(dominant_left)
            local_left = tuple(graph.left[i] for i in dominant_left)
            local_right = tuple(graph.right[j] for j in dominant_right)

            def choose(candidates):
                fibers = tuple(local_right[i - a] for i in candidates)
                return candidates[self.choices.point(fibers)]

            result = coherent_split(local, range(a), range(a, len(vertices)), choose=choose)
            self.parts.extend(tuple(tuple(local_left[i] for i in block) for block in blocks) for blocks in result.parts)
            if result.bipartite is None:
                return
            reduced = result.bipartite
            graph = BipartiteGraph(tuple(local_left[i] for i in reduced.left), tuple(frozenset().union(*(local_right[i - a] for i in fiber)) for fiber in reduced.right), reduced.neighborhoods)

    def finish(self):
        self.choices.finish()
        result = SplitResult(self.domain, tuple(self.parts), self.johnson, tuple(self.choices.recorded), self.threshold)
        result.colored_partition()  # Check coverage, disjointness, and block sizes.
        johnson_vertices = set() if self.johnson is None else set(self.johnson.vertices)
        for blocks in result.parts:
            union = {v for block in blocks for v in block}
            if union == johnson_vertices:
                continue
            if any(len(block) > self.threshold for block in blocks) or (len(union) > self.threshold and len(blocks[0]) < 2):
                raise AssertionError("output is not an admissible alpha-partition")
        if self.johnson is not None and len(self.johnson.vertices) <= self.threshold:
            raise AssertionError("embedded Johnson scheme is not dominant")
        return result


def bipartite_split(graph: BipartiteGraph, beta: Fraction = Fraction(2, 3), *, plan: Iterable[Choice] | None = None) -> SplitResult:
    if not Fraction(2, 3) <= beta < 1:
        raise ValueError("require 2/3 <= beta < 1")
    worker = _Split(graph.left, beta * len(graph.left), plan)
    worker.bipartite(graph)
    return worker.finish()


def split_or_johnson(configuration: Configuration, alpha: Fraction = Fraction(2, 3), *, plan: Iterable[Choice] | None = None) -> SplitResult:
    if not Fraction(2, 3) <= alpha < 1:
        raise ValueError("require 2/3 <= alpha < 1")
    if configuration.arity != 2 or not is_coherent(configuration) or len(configuration.vertex_classes()) != 1 or is_clique(configuration):
        raise ValueError("a nontrivial homogeneous coherent configuration is required")
    n = configuration.degree
    worker = _Split(tuple(range(n)), alpha * n, plan)
    partition = homogeneous_partition(configuration, range(n))
    if partition is not None:
        worker.parts.append(partition)
        return worker.finish()
    model = recognize_johnson(configuration)
    if model is not None:
        worker.johnson = EmbeddedJohnson(tuple(range(n)), model)
        worker.add_cell(range(n))
        return worker.finish()
    pivot = worker.choices.point(tuple(frozenset((v,)) for v in range(n)))
    cells = {}
    for v in range(n):
        cells.setdefault(configuration.color((pivot, v)), []).append(v)
    dominant = next((tuple(vertices) for vertices in cells.values() if len(vertices) > worker.threshold), None)
    if dominant is None:
        for _, vertices in sorted(cells.items()):
            worker.add_cell(vertices)
        return worker.finish()
    worker.add_cell(v for v in range(n) if v not in dominant)
    for _, right in sorted(cells.items()):
        if right[0] in dominant:
            continue
        colors = sorted({configuration.color((v, w)) for v in dominant for w in right})
        for color in colors:
            graph = BipartiteGraph(dominant, tuple(frozenset((w,)) for w in right), tuple(frozenset(i for i, w in enumerate(right) if configuration.color((v, w)) == color) for v in dominant))
            if graph.is_nontrivial_semiregular():
                worker.bipartite(graph)
                return worker.finish()
    raise AssertionError("uniprimitive configuration has no nontrivial bipartite reduction")
