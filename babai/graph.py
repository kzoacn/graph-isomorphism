"""Simple undirected graphs and the exact GI-to-string-isomorphism reduction."""

from dataclasses import dataclass
from itertools import combinations
from collections.abc import Iterable

from .action import Homomorphism
from .group import Coset, symmetric_group
from .permutation import Permutation, identity, validate
from .string_iso import LuksSolver
from .solver import BabaiSolver
from .standard_blocks import subset_action


@dataclass(frozen=True)
class Graph:
    n: int
    edges: frozenset[tuple[int, int]]

    def __init__(self, n: int, edges: Iterable[Iterable[int]] = ()):
        identity(n)
        normalized = set()
        for edge in edges:
            edge = tuple(edge)
            if len(edge) != 2 or any(type(v) is not int or not 0 <= v < n for v in edge) or edge[0] == edge[1]:
                raise ValueError("each edge must join two distinct vertices in range(n)")
            normalized.add(tuple(sorted(edge)))
        object.__setattr__(self, "n", n)
        object.__setattr__(self, "edges", frozenset(normalized))

    def relabel(self, permutation: Permutation) -> "Graph":
        permutation = validate(permutation, self.n)
        return Graph(self.n, ((permutation[a], permutation[b]) for a, b in self.edges))

    def validates_isomorphism(self, target: "Graph", permutation: Permutation) -> bool:
        return self.n == target.n and self.relabel(permutation) == target


@dataclass(frozen=True)
class GraphStringProblem:
    action: Homomorphism
    source: tuple[int, ...]
    target: tuple[int, ...]


def graph_string_problem(source: Graph, target: Graph) -> GraphStringProblem:
    if source.n != target.n:
        raise ValueError("graph orders differ")
    group = symmetric_group(source.n)
    pairs = tuple(combinations(range(source.n), 2))
    indices = {pair: i for i, pair in enumerate(pairs)}
    action = Homomorphism(group, len(pairs), (tuple(indices[tuple(sorted((g[a], g[b])))] for a, b in pairs) for g in group.generators))
    return GraphStringProblem(action, tuple(int(pair in source.edges) for pair in pairs), tuple(int(pair in target.edges) for pair in pairs))


def graph_isomorphisms(source: Graph, target: Graph, *, solver: LuksSolver | None = None) -> Coset | None:
    """Return the full vertex-isomorphism coset, or None if disproven.

    The default BabaiSolver starts with the known symmetric vertex action
    on unordered pairs. Passing a plain LuksSolver explicitly keeps that
    component's limited continuation behavior.
    """
    if source.n != target.n or len(source.edges) != len(target.edges):
        return None
    problem = graph_string_problem(source, target)
    solver = solver or BabaiSolver()
    if isinstance(solver, BabaiSolver) and source.n >= 3:
        # The pair action is faithful for n >= 3, so its inverse gives the
        # giant representation on the original vertex domain immediately.
        inverse_action = Homomorphism(problem.action.image, source.n, (problem.action.lift(g) for g in problem.action.image.generators))
        supports = tuple(frozenset(pair) for pair in combinations(range(source.n), 2))
        model = subset_action(inverse_action, supports)
        result = solver.solve_model(model, problem.source, problem.target)
    else:
        result = solver.solve(problem.action.image, problem.source, problem.target)
    lifted = problem.action.preimage_coset(result)
    if lifted is not None:
        if not source.validates_isomorphism(target, lifted.representative):
            raise AssertionError("invalid isomorphism returned by string solver")
        if not all(source.validates_isomorphism(source, g) for g in lifted.subgroup.generators):
            raise AssertionError("invalid automorphism returned by string solver")
    return lifted
