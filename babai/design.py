"""The Design Lemma: produce a witness, exposing every individualized point.

Helfgott, Proposition 5.1. A caller searching for isomorphisms must pay for
the images of the witness tuple; a lexicographically first witness is not
an invariant of an unlabeled structure.
"""

from dataclasses import dataclass
from fractions import Fraction
from itertools import permutations

from .configuration import Configuration, is_clique, is_coherent


@dataclass(frozen=True)
class DesignWitness:
    fixed: tuple[int, ...]
    vertex_colors: tuple[int, ...]
    dominant: tuple[int, ...] | None
    configuration: Configuration | None


def design_at(configuration: Configuration, fixed: tuple[int, ...], alpha: Fraction = Fraction(1, 2)) -> DesignWitness | None:
    """Evaluate the Design Lemma alternatives at a prescribed ordered tuple."""
    if not Fraction(1, 2) <= alpha < 1:
        raise ValueError("require 1/2 <= alpha < 1")
    conditioned = configuration.condition(fixed)
    vertex_colors = conditioned.skeleton(1).colors
    dominant = next((vertices for _, vertices in conditioned.vertex_classes() if len(vertices) > alpha * configuration.degree), None)
    if dominant is None:
        return DesignWitness(fixed, vertex_colors, None, None)
    if conditioned.arity >= 2:
        binary = conditioned.skeleton(2).induced(dominant)
        if not is_clique(binary):
            return DesignWitness(fixed, vertex_colors, dominant, binary)
    return None


def design_lemma(configuration: Configuration, alpha: Fraction = Fraction(1, 2)) -> DesignWitness:
    n, k = configuration.degree, configuration.arity
    if not 2 <= k <= n // 2:
        raise ValueError("Design Lemma requires 2 <= arity <= degree/2")
    if not Fraction(1, 2) <= alpha < 1:
        raise ValueError("require 1/2 <= alpha < 1")
    if not is_coherent(configuration):
        raise ValueError("Design Lemma requires a coherent configuration")
    if any(len(twins) > alpha * n for twins in configuration.twin_classes()):
        raise ValueError("a twin class exceeds alpha * degree")
    for length in range(k):
        for fixed in permutations(range(n), length):
            result = design_at(configuration, fixed, alpha)
            if result is not None:
                return result
    raise AssertionError("Design Lemma hypotheses hold but no witness exists")

