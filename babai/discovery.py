"""Assemble canonical structure discovery and replay its choices in a pair."""

from fractions import Fraction
from itertools import product

from .configuration import Configuration
from .design import design_at, design_lemma
from .split_johnson import Choice, EmbeddedJohnson, IncompatibleChoices, SplitResult, split_or_johnson
from .standard_blocks import StandardBlocks


def fiber_configurations(model: StandardBlocks, source: tuple, target: tuple) -> tuple[Configuration, Configuration] | None:
    source_colors, target_colors = model.fiber_colors(source), model.fiber_colors(target)
    if source_colors is None or target_colors is None:
        return None
    m, k = model.action.image_degree, model.subset_size
    index = {subset: i for i, subset in enumerate(model.subsets)}
    palette = {}
    for c in source_colors + target_colors:
        if c not in palette:
            palette[c] = len(palette) + 1
    result = []
    for block_colors in (source_colors, target_colors):
        colors = []
        for points in product(range(m), repeat=k):
            subset = frozenset(points)
            colors.append(0 if len(subset) != k else palette[block_colors[index[subset]]])
        result.append(Configuration(m, k, tuple(colors)))
    return tuple(result)


def discover_relation(configuration: Configuration, alpha: Fraction = Fraction(2, 3), *, fixed: tuple[int, ...] | None = None, split_plan: tuple[Choice, ...] | None = None) -> SplitResult:
    """Design Lemma followed, when needed, by Split-or-Johnson.

    Choice labels refer to the original relational domain. A target call
    receives images of the source prefix and of its subsequent SoJ plan.
    """
    if not Fraction(2, 3) <= alpha < 1:
        raise ValueError("require 2/3 <= alpha < 1")
    witness = design_lemma(configuration, alpha) if fixed is None else design_at(configuration, fixed, alpha)
    if witness is None:
        raise IncompatibleChoices("the proposed Design Lemma tuple is not a witness")
    conditioned = configuration.condition(witness.fixed)
    choices = (Choice("design", witness.fixed),)
    m = configuration.degree
    if witness.dominant is None:
        if split_plan:
            raise IncompatibleChoices("a coloring-only branch cannot consume a SoJ plan")
        parts = tuple((vertices,) for _, vertices in conditioned.vertex_classes())
        return SplitResult(tuple(range(m)), parts, None, choices, alpha * m)
    dominant = witness.dominant
    positions = {v: i for i, v in enumerate(dominant)}
    local_plan = None
    if split_plan is not None:
        try:
            local_plan = tuple(Choice(choice.kind, tuple(positions[v] for v in choice.vertices)) for choice in split_plan)
        except KeyError as exc:
            raise IncompatibleChoices("SoJ choice is outside the target dominant class") from exc
    beta = alpha * m / len(dominant)
    local = split_or_johnson(witness.configuration, beta, plan=local_plan)
    parts = tuple((vertices,) for _, vertices in conditioned.vertex_classes() if vertices != dominant)
    parts += tuple(tuple(tuple(dominant[v] for v in block) for block in blocks) for blocks in local.parts)
    choices += tuple(Choice(choice.kind, tuple(dominant[v] for v in choice.vertices)) for choice in local.choices)
    embedded = None if local.johnson is None else EmbeddedJohnson(tuple(dominant[v] for v in local.johnson.vertices), local.johnson.model)
    result = SplitResult(tuple(range(m)), parts, embedded, choices, alpha * m)
    result.colored_partition()
    return result


def embed_split(local: SplitResult, domain: tuple[int, ...], degree: int, *, prefix: tuple[Choice, ...] = (), alpha: Fraction = Fraction(2, 3)) -> SplitResult:
    """Embed a result on an invariant subset and color its complement."""
    if len(domain) != len(local.domain) or len(set(domain)) != len(domain) or any(not 0 <= v < degree for v in domain):
        raise ValueError("invalid embedding domain")
    positions = {v: domain[i] for i, v in enumerate(local.domain)}
    complement = tuple(v for v in range(degree) if v not in domain)
    parts = ((complement,),) if complement else ()
    parts += tuple(tuple(tuple(positions[v] for v in block) for block in blocks) for blocks in local.parts)
    embedded = None
    if local.johnson is not None:
        embedded = EmbeddedJohnson(tuple(positions[v] for v in local.johnson.vertices), local.johnson.model)
    choices = prefix + tuple(Choice(choice.kind, tuple(positions[v] for v in choice.vertices)) for choice in local.choices)
    result = SplitResult(tuple(range(degree)), parts, embedded, choices, alpha * degree)
    result.colored_partition()
    return result


def coloring_split(degree: int, cells: tuple[tuple[int, ...], ...], *, choices: tuple[Choice, ...] = (), alpha: Fraction = Fraction(2, 3)) -> SplitResult:
    result = SplitResult(tuple(range(degree)), tuple((cell,) for cell in cells if cell), None, choices, alpha * degree)
    result.colored_partition()
    return result

