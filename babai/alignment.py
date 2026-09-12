"""Constructive alignment of Split-or-Johnson outputs in S_m or A_m.

These are the polynomial operations used before pulling a discovered
structure back through a giant representation (Babai section 11).
"""

from .group import Coset, PermutationGroup, symmetric_group
from .partition import even_subgroup
from .permutation import compose, from_cycles, identity, inverse, sign
from .split_johnson import SplitResult


def split_stabilizer(result: SplitResult) -> PermutationGroup:
    """Automorphisms of the discovered structure, in domain-index coordinates."""
    degree = len(result.domain)
    positions = {v: i for i, v in enumerate(result.domain)}
    johnson_vertices = set() if result.johnson is None else set(result.johnson.vertices)
    generators = []
    for blocks in result.parts:
        vertices = {v for block in blocks for v in block}
        if vertices == johnson_vertices:
            embedded = result.johnson
            for atom_permutation in symmetric_group(embedded.model.atoms).generators:
                local = embedded.model.induce(atom_permutation)
                p = list(range(degree))
                for i, v in enumerate(embedded.vertices):
                    p[positions[v]] = positions[embedded.vertices[local[i]]]
                generators.append(tuple(p))
        else:
            local_blocks = tuple(tuple(positions[v] for v in block) for block in blocks)
            for block in local_blocks:
                generators.extend(from_cycles(degree, (block[0], v)) for v in block[1:])
            for block in local_blocks[1:]:
                generators.append(from_cycles(degree, *zip(local_blocks[0], block)))
    return PermutationGroup(degree, generators)


def align_splits(source: SplitResult, target: SplitResult, *, even: bool = False) -> Coset | None:
    """Align structures and every corresponding individualized point.

    Domain vertices can have different labels; the returned permutation
    maps *indices* in source.domain to indices in target.domain. Point
    choices must lie in these domains (the full, non-bipartite setting).
    """
    if len(source.domain) != len(target.domain):
        return None
    if (source.johnson is None) != (target.johnson is None):
        return None
    ordinary = source.colored_partition().align(target.colored_partition())
    if ordinary is None:
        return None
    source_positions = {v: i for i, v in enumerate(source.domain)}
    target_positions = {v: i for i, v in enumerate(target.domain)}
    representative = list(ordinary.representative)
    if source.johnson is not None:
        left, right = source.johnson, target.johnson
        if (left.model.atoms, left.model.subset_size, left.model.distance_colors) != (right.model.atoms, right.model.subset_size, right.model.distance_colors):
            return None
        left_color = next(i for i, blocks in enumerate(source.parts) if set(left.vertices) == {v for block in blocks for v in block})
        right_color = next(i for i, blocks in enumerate(target.parts) if set(right.vertices) == {v for block in blocks for v in block})
        if left_color != right_color:
            return None
        index = {s: i for i, s in enumerate(right.model.subsets)}
        for v, subset in zip(left.vertices, left.model.subsets):
            representative[source_positions[v]] = target_positions[right.vertices[index[subset]]]
    representative = tuple(representative)
    group = split_stabilizer(source)
    if even:
        if sign(representative) == -1:
            odd = next((g for g in group.generators if sign(g) == -1), None)
            if odd is None:
                return None
            representative = compose(odd, representative)
        group = even_subgroup(group)

    if len(source.choices) != len(target.choices):
        return None
    constraints = {}
    for left, right in zip(source.choices, target.choices):
        if left.kind != right.kind or len(left.vertices) != len(right.vertices):
            return None
        for a, b in zip(left.vertices, right.vertices):
            if a not in source_positions or b not in target_positions:
                raise ValueError("alignment requires choices within the aligned domains")
            if a in constraints and constraints[a] != b:
                return None
            constraints[a] = b
    if len(set(constraints.values())) != len(constraints):
        return None
    inverse_representative = inverse(representative)
    alignment = group.transporter(tuple(source_positions[a] for a in constraints), tuple(inverse_representative[target_positions[b]] for b in constraints.values()))
    return None if alignment is None else alignment.multiply_right(representative)

