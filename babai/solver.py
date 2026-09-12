"""Babai recursion through complete subset actions and local certificates.

The modelled recursion also handles non-giant image groups, using the
orbital-configuration route of Babai section 13. Every costly structural
step ends with a smaller actual window or a smaller ideal domain.
"""

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from math import factorial

from .action import Homomorphism, block_action, restriction_action
from .aggregation import certificate_table, nonfull_relations
from .alignment import align_splits, split_stabilizer
from .blocks import maximal_block_system, nontrivial_block_system
from .configuration import individualization_compatible, refine_joint
from .discovery import coloring_split, discover_relation, embed_split, fiber_configurations
from .group import Coset, PermutationGroup, merge_isomorphism_cosets, symmetric_group
from .local_certificates import LocalCertificates, contains_alternating
from .orbitals import break_double_transitivity, is_doubly_transitive, orbital_configuration, orbital_graphs
from .permutation import compose, identity, pullback
from .split_johnson import Choice, IncompatibleChoices, SplitResult, split_or_johnson
from .standard_blocks import StandardBlocks, descend_action, subset_action
from .string_iso import GiantActionContext, LuksSolver, small_quotient_bound
from .top_action import giant_top_action
from .tower import BlockTower

ALPHA = Fraction(2, 3)


@dataclass(frozen=True)
class RecursionEvent:
    kind: str
    actual_degree: int
    ideal_degree: int
    child_actual_degree: int | None = None
    child_ideal_degree: int | None = None
    branches: int = 1


def locality(actual_degree: int) -> int:
    return max(9, actual_degree.bit_length() + 2)


class BabaiSolver(LuksSolver):
    """String isomorphism with the Babai continuation of Luks reduction.

    ``events`` exposes the actual/ideal degrees and enumeration indices
    used by the implementation; it is useful for auditing the recurrence.
    """

    def __init__(self, *, prefer_structural: bool = False):
        super().__init__(self._large_quotient)
        self.events: list[RecursionEvent] = []
        self.prefer_structural = prefer_structural

    def _large_quotient(self, context: GiantActionContext) -> Coset | None:
        supports = [frozenset()] * context.group.degree
        for i, block in enumerate(context.blocks):
            for v in block:
                supports[v] = frozenset((i,))
        model = subset_action(context.quotient, tuple(supports))
        return self.solve_model(model, context.source, context.target, tower=context.tower)

    def solve_model(self, model: StandardBlocks, source: tuple, target: tuple, *, tower: BlockTower | None = None) -> Coset | None:
        source, target = tuple(source), tuple(target)
        group = model.action.domain
        if len(source) != group.degree or len(target) != group.degree:
            raise ValueError("strings have the wrong degree")
        if tower is None:
            tower = self.current_tower
            if tower is None or tower.degree != group.degree:
                tower = BlockTower(group.degree)
        tower.validate_for(group)
        self._towers.append(tower)
        try:
            return self._modelled(model, source, target)
        finally:
            self._towers.pop()

    def _modelled(self, model: StandardBlocks, source: tuple, target: tuple) -> Coset | None:
        self.calls += 1
        action, group = model.action, model.action.domain
        n, m = group.degree, action.image_degree
        if Counter(source) != Counter(target):
            return None
        if source == target and all(all(source[g[v]] == source[v] for v in range(n)) for g in group.generators):
            return Coset(group, identity(n))
        if not group.generators:
            return None

        source_colors, target_colors = model.fiber_colors(source), model.fiber_colors(target)
        if (source_colors is None) != (target_colors is None):
            return None
        if source_colors is not None and len(model.blocks) < n:
            # A polynomial reduction of actual degree by at least two.
            projection = block_action(group, model.blocks)
            images = []
            for g in projection.image.generators:
                lift = projection.lift(g)
                if lift is None:
                    raise AssertionError("block generator has no lift")
                images.append(action(lift))
            descended = Homomorphism(projection.image, m, images)
            child = subset_action(descended, model.subsets)
            self.events.append(RecursionEvent("constant_fibers", n, m, len(model.blocks), m))
            result = self.solve_model(child, source_colors, target_colors, tower=BlockTower(len(model.blocks)))
            return projection.preimage_coset(result)

        giant = contains_alternating(action.image)
        if giant and source_colors is not None:
            if model.subset_size == 1:
                result = self.solve(action.image, source_colors, target_colors, tower=BlockTower(m))
                return action.preimage_coset(result)
            raw_source, raw_target = fiber_configurations(model, source, target)
            source_twins = next((c for c in raw_source.twin_classes() if len(c) > ALPHA * m), None)
            target_twins = next((c for c in raw_target.twin_classes() if len(c) > ALPHA * m), None)
            if (source_twins is None) != (target_twins is None):
                return None
            if source_twins is not None:
                if len(source_twins) != len(target_twins):
                    return None
                return self._large_symmetry(model, source, target, source_twins, target_twins)
            refined_source, refined_target = refine_joint(raw_source, raw_target).configurations
            return self._relation_family(model, source, target, refined_source, refined_target, tuple(range(m)), tuple(range(m)))

        if giant:
            top = giant_top_action(action, source, target, self, max_orbit_size=len(model.blocks[0]))
            self.events.append(RecursionEvent("top_action", n, m, len(model.blocks[0])))
            if top.has_alternating_automorphisms:
                return top.isomorphisms

        if m <= 10 * locality(n) and (giant or not self.prefer_structural):
            return self._enumerate_action(action, source, target, factorial(m), n // m, "small_ideal")
        if action.image.order <= small_quotient_bound(m) and not self.prefer_structural:
            return self._enumerate_action(action, source, target, small_quotient_bound(m), n // m, "small_image")
        if giant:
            return self._aggregate(model, source, target)

        # The actual domain remains a complete subset family even if its
        # group is intransitive. Decompose the *ideal* action first, so the
        # cost of breaking double transitivity is followed by real progress.
        ideal_orbits = action.image.orbits()
        if len(ideal_orbits) > 1:
            structure = coloring_split(m, ideal_orbits)
            return self._solve_structure(model, Coset(group, identity(n)), source, target, structure, allow_orbit_descent=True)
        blocks = nontrivial_block_system(action.image)
        if blocks is not None:
            blocks = maximal_block_system(action.image, blocks)
            structure = SplitResult(tuple(range(m)), (blocks,), None, (), ALPHA * m)
            return self._solve_structure(model, Coset(group, identity(n)), source, target, structure)
        if is_doubly_transitive(action.image):
            fixed = break_double_transitivity(action.image)
            remaining = tuple(v for v in range(m) if v not in fixed)
            structure = coloring_split(m, (remaining,) + tuple((v,) for v in fixed), choices=(Choice("transitivity", fixed),))
            return self._individualize_group_structure(model, source, target, structure, allow_orbit_descent=True)
        structure = split_or_johnson(orbital_configuration(action.image), ALPHA)
        return self._individualize_group_structure(model, source, target, structure)

    def _enumerate_action(self, action: Homomorphism, source: tuple, target: tuple, bound: int, orbit_bound: int, kind: str) -> Coset | None:
        kernel = action.kernel
        orbits = kernel.orbits()
        if any(len(orbit) > orbit_bound for orbit in orbits):
            raise AssertionError("kernel recursion exceeds its proved orbit bound")
        self.events.append(RecursionEvent(kind, action.domain.degree, action.image_degree, max(map(len, orbits), default=0), branches=action.image.order))

        def branches():
            for image in action.image.elements(max_order=bound):
                lift = action.lift(image)
                if lift is None:
                    raise AssertionError("image element has no lift")
                result = Coset(kernel, lift)
                for orbit in orbits:
                    result = self.solve_window(result, source, target, orbit)
                    if result is None:
                        break
                yield result

        return merge_isomorphism_cosets(branches())

    def _individualize_group_structure(self, model, source, target, structure, *, allow_orbit_descent=False):
        action = model.action
        image_stabilizer = action.image.pointwise_stabilizer(structure.fixed_points)
        subgroup = action.preimage(image_stabilizer)
        if not image_stabilizer.is_subgroup_of(split_stabilizer(structure)):
            raise AssertionError("group-derived structure is not invariant after its declared choices")
        bound = action.image_degree ** len(structure.fixed_points)
        representatives = action.image.coset_representatives(image_stabilizer, max_index=bound)
        self.events.append(RecursionEvent("group_individualization", action.domain.degree, action.image_degree, branches=len(representatives)))
        child_model = model.restrict_group(subgroup)

        def branches():
            for image in representatives:
                lift = action.lift(image)
                if lift is None:
                    raise AssertionError("individualization representative has no lift")
                result = self._solve_structure(child_model, Coset(subgroup, identity(subgroup.degree)), source, pullback(target, lift), structure, allow_orbit_descent=allow_orbit_descent)
                yield None if result is None else result.multiply_right(lift)

        return merge_isomorphism_cosets(branches())

    def _aligned_structure(self, model, source, target, left, right):
        even = model.action.image.order != factorial(model.action.image_degree)
        alignment = align_splits(left, right, even=even)
        alignment = model.action.preimage_coset(alignment)
        if alignment is None:
            return None
        return self._solve_structure(model, alignment, source, target, left)

    def _solve_structure(self, model, candidates, source, target, structure, *, allow_orbit_descent=False):
        n, m = model.action.domain.degree, model.action.image_degree
        if structure.domain != tuple(range(m)):
            raise ValueError("structure is not on the model's ideal domain")
        subgroup_action = model.action.restrict(candidates.subgroup)
        if not subgroup_action.image.is_subgroup_of(split_stabilizer(structure)):
            raise AssertionError("candidate group does not preserve the discovered structure")
        colors = [0] * m
        classes = []
        for c, blocks in enumerate(structure.parts):
            vertices = tuple(sorted(v for block in blocks for v in block))
            classes.append(vertices)
            for v in vertices:
                colors[v] = c
        windows = sorted(model.color_windows(tuple(colors)), key=lambda item: (-len(item[1]), item[0]))
        result = candidates
        for signature, window in windows:
            if len(window) <= ALPHA * n:
                self.events.append(RecursionEvent("color_window", n, m, len(window)))
                result = self.solve_window(result, source, target, window)
            else:
                if len(signature) != 1 or signature[0][1] != model.subset_size:
                    raise AssertionError("a dominant induced window must use one ideal color")
                color = signature[0][0]
                vertices = classes[color]
                if len(vertices) <= ALPHA * m:
                    raise AssertionError("dominant window has no dominant ideal color")
                result = self._dominant_structure_window(model, result, source, target, window, vertices, structure.parts[color], structure, allow_orbit_descent)
            if result is None:
                return None
        return result

    def _dominant_structure_window(self, model, candidates, source, target, window, vertices, blocks, structure, allow_orbit_descent):
        action = model.action.restrict(candidates.subgroup)
        n, m = action.domain.degree, action.image_degree
        johnson = structure.johnson if structure.johnson is not None and set(structure.johnson.vertices) == set(vertices) else None
        images = []
        if johnson is not None:
            local_positions = {v: i for i, v in enumerate(johnson.vertices)}
            for g in candidates.subgroup.generators:
                image = action(g)
                local = tuple(local_positions[image[v]] for v in johnson.vertices)
                atoms = johnson.model.recover(local)
                if atoms is None:
                    raise AssertionError("candidate is not a Johnson automorphism")
                images.append(atoms)
            new_degree = johnson.model.atoms
            labels = dict(zip(johnson.vertices, johnson.model.subsets))
            kind = "johnson_descent"
        elif len(blocks) >= 2:
            block_index = {v: i for i, block in enumerate(blocks) for v in block}
            for g in candidates.subgroup.generators:
                image = action(g)
                images.append(tuple(block_index[image[block[0]]] for block in blocks))
            new_degree = len(blocks)
            labels = {v: frozenset((block_index[v],)) for v in vertices}
            kind = "block_descent"
        else:
            if not allow_orbit_descent or len(vertices) >= m:
                raise AssertionError("the structural recursion makes no progress")
            positions = {v: i for i, v in enumerate(vertices)}
            for g in candidates.subgroup.generators:
                image = action(g)
                images.append(tuple(positions[image[v]] for v in vertices))
            new_degree = len(vertices)
            labels = {v: frozenset((positions[v],)) for v in vertices}
            kind = "ideal_orbit"
        if kind != "ideal_orbit" and 2 * new_degree > m:
            raise AssertionError("structural descent did not halve the ideal domain")
        new_action = Homomorphism(candidates.subgroup, new_degree, images)
        self.events.append(RecursionEvent(kind, n, m, len(window), new_degree))
        supports = [frozenset()] * n
        for v in window:
            supports[v] = frozenset().union(*(labels[p] for p in model.supports[v]))
        trivial_support = any(not 0 < len(supports[v]) < new_degree for v in window)
        if new_degree <= 10 * locality(len(window)) and (not self.prefer_structural or trivial_support):
            return self._auxiliary_window(new_action, candidates, source, target, window, None)
        by_size = {}
        for v in window:
            size = len(supports[v])
            if not 0 < size < new_degree:
                raise AssertionError("large reduced ideal action has a trivial support")
            by_size.setdefault(size, []).append(v)
        result = candidates
        for _, vertices in sorted(by_size.items(), key=lambda item: (-len(item[1]), item[0])):
            result = self._auxiliary_window(new_action, result, source, target, tuple(vertices), tuple(supports))
            if result is None:
                return None
        return result

    def _auxiliary_window(self, action, candidates, source, target, window, supports):
        action = action.restrict(candidates.subgroup)
        projection, descended = descend_action(action, window)
        x = tuple(source[v] for v in window)
        y = tuple(target[candidates.representative[v]] for v in window)
        tower = self.current_tower.restrict(window)
        if supports is None:
            # Small reduced top group. Its kernel fixes either all former
            # Johnson vertices or every block of an equipartition.
            # In either case every kernel orbit is at most 2/3 of the window.
            self._towers.append(tower)
            try:
                result = self._enumerate_action(descended, x, y, factorial(descended.image_degree), 2 * len(window) // 3, "small_structural_image")
            finally:
                self._towers.pop()
        else:
            child_model = subset_action(descended, tuple(supports[v] for v in window))
            result = self.solve_model(child_model, x, y, tower=tower)
        result = projection.preimage_coset(result)
        return None if result is None else result.multiply_right(candidates.representative)

    def _relation_family(self, model, source, target, left_configuration, right_configuration, left_domain, right_domain):
        m, d = model.action.image_degree, len(left_domain)
        if d != len(right_domain):
            return None
        alpha = ALPHA * m / d
        left_local = discover_relation(left_configuration, alpha)
        left = embed_split(left_local, left_domain, m)
        ambient = symmetric_group(d)
        fixed = left_local.fixed_points
        stabilizer = ambient.pointwise_stabilizer(fixed)
        representatives = ambient.coset_representatives(stabilizer, max_index=d ** len(fixed))
        self.events.append(RecursionEvent("relation_individualization", model.action.domain.degree, m, branches=len(representatives)))

        def branches():
            for p in representatives:
                if not individualization_compatible(left_configuration, right_configuration, fixed, tuple(p[v] for v in fixed)):
                    continue
                try:
                    right_local = discover_relation(right_configuration, alpha, fixed=tuple(p[v] for v in left_local.choices[0].vertices), split_plan=tuple(Choice(c.kind, tuple(p[v] for v in c.vertices)) for c in left_local.choices[1:]))
                    right = embed_split(right_local, right_domain, m)
                except IncompatibleChoices:
                    continue
                yield self._aligned_structure(model, source, target, left, right)

        return merge_isomorphism_cosets(branches())

    def _large_symmetry(self, model, source, target, left_domain, right_domain):
        m, n = model.action.image_degree, model.action.domain.degree
        if len(left_domain) != len(right_domain):
            return None
        left = coloring_split(m, (tuple(left_domain), tuple(v for v in range(m) if v not in left_domain)))
        right = coloring_split(m, (tuple(right_domain), tuple(v for v in range(m) if v not in right_domain)))
        even = model.action.image.order != factorial(m)
        candidates = model.action.preimage_coset(align_splits(left, right, even=even))
        if candidates is None:
            return None
        action = model.action.restrict(candidates.subgroup)
        action = action.then(restriction_action(action.image, tuple(left_domain)))
        domain = set(left_domain)
        window = tuple(v for v, support in enumerate(model.supports) if support.issubset(domain))
        projection, descended = descend_action(action, window)
        tower = self.current_tower.restrict(window)
        self._towers.append(tower)
        try:
            top = giant_top_action(descended, tuple(source[v] for v in window), tuple(target[candidates.representative[v]] for v in window), self, max_orbit_size=len(window) // len(left_domain))
        finally:
            self._towers.pop()
        if not top.has_alternating_automorphisms:
            raise AssertionError("the claimed large top symmetry did not lift to automorphisms")
        result = projection.preimage_coset(top.isomorphisms)
        if result is None:
            return None
        result = result.multiply_right(candidates.representative)
        self.events.append(RecursionEvent("large_symmetry", n, m, len(window) // len(left_domain), len(left_domain)))
        colors = tuple(0 if v in domain else 1 for v in range(m))
        for _, other in model.color_windows(colors):
            if other == window:
                continue
            if len(other) > ALPHA * n:
                raise AssertionError("large-symmetry remainder exceeds 2/3 of the actual domain")
            result = self.solve_window(result, source, target, other)
            if result is None:
                return None
        return result

    def _aggregate(self, model, source, target):
        action, m = model.action, model.action.image_degree
        k = locality(action.domain.degree)
        if not 10 * k < m:
            raise AssertionError("local-certificate aggregation requires k < m/10")
        engine = LocalCertificates(action, self)
        left = certificate_table(engine, source, k)
        right = certificate_table(engine, target, k)
        return self._aggregate_tables(model, left, right)

    def _aggregate_tables(self, model, left, right):
        source, target = left.source, right.source
        action, m = model.action, model.action.image_degree
        a, b = left.structure, right.structure
        if a.kind != b.kind or len(a.support) != len(b.support):
            return None
        if a.kind == "orbit_partition":
            make = lambda p: SplitResult(tuple(range(m)), tuple(blocks for _, blocks in p.classes), None, (), ALPHA * m)
            return self._aligned_structure(model, source, target, make(a.partition), make(b.partition))
        if a.large_orbit is not None:
            if len(a.large_orbit) != len(b.large_orbit):
                return None
            if len(a.large_orbit) <= ALPHA * m:
                make = lambda orbit: coloring_split(m, (orbit, tuple(v for v in range(m) if v not in orbit)))
                return self._aligned_structure(model, source, target, make(a.large_orbit), make(b.large_orbit))
            if a.kind == "large_symmetry":
                return self._large_symmetry(model, source, target, a.large_orbit, b.large_orbit)
            return self._nongiant_certificate_orbit(model, source, target, a, b)
        left_domain = tuple(v for v in range(m) if v not in a.support)
        right_domain = tuple(v for v in range(m) if v not in b.support)
        if len(left_domain) <= ALPHA * m:
            return self._aligned_structure(model, source, target, coloring_split(m, (left_domain, a.support)), coloring_split(m, (right_domain, b.support)))
        x, y = nonfull_relations(left, right, left_domain, right_domain)
        x, y = refine_joint(x, y).configurations
        return self._relation_family(model, source, target, x, y, left_domain, right_domain)

    def _nongiant_certificate_orbit(self, model, source, target, left, right):
        m, n = model.action.image_degree, model.action.domain.degree
        left_domain, right_domain = left.large_orbit, right.large_orbit
        left_group = restriction_action(left.group, left_domain).image
        right_group = restriction_action(right.group, right_domain).image
        fixed = break_double_transitivity(left_group)
        remaining_indices = tuple(v for v in range(len(left_domain)) if v not in fixed)
        remaining = tuple(left_domain[i] for i in remaining_indices)
        prefix = (Choice("transitivity", tuple(left_domain[i] for i in fixed)),)
        if len(remaining) <= ALPHA * m:
            left_structure = coloring_split(m, (tuple(v for v in range(m) if v not in left_domain), remaining) + tuple((left_domain[i],) for i in fixed), choices=prefix)
            local_structure = None
        else:
            residual = restriction_action(left_group.pointwise_stabilizer(fixed), remaining_indices).image
            raw = orbital_graphs(residual)[0]
            coherent = refine_joint(raw).configurations[0]
            local_structure = split_or_johnson(coherent, ALPHA * m / len(remaining))
            left_structure = embed_split(local_structure, remaining, m, prefix=prefix)
        positions = {v: i for i, v in enumerate(left_domain)}
        all_fixed = tuple(positions[v] for v in left_structure.fixed_points)
        ambient = symmetric_group(len(left_domain))
        representatives = ambient.coset_representatives(ambient.pointwise_stabilizer(all_fixed), max_index=len(left_domain) ** len(all_fixed))
        self.events.append(RecursionEvent("certificate_orbital_choices", n, m, branches=len(representatives) * max(1, len(right_domain) - 1)))

        def branches():
            for p in representatives:
                target_fixed = tuple(p[i] for i in fixed)
                target_remaining_indices = tuple(v for v in range(len(right_domain)) if v not in target_fixed)
                target_remaining = tuple(right_domain[i] for i in target_remaining_indices)
                target_prefix = (Choice("transitivity", tuple(right_domain[i] for i in target_fixed)),)
                if local_structure is None:
                    right_structure = coloring_split(m, (tuple(v for v in range(m) if v not in right_domain), target_remaining) + tuple((right_domain[i],) for i in target_fixed), choices=target_prefix)
                    yield self._aligned_structure(model, source, target, left_structure, right_structure)
                    continue
                residual = restriction_action(right_group.pointwise_stabilizer(target_fixed), target_remaining_indices).image
                if len(residual.orbits()) != 1 or is_doubly_transitive(residual):
                    continue
                target_positions = {v: i for i, v in enumerate(target_remaining)}
                plan = tuple(Choice(c.kind, tuple(target_positions[right_domain[p[remaining_indices[v]]]] for v in c.vertices)) for c in local_structure.choices)
                for raw in orbital_graphs(residual):
                    coherent = refine_joint(raw).configurations[0]
                    try:
                        target_local = split_or_johnson(coherent, ALPHA * m / len(remaining), plan=plan)
                    except IncompatibleChoices:
                        continue
                    right_structure = embed_split(target_local, target_remaining, m, prefix=target_prefix)
                    yield self._aligned_structure(model, source, target, left_structure, right_structure)

        return merge_isomorphism_cosets(branches())
