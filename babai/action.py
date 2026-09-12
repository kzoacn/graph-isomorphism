"""Constructive homomorphisms between permutation groups.

A generator and its image act together on a disjoint union of domains.
Stabilizer chains of this joint action give kernels and constructive lifts.
This avoids enumerating either group, even when the image is alternating.
"""

from collections.abc import Iterable, Sequence
from functools import cached_property

from .group import Coset, PermutationGroup, StabilizerChain, _points
from .permutation import Permutation, compose, identity, inverse, validate


def _projection_lift(target: Permutation, coordinates: tuple[int, ...], chain: StabilizerChain) -> Permutation | None:
    representative = identity(len(chain.base))
    for i, point in enumerate(coordinates):
        target_point = coordinates[target[i]]
        preimage = inverse(representative)[target_point]
        t = chain.transversals[i].get(preimage)
        if t is None:
            return None
        representative = compose(t, representative)
    return representative


class Homomorphism:
    """A homomorphism G -> Sym(m), given on ``G.generators`` in that order.

    Invalid assignments (those not respecting relations in G) are rejected.
    ``preimage(H)`` requires H <= image(phi); it does not hide an arbitrary
    subgroup intersection, which is itself an isomorphism-hard problem.
    """

    def __init__(self, group: PermutationGroup, image_degree: int, images: Iterable[Sequence[int]]):
        identity(image_degree)
        images = tuple(validate(p, image_degree) for p in images)
        if len(images) != len(group.generators):
            raise ValueError("provide exactly one image per group generator")
        self.domain = group
        self.image_degree = image_degree
        self.generator_images = images
        n = group.degree
        self._joint = PermutationGroup(n + image_degree, (g + tuple(n + j for j in p) for g, p in zip(group.generators, images)))
        if self._joint.order != group.order:
            raise ValueError("generator images do not define a homomorphism")

    @cached_property
    def image(self) -> PermutationGroup:
        return PermutationGroup(self.image_degree, self.generator_images)

    @cached_property
    def _image_chain(self) -> StabilizerChain:
        return self._joint.chain_with_prefix(range(self.domain.degree, self.domain.degree + self.image_degree))

    def __call__(self, element: Sequence[int]) -> Permutation:
        element = validate(element, self.domain.degree)
        lift = _projection_lift(element, tuple(range(self.domain.degree)), self._joint.chain)
        if lift is None:
            raise ValueError("element is not in the domain group")
        n = self.domain.degree
        return tuple(p - n for p in lift[n:])

    def lift(self, element: Sequence[int]) -> Permutation | None:
        """Return a preimage of an image permutation, or None outside the image."""
        element = validate(element, self.image_degree)
        n = self.domain.degree
        coordinates = tuple(range(n, n + self.image_degree))
        joint_lift = _projection_lift(element, coordinates, self._image_chain)
        return None if joint_lift is None else joint_lift[:n]

    @cached_property
    def kernel(self) -> PermutationGroup:
        n = self.domain.degree
        return PermutationGroup(n, (g[:n] for g in self._image_chain.generators_from(self.image_degree))).compact()

    def preimage(self, subgroup: PermutationGroup) -> PermutationGroup:
        if not subgroup.is_subgroup_of(self.image):
            raise ValueError("preimage requires a subgroup of the image")
        generators = list(self.kernel.generators)
        for g in subgroup.generators:
            lift = self.lift(g)
            if lift is None:
                raise AssertionError("image generator has no lift")
            generators.append(lift)
        return PermutationGroup(self.domain.degree, generators).compact()

    def preimage_coset(self, coset: Coset | None) -> Coset | None:
        if coset is None:
            return None
        representative = self.lift(coset.representative)
        if representative is None:
            return None
        return Coset(self.preimage(coset.subgroup), representative)

    def restrict(self, subgroup: PermutationGroup) -> "Homomorphism":
        if not subgroup.is_subgroup_of(self.domain):
            raise ValueError("restriction requires a subgroup of the domain")
        return Homomorphism(subgroup, self.image_degree, (self(g) for g in subgroup.generators))

    def then(self, other: "Homomorphism") -> "Homomorphism":
        if not self.image.is_subgroup_of(other.domain):
            raise ValueError("the first image is not contained in the second domain")
        return Homomorphism(self.domain, other.image_degree, (other(p) for p in self.generator_images))


def restriction_action(group: PermutationGroup, points: Iterable[int]) -> Homomorphism:
    """Action on an invariant ordered subset, relabeled as range(len(points))."""
    points = _points(points, group.degree)
    position = {p: i for i, p in enumerate(points)}
    if any(g[p] not in position for g in group.generators for p in points):
        raise ValueError("subset is not invariant under the group")
    images = (tuple(position[g[p]] for p in points) for g in group.generators)
    return Homomorphism(group, len(points), images)


def block_action(group: PermutationGroup, blocks: Iterable[Iterable[int]]) -> Homomorphism:
    """Action on an invariant partition of the entire domain."""
    blocks = tuple(tuple(block) for block in blocks)
    flattened = tuple(p for block in blocks for p in block)
    _points(flattened, group.degree)
    if len(flattened) != group.degree or any(not block for block in blocks):
        raise ValueError("blocks must be nonempty and partition the entire domain")
    location = {p: i for i, block in enumerate(blocks) for p in block}
    images = []
    for g in group.generators:
        image = []
        for block in blocks:
            targets = {location[g[p]] for p in block}
            if len(targets) != 1:
                raise ValueError("partition is not invariant under the group")
            image.append(targets.pop())
        images.append(tuple(image))
    return Homomorphism(group, len(blocks), images)

