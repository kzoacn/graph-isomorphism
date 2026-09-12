"""Local-certificate tables, fullness structure, and local-guide relations.

These provide the data used by the giant-action recursion. In particular,
the relation colors are equivalence classes in a shared two-object local
guide, not independent hashes of certificates from the two strings.
"""

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from itertools import combinations, permutations, product

from .action import restriction_action
from .configuration import Configuration
from .group import Coset, PermutationGroup
from .local_certificates import LocalCertificate, LocalCertificates, contains_alternating
from .partition import ColoredPartition
from .permutation import inverse


@dataclass(frozen=True)
class FullnessStructure:
    group: PermutationGroup
    support: tuple[int, ...]
    kind: str
    partition: ColoredPartition | None
    large_orbit: tuple[int, ...] | None


def fullness_structure(image: PermutationGroup) -> FullnessStructure:
    """Helfgott section 6.2's cases; no isomorphism answer is implied."""
    orbits = image.orbits()
    support = tuple(sorted(v for orbit in orbits if len(orbit) > 1 for v in orbit))
    if 2 * len(support) < image.degree:
        return FullnessStructure(image, support, "small_support", None, None)
    large_orbit = next((orbit for orbit in orbits if 2 * len(orbit) > image.degree), None)
    if large_orbit is None:
        classes = {}
        for orbit in orbits:
            classes.setdefault(len(orbit), []).append(orbit)
        partition = ColoredPartition(image.degree, tuple((size, tuple(parts)) for size, parts in sorted(classes.items())))
        return FullnessStructure(image, support, "orbit_partition", partition, None)
    giant = contains_alternating(restriction_action(image, large_orbit).image)
    return FullnessStructure(image, support, "large_symmetry" if giant else "large_nongiant_orbit", None, large_orbit)


@dataclass(frozen=True)
class CertificateTable:
    engine: LocalCertificates
    source: tuple
    arity: int
    certificates: dict[tuple[int, ...], LocalCertificate]
    full_group: PermutationGroup
    structure: FullnessStructure


def certificate_table(engine: LocalCertificates, source: Sequence, arity: int) -> CertificateTable:
    """Compute every k-subset certificate and generate F from the full ones.

    This costs m^O(k), with recursive work as documented by LocalCertificates.
    The main algorithm must choose k=O(log n) and use the small-m reduction
    unless k < m/10. This utility itself also accepts larger k for direct
    certificate computations; it does not claim the main recurrence then.
    """
    source = tuple(source)
    if type(arity) is not int or not 1 <= arity <= engine.action.image_degree:
        raise ValueError("invalid certificate arity")
    certificates = {}
    generators = []
    degree = engine.action.domain.degree
    for test_set in combinations(range(engine.action.image_degree), arity):
        certificate = engine.build(source, test_set)
        certificates[test_set] = certificate
        if certificate.full:
            generators.extend(certificate.certificate_group.generators)
            if len(generators) > max(1, degree**2):
                generators = list(PermutationGroup(degree, generators).compact().generators)
    group = PermutationGroup(degree, generators).compact()
    image = engine.action.restrict(group).image
    return CertificateTable(engine, source, arity, certificates, group, fullness_structure(image))


def local_guide_relations(domains: tuple[tuple[int, ...], tuple[int, ...]], arity: int, equivalent: Callable[[int, tuple[int, ...], int, tuple[int, ...]], bool]) -> tuple[Configuration, Configuration]:
    """Color ordered tuples by a local guide's equivalence classes.

    ``equivalent(i,u,j,v)`` must be an equivalence relation saying whether
    an allowed local isomorphism maps the ordered tuple u to v. Tuple
    entries are original domain labels. Output domains use local indices.
    Repeated tuples have color 0; all other colors have joint meanings.
    """
    if len(domains) != 2 or len(domains[0]) != len(domains[1]) or not 1 <= arity <= len(domains[0]):
        raise ValueError("require two equal-size domains and a valid arity")
    if any(len(set(domain)) != len(domain) for domain in domains):
        raise ValueError("domain labels must be distinct")
    representatives = []
    tuple_colors = []
    for obj, domain in enumerate(domains):
        colors = {}
        for points in permutations(domain, arity):
            color = next((i + 1 for i, (other_obj, other_points) in enumerate(representatives) if equivalent(other_obj, other_points, obj, points)), None)
            if color is None:
                representatives.append((obj, points))
                color = len(representatives)
            colors[points] = color
        tuple_colors.append(colors)
    return tuple(Configuration(len(domain), arity, tuple(colors.get(points, 0) for points in product(domain, repeat=arity))) for domain, colors in zip(domains, tuple_colors))


def nonfull_relations(left: CertificateTable, right: CertificateTable, left_domain: tuple[int, ...], right_domain: tuple[int, ...]) -> tuple[Configuration, Configuration]:
    """Build the small-support case's canonical relations from nonfull tests."""
    if left.engine is not right.engine or left.arity != right.arity:
        raise ValueError("tables must use the same certificate engine and arity")
    tables = (left, right)
    domains = (left_domain, right_domain)
    for table, domain in zip(tables, domains):
        for test_set in combinations(sorted(domain), table.arity):
            if test_set not in table.certificates or table.certificates[test_set].full:
                raise ValueError("every test set in the domain must have a nonfull certificate")
    action = left.engine.action
    cache = {}

    def equivalent(i, source_tuple, j, target_tuple):
        source_set, target_set = tuple(sorted(source_tuple)), tuple(sorted(target_tuple))
        key = (i, source_set, j, target_set)
        if key not in cache:
            result = left.engine.compare(tables[i].source, tables[j].source, tables[i].certificates[source_set], tables[j].certificates[target_set])
            cache[key] = None if result is None else Coset(action.restrict(result.subgroup).image, action(result.representative))
        projected = cache[key]
        if projected is None:
            return False
        representative_inverse = inverse(projected.representative)
        return projected.subgroup.transporter(source_tuple, tuple(representative_inverse[p] for p in target_tuple)) is not None

    return local_guide_relations(domains, left.arity, equivalent)

