"""Independent finite checks of subset-model cosets and structural choice replay.

Run with ``python3 -m scripts.verify_models``. The coset oracle uses its own
breadth-first closure of permutation generators, not the solver's enumerator,
membership tests, or stabilizer chain. These finite checks do not establish an
asymptotic bound. Only the replay checks lower the small-right cutoff.
"""

import argparse
from collections import deque
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
from random import Random
from time import perf_counter
from unittest.mock import patch

from babai import BabaiSolver, PermutationGroup
from babai.action import Homomorphism
from babai.coherent_split import BipartiteGraph
from babai.permutation import inverse, pullback
from babai.split_johnson import Choice, bipartite_split
from babai.standard_blocks import subset_action


def closure(generators, degree):
    identity = tuple(range(degree))
    seen = {identity}
    pending = deque([identity])
    while pending:
        p = pending.popleft()
        for g in generators:
            q = tuple(g[p[v]] for v in range(degree))
            if q not in seen:
                seen.add(q)
                pending.append(q)
    return tuple(sorted(seen))


def make_model(m, rank, fibers, kind):
    cycle = tuple((v + 1) % m for v in range(m))
    if kind == "symmetric":
        generators = (cycle, (1, 0) + tuple(range(2, m)))
    elif kind == "alternating":
        generators = tuple(tuple(1 if v == 0 else j if v == 1 else 0 if v == j else v for v in range(m)) for j in range(2, m))
    elif kind == "cyclic":
        generators = (cycle,)
    elif kind == "dihedral":
        generators = (cycle, tuple((-v) % m for v in range(m)))
    elif kind == "intransitive":
        generators = ((1, 2, 0) + tuple(range(3, m)),)
    else:
        raise ValueError("unknown fixture group")
    subsets = tuple(map(frozenset, combinations(range(m), rank)))
    positions = {s: i for i, s in enumerate(subsets)}
    images = {}
    for g in generators:
        induced = tuple(fibers * positions[frozenset(g[v] for v in s)] + f for s in subsets for f in range(fibers))
        images[induced] = g
    degree = fibers * len(subsets)
    if fibers == 2:
        images[(1, 0) + tuple(range(2, degree))] = tuple(range(m))
    group = PermutationGroup(degree, images)
    action = Homomorphism(group, m, (images[g] for g in group.generators))
    model = subset_action(action, tuple(s for s in subsets for _ in range(fibers)))
    ambient = closure(group.generators, degree)
    if len(ambient) != group.order:
        raise AssertionError("group order disagrees with independent closure")
    return model, ambient


def verify_cosets(seed):
    rng = Random(seed)
    checks, positive, negative = 0, 0, 0
    branches = set()
    for m, rank, fibers in ((4, 1, 1), (5, 2, 1), (6, 3, 1), (4, 1, 2)):
        for kind in ("symmetric", "alternating", "cyclic", "dihedral", "intransitive"):
            model, ambient = make_model(m, rank, fibers, kind)
            degree = len(model.supports)
            for iteration in range(6):
                source = tuple(rng.randrange(2 + iteration % 2) for _ in range(degree))
                if iteration % 2:
                    target = list(source)
                    rng.shuffle(target)
                    target = tuple(target)
                else:
                    target = pullback(source, inverse(rng.choice(ambient)))
                expected = {g for g in ambient if all(source[v] == target[g[v]] for v in range(degree))}
                for structural in (False, True):
                    solver = BabaiSolver(prefer_structural=structural)
                    result = solver.solve_model(model, source, target)
                    actual = set() if result is None else {
                        tuple(result.representative[h[v]] for v in range(degree))
                        for h in closure(result.subgroup.generators, degree)
                    }
                    if actual != expected:
                        raise AssertionError((m, rank, fibers, kind, iteration, structural, len(expected), len(actual)))
                    checks += 1
                    positive += bool(expected)
                    negative += not expected
                    branches.update(event.kind for event in solver.events)
    return {"comparisons": checks, "nonempty_cosets": positive, "empty_cosets": negative, "events": sorted(branches)}


def mapped_parts(result, permutation):
    return tuple(frozenset(frozenset(permutation[v] for v in block) for block in blocks) for blocks in result.parts)


def verify_replay(seed):
    rng = Random(seed)
    checks = 0
    for m in range(5, 9):
        for rank in range(2, min(3, m // 2) + 1):
            universe = list(combinations(range(m), rank))
            for iteration in range(12):
                kept = rng.sample(universe, rng.randrange(3 * m // 2 + 1, len(universe) + 1))
                left = tuple(range(len(kept)))
                right = tuple(frozenset((len(kept) + v,)) for v in range(m))
                graph = BipartiteGraph(left, right, tuple(map(frozenset, kept)))
                p = list(range(len(kept) + m))
                rng.shuffle(p)
                left_order, right_order = list(range(len(kept))), list(range(m))
                rng.shuffle(left_order)
                rng.shuffle(right_order)
                positions = {v: i for i, v in enumerate(right_order)}
                mapped = BipartiteGraph(
                    tuple(p[left[i]] for i in left_order),
                    tuple(frozenset(p[v] for v in right[i]) for i in right_order),
                    tuple(frozenset(positions[v] for v in kept[i]) for i in left_order),
                )
                with patch("babai.split_johnson._small_right_limit", lambda _: 0):
                    result = bipartite_split(graph)
                    plan = tuple(Choice(c.kind, tuple(p[v] for v in c.vertices)) for c in result.choices)
                    other = bipartite_split(mapped, plan=plan)
                if mapped_parts(result, p) != mapped_parts(other, tuple(range(len(p)))):
                    raise AssertionError((m, rank, iteration, "partition replay"))
                if (result.johnson is None) != (other.johnson is None):
                    raise AssertionError((m, rank, iteration, "Johnson presence"))
                if result.johnson is not None:
                    a, b = result.johnson, other.johnson
                    indices = {v: i for i, v in enumerate(b.vertices)}
                    for i, u in enumerate(a.vertices):
                        for j, v in enumerate(a.vertices):
                            if len(a.model.subsets[i] & a.model.subsets[j]) != len(b.model.subsets[indices[p[u]]] & b.model.subsets[indices[p[v]]]):
                                raise AssertionError((m, rank, iteration, "Johnson replay"))
                checks += 1
    return {"hypergraphs": checks, "small_right_cutoff_override": 0}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=20260913)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    started = perf_counter()
    report = {"seed": args.seed, "cosets": verify_cosets(args.seed), "replay": verify_replay(args.seed)}
    digest = sha256()
    for path in sorted((Path(__file__).resolve().parents[1] / "babai").glob("*.py")):
        digest.update(path.name.encode())
        digest.update(path.read_bytes())
    report.update(passed=True, code_sha256=digest.hexdigest(), verifier_sha256=sha256(Path(__file__).read_bytes()).hexdigest(), elapsed_seconds=round(perf_counter() - started, 3))
    output = json.dumps(report, indent=2) + "\n"
    if args.report is not None:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(output)
    print(output, end="")


if __name__ == "__main__":
    main()
