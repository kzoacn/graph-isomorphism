"""Independent atlas regression: python3 -m scripts.verify_atlas --max-n 6.

NetworkX is used only for the reference graph atlas and VF2 automorphism
counts, never by the implementation. The report records the tested code
digest so it cannot silently stand in for a later modified implementation.
"""

import argparse
from collections import defaultdict
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import platform
from random import Random
from time import perf_counter

import networkx as nx

from babai.graph import Graph, graph_isomorphisms


def code_digest():
    digest = sha256()
    root = Path(__file__).resolve().parents[1]
    for path in sorted((root / "babai").glob("*.py")):
        digest.update(path.name.encode())
        digest.update(path.read_bytes())
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-n", type=int, default=6, choices=range(8))
    parser.add_argument("--report", type=Path, default=Path("docs/verification-atlas.json"))
    args = parser.parse_args()
    rng = Random(17014574)
    started = perf_counter()
    graphs = [(i, g) for i, g in enumerate(nx.graph_atlas_g()) if len(g) <= args.max_n]
    positives = 0
    negatives = 0
    buckets = defaultdict(list)
    for atlas_index, original in graphs:
        n = len(original)
        graph = Graph(n, original.edges())
        p = list(range(n))
        rng.shuffle(p)
        renamed = graph.relabel(tuple(p))
        result = graph_isomorphisms(graph, renamed)
        expected = sum(1 for _ in nx.algorithms.isomorphism.GraphMatcher(original, original).isomorphisms_iter())
        if result is None or result.order != expected or not result.contains(tuple(p)):
            raise AssertionError(f"positive/counting failure at atlas entry {atlas_index}")
        if not graph.validates_isomorphism(renamed, result.representative):
            raise AssertionError(f"invalid witness at atlas entry {atlas_index}")
        positives += 1
        buckets[n, len(graph.edges)].append((atlas_index, graph, original))
        if positives % 100 == 0:
            print(f"verified {positives} relabelings and full automorphism counts", flush=True)
    for bucket in buckets.values():
        for (i, a, nx_a), (j, b, nx_b) in combinations(bucket, 2):
            # Check the reference, rather than inferring distinctness only
            # from atlas indices or its documented ordering.
            if nx.is_isomorphic(nx_a, nx_b):
                raise AssertionError(f"duplicate atlas isomorphism classes {i}, {j}")
            if graph_isomorphisms(a, b) is not None:
                raise AssertionError(f"false positive between atlas entries {i}, {j}")
            negatives += 1
            if negatives % 1000 == 0:
                print(f"verified {negatives} equal-order, equal-edge-count nonisomorphisms", flush=True)
    report = {
        "passed": True,
        "max_vertices": args.max_n,
        "positive_relabelings_and_full_counts": positives,
        "negative_pairs_equal_vertices_and_edges": negatives,
        "source": "https://networkx.org/documentation/stable/reference/generated/networkx.generators.atlas.graph_atlas_g.html",
        "networkx": nx.__version__,
        "python": platform.python_version(),
        "code_sha256": code_digest(),
        "elapsed_seconds": round(perf_counter() - started, 3),
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
