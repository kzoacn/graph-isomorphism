# Babai Graph Isomorphism: Python Reference Implementation

This project tests isomorphism of general finite simple undirected graphs and
returns an isomorphism together with generators of the automorphism group.
It uses Babai's Luks framework, local certificates, the Design Lemma, the
corrected Split-or-Johnson procedure, and recursion on group actions.
The implementation requires only the Python standard library at runtime.

To make the implementation inspectable, it uses explicit higher-arity tuples,
deterministic group operations, and generous enumeration thresholds. The
[handwritten complexity analysis](docs/COMPLEXITY.md) derives a conservative
bound of `exp(O((log n)^5))`, where `n` is the number of graph vertices.
This is a research reference implementation with large constants: highly
symmetric or larger instances can be slow.

## AI Disclosure

This project was developed with extensive assistance from OpenAI Codex.
AI generated and revised implementation code, tests, examples, and documentation,
including the complexity analysis and English translations.

Automated tests and comparisons against a separate implementation were run during
development. These checks do not constitute independent expert review or formal
verification of the implementation or its claimed complexity bound. The underlying
mathematical results are attributed to the research papers listed in References.

## Usage

Python 3.11 or later is required. Run the example from the project directory:

```bash
python3 -m examples.demo
```

Alternatively, install the project as a local package:

```bash
python3 -m pip install -e .
```

```python
from babai import Graph, graph_isomorphisms

g = Graph(5, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)])
h = g.relabel((2, 4, 1, 0, 3))
result = graph_isomorphisms(g, h)
if result is None:
    print("Not isomorphic")
else:
    print("An isomorphism:", result.representative)
    print("Automorphism group generators:", result.subgroup.generators)
    print("Number of isomorphisms:", result.order)  # 10 for a 5-cycle
    assert g.validates_isomorphism(h, result.representative)
```

Vertices are numbered `0, ..., n-1`. Self-loops and out-of-range vertices are
rejected; duplicate edges are merged. Empty graphs, isolated vertices, and
graphs with different numbers of vertices are handled.

The result `Coset(H, r)` describes all isomorphisms `H * r`, where `H = Aut(g)`.
Permutations act on the right: `compose(p, q)[v] == q[p[v]]`, meaning that `p`
is applied first and `q` second. `None` means nonisomorphism has been established.

## String Isomorphism and Recursion Records

```python
from babai import BabaiSolver, symmetric_group

solver = BabaiSolver()
result = solver.solve(symmetric_group(4), "aabb", "baba")
assert result is not None and result.order == 4

# Pass this solver to the graph interface to inspect model recursion events.
# Events record actual and ideal degrees, subproblem sizes, and branch counts.
print(solver.events)
```

`BabaiSolver(prefer_structural=True)` prioritizes structural reductions for
non-giant image groups, making those branches easier to study. The default
strategy may instead enumerate an image group whose order is provably small enough.

`LuksSolver` remains available as a separate component. Used on its own, it
requires an external continuation for large quotient groups. The default graph
interface uses `BabaiSolver`, which includes the main recursion for large groups.

## Tests and Independent Validation

```bash
python3 -m pip install -e '.[test]'
python3 -m pytest -q
```

Tests cover full isomorphism cosets, group homomorphisms and stabilizers, parity,
Johnson representations, local certificates and their aggregation, replayable
individualization choices, inherited block systems, and subset representations
across recursive levels. They include the rook/Shrikhande strongly regular graph
pair: ordinary 2-WL color statistics agree, but the solver correctly rejects
isomorphism.

The independent atlas validation uses the NetworkX graph atlas as input and VF2
as a reference algorithm:

```bash
python3 -m pip install -e '.[validation]'
python3 -m scripts.verify_atlas --max-n 6
```

Results are recorded in [verification-atlas.json](docs/verification-atlas.json),
including the SHA-256 digest of the tested code. This validation checks relabelings
and full automorphism counts for 209 graphs, plus 1,340 nonisomorphic pairs with
equal vertex and edge counts. NetworkX is not used by the implementation at runtime.

Additional subset-model and structural replay checks use an independent explicit
group-closure oracle:

```bash
python3 -m scripts.verify_models
```

See the [correctness and complexity review](docs/REVIEW.md) for findings,
reproduction commands, and limits of this validation.

## Implementation Map

| File | Contents |
| --- | --- |
| `babai/solver.py` | Unified main recursion, giant and non-giant image groups, certificate aggregation branches |
| `babai/graph.py` | Graph-to-string isomorphism reduction and lifting of results |
| `babai/string_iso.py`, `babai/tower.py` | Luks recursion, the window chain rule, and block towers |
| `babai/group.py`, `babai/action.py` | Deterministic Schreier-Sims, cosets, homomorphisms, and kernels |
| `babai/standard_blocks.py` | Standard blocks and complete subset-action models with uniform fibers |
| `babai/local_certificates.py`, `babai/aggregation.py` | Local certificates, comparisons, fullness structures, and local relations |
| `babai/configuration.py`, `babai/design.py` | Higher-arity WL, coherent configurations, and the Design Lemma |
| `babai/split_johnson.py`, `babai/coherent_split.py` | Split-or-Johnson with the 2017 correction |
| `babai/johnson.py` | Constructive recognition of Johnson schemes |
| `babai/discovery.py`, `babai/alignment.py` | Structure discovery, choice replay, and structural alignment |
| `babai/orbitals.py`, `babai/top_action.py` | Orbital configurations, transitivity reduction, and top-action lifting |

See the [algorithm audit](docs/ALGORITHM.md) for implementation invariants and
validation coverage.

## References

- L. Babai, [Graph Isomorphism in Quasipolynomial Time](https://arxiv.org/abs/1512.03547).
- H. A. Helfgott, J. Bajpai, D. Dona,
  [Graph Isomorphisms in quasi-polynomial time](https://arxiv.org/html/1710.04574).
- L. Babai, [2017 correction announcement](https://people.cs.uchicago.edu/~laci/update.html).
