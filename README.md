# Babai Graph Isomorphism in Quasipolynomial Time

A Python reference implementation of Babai's graph isomorphism framework for
finite simple undirected graphs. Returns an isomorphism and generators of the
automorphism group. Requires Python 3.11+ and only the standard library at runtime.

This research implementation has large constants and can be slow. Its
[handwritten analysis](docs/COMPLEXITY.md) derives a conservative
quasipolynomial time bound of `exp(O((log n)^5))` for graphs with `n` vertices.

## Quick Start

```bash
python3 -m pip install -e .
```

```python
from babai import Graph, graph_isomorphisms

g = Graph(5, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)])
h = g.relabel((2, 4, 1, 0, 3))
result = graph_isomorphisms(g, h)
assert result is not None and result.order == 10
print(result.representative)
print(result.subgroup.generators)
```

Vertices are numbered `0, ..., n-1`. The result describes all isomorphisms;
`result.order` is their count. `None` means the graphs are not isomorphic.

## Validation

```bash
python3 -m pip install -e '.[validation]'
python3 -m pytest -q
python3 -m scripts.verify_models
python3 -m scripts.verify_atlas --max-n 6
```

## Documentation

- [Algorithm and invariants](docs/ALGORITHM.md)
- [Correctness review and validation results](docs/REVIEW.md)
- Mathematical sources: [Babai's paper](https://arxiv.org/abs/1512.03547)
  and [Helfgott's corrected exposition](https://arxiv.org/html/1710.04574).

## AI Disclosure

OpenAI Codex (GPT-6-Astra) generated and revised code, tests, examples, and
documentation, including the complexity analysis. Automated checks and comparisons
against independent algorithms do not constitute independent expert review or
formal verification of the implementation or its claimed complexity bound.
