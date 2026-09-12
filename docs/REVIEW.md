# Correctness and Complexity Review

Review date: 2026-09-13.

No solver correctness defect or counterexample was found in this audit. The
`exp(O((log n)^5))` upper bound is consistent with the reviewed recursion,
subject to the stated input-cost assumptions and the cited mathematical lemmas.
This is a code and mathematical-argument review by Codex, not independent expert
review or formal verification. Passing tests does not prove the asymptotic bound.

The reviewed algorithm files have SHA-256 digest
`15c0f489ea85a549ae7373d3098a91ba582a5ba00aaa03ababd22458476ae986`, using the
filename-and-content convention in `scripts/verify_atlas.py`. No file under
`babai/` was changed during the review.

## Findings and Documentation Changes

The previous analysis omitted several intermediate arguments needed to assess
the implementation's modifications to the published algorithm. They are now
explicit in [COMPLEXITY.md](COMPLEXITY.md):

- Section 2 proves the small structural kernel's `2W/3` orbit bound, including
  complemented supports, and explains why its action factors through the whole
  pure window even when some new supports are trivial. The Johnson support
  calculation now explicitly uses the threshold to obtain `t<m'/10<sqrt(c)`.
- Section 3 explains why non-full certificate relations satisfy the Design
  Lemma's twin-class hypothesis.
- Section 4 proves that the enlarged `64 ceil(log2 L)^2` cutoff supports the
  design-counting argument. It accounts for terminal choices, constant-factor
  reductions after other choices, and canonical reductions by just one point.
- Section 6 gives an explicit phase recurrence. It uses the implementation's
  decreasing window-size order, pairs the double-transitivity step with its
  next structural reduction, and counts recursive calls made by unsuccessful
  TopAction tests.

The general string input contract also needed clarification: comparability
alone is insufficient, because colors are dictionary keys. Hashability, stable
equality, input generator-list size, and symbol operation costs are now stated.
The graph interface's binary colors and initial generators satisfy these
conditions. README now calls the complexity derivation a handwritten analysis.

## Correctness Arguments Checked

The audit traced right-action composition, coset transport, window alignment,
nonfaithful action kernels, generator-defined homomorphisms, and coset merging.
It checked that structural restrictions are followed by actual string checks,
and that all images of recorded choices are included when completeness requires
them. Numerical orbital colors from different groups are matched through
individual relations rather than presumed to have common meanings.

The local-certificate threshold and short kernel windows agree with Theorem
8.3.5 and Corollary 8.3.7; the non-giant transitivity bound agrees with Theorem
13.1.1. These are from [Babai's paper](https://arxiv.org/abs/1512.03547).
The coherent primitive-right case was checked against Proposition 5.8 in the
PDF version of [Helfgott's corrected exposition](https://arxiv.org/html/1710.04574#S5.SS2).

For the complexity calculation, the reviewed bounds are:

| Work | Bound or progress |
| --- | --- |
| Group operations | Polynomial in degree and supplied generator-list length |
| Higher-arity configurations and local certificate tables | `N^O(log N)` local work and bounded short recursive calls |
| All choices and their matching in one structural stage | `exp(O((log N)^3))` |
| Continuing block or Johnson stages before actual-size descent | `O(log N)` along a branch |
| One phase's overhead and short subproblems | `A(N)=exp(O((log N)^4))` |
| Combined recurrence | `T(N)<=A(N)*(1+T(floor(2N/3)))` |

Summing the logarithmic phase costs gives `exp(O((log N)^5))`. Substituting
`N=binom(n,2)<=n^2` preserves this exponent in terms of graph vertices. This
derivation does not claim an optimized exponent or practical scalability.

## Executed Validation

| Check | Result |
| --- | --- |
| Existing pytest suite | 98 passed |
| NetworkX atlas/VF2, graphs through six vertices | 209 relabelings and complete automorphism counts passed |
| Atlas pairs with equal vertex and edge counts | 1,340 nonisomorphic pairs correctly rejected |
| Independent subset-model coset oracle | 240 comparisons passed: 152 nonempty and 88 empty cosets |
| Relabeling and choice replay on uniform hypergraphs | 84 passed |

The coset oracle closes generators by its own breadth-first traversal and
compares the entire returned set with the expected set. It does not use the
implementation's group enumeration or membership algorithms for this comparison.
Fixtures include symmetric, alternating, cyclic, dihedral, and intransitive
groups; subset ranks one through three; the balanced `t=m/2` case; and a
nontrivial fiber kernel. Both solver strategies are checked without changing
their thresholds. The separate replay checks lower only the small-right cutoff
to exercise deeper reductions on tractable examples.

Reproduction commands:

```bash
python3 -m pytest -q
python3 -m scripts.verify_atlas --max-n 6 --report /tmp/babai-review-atlas.json
python3 -m scripts.verify_models --report /tmp/babai-review-models.json
```

The atlas rerun reproduced the counts and code digest in
[verification-atlas.json](verification-atlas.json). The additional checks are
recorded in [verification-models.json](verification-models.json), including the
seed and digests of both the algorithm and verifier.

## Limits of the Review

There was no end-to-end aggregation run at the unmodified large-parameter
condition `10*max(9, N.bit_length()+2)<m`. Some existing aggregation and
combinatorial tests use lowered thresholds or directly constructed branch
fixtures. These establish finite-case behavior of those branches, not the
theorems governing arbitrarily large instances.

Validating a returned graph witness and its automorphism generators checks
soundness of that output. It cannot by itself detect a false negative or an
incomplete coset; the independent finite oracles provide that additional check
only for their tested inputs. No empirical timings or sampled recursion events
are used as evidence of the worst-case time bound.
