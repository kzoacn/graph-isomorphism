# Implementation Audit

## Scope and Entry Points

The default entry point, `graph_isomorphisms`, handles general finite simple
undirected graphs and returns the full isomorphism coset. The general string
isomorphism entry point is `BabaiSolver.solve`, which proceeds from Luks reductions
to the complete subset-model recursion.

The implementation follows the orbital-configuration route in Section 13 of
Babai's paper to reduce non-giant image groups, avoiding a dependency on an
unimplemented Cameron group recognizer. Local certificates use the paper's
Unaffected Stabilizer Theorem, and the combinatorial recursion incorporates the
2017 correction. See [COMPLEXITY.md](COMPLEXITY.md) for the worst-case analysis
and the more generous thresholds used in this implementation.

## Key Invariants

### Permutations and Cosets

`compose(p, q)[v] = q[p[v]]`. String isomorphisms satisfy `x[v] == y[g[v]]`.
`Coset(H, r)` denotes `{compose(h, r): h in H}`.

The window chain rule uses `y[r[v]]` as the aligned target. It first solves the
restricted action, then takes the homomorphism preimage to restore the full
domain and retain the kernel. `merge_isomorphism_cosets` requires that the union
of branches is known to be an isomorphism coset. `coset_hull` explicitly constructs
the smallest containing coset and is used for the generator fibers in TopAction;
the union of those fibers need not itself be a coset.

### Group Operations

Schreier-Sims constructs a deterministic stabilizer chain without enumerating the
entire group. The group order is the product of the basic orbit sizes, and
membership is tested by sifting. Group objects are immutable so that changing
generators cannot leave a stale stabilizer chain in use.

Homomorphisms use the joint action of the original and image representations,
checking that the generator images respect the group relations. Stabilizer chains
of the joint group provide kernels and constructive lifts. `preimage(H)` accepts
only a subgroup of the image; it does not hide a general subgroup-intersection
operation. Every enumeration of group elements or cosets requires an explicit bound.

### Block Systems and Subset Models

`BlockTower` retains existing block systems. New Luks block systems must coarsen
the existing tower, and the whole tower is restricted when descending to a window.
A fresh initial tower is created for a new representation only when complete
constant fibers are compressed and the actual degree decreases by at least half.

`subset_action` verifies that the actual domain consists of uniform fibers over
all `t`-subsets of the ideal domain:

- Each actual point has a nonempty proper support subset. Complements are taken
  when necessary to ensure `t <= m/2`.
- Every `t`-subset occurs, with the same fiber size.
- The support map commutes with every group generator and its image.

This model does not require the current group to be transitive or its image to be
giant. After a structural reduction, the new support is the set of blocks touched
or the union of Johnson subsets. The actual domain is partitioned by support size.
Each new model is checked again for completeness, uniformity, and equivariance.
`descend_action` checks that the ideal action factors through the restriction to
the actual window.

### Configurations and Noncanonical Choices

F2 retains equality patterns and all coordinate maps, including noninjective maps.
Higher-arity WL uses the joint color distribution for replacements at the same
point in every coordinate. Signatures are compared by their complete values;
hash collisions do not change equality decisions.

Local relations use shared color meanings. Every noncanonical choice in the
Design Lemma and Split-or-Johnson is recorded in `Choice`. The main recursion
enumerates its possible images and replays the choices on the target structure.
The corresponding point-mapping constraints are explicitly enforced after
structural alignment. Binary WL after individualization only prunes this already
bounded choice family; it introduces no new branches.

Orbital color integers from different certificate groups do not automatically
have shared meanings. `_nongiant_certificate_orbit` fixes one source orbital
relation and compares it with every possible target orbital relation.

### Johnson and Split-or-Johnson

Johnson recognition reconstructs the underlying set and then verifies every pair
color; intersection numbers alone are insufficient. For a complete uniform
hypergraph with `m = 2k`, complementary pairs give blocks of size two.

The corrected coherent primitive-right branch produces a single recursive
instance whose right side is at most half as large. Other combinatorial steps
that introduce noncanonical choices also shrink the right side by a constant
factor. Only canonical steps that add no choices may decrease it by just one
point. The enlarged small-right threshold and its index cost are included in
the complexity analysis.

### Local Certificates and TopAction

The default locality parameter is `max(9, N.bit_length() + 2)`, satisfying the
strict Unaffected Stabilizer threshold. Certificate aggregation is used only when
`10*k < m`; other cases fall under enumeration of a small ideal action. General
local updates recurse on windows of size at most `floor(N/k)`. Natural symmetric
and alternating actions can be solved directly in polynomial time.

A full certificate fixes the complement of its window. Its generators are checked
to be global string automorphisms, and its image must contain the alternating
group on the test set. Certificate comparison replays both windows together and
can impose an ordered test-tuple constraint. Local-guide colors are equivalence
classes of local isomorphisms across the two inputs.

TopAction checks that image-group generators lift to actual isomorphisms. Failure
to lift a generator establishes only a failure of surjectivity, not nonisomorphism.
The source automorphism projection must first be verified to contain the
alternating group before a subsequent failure can establish nonisomorphism.

## Main Recursion Branches

| Branch | Progress or result |
| --- | --- |
| Natural `S_m` / `A_m` action | Solve by color classes and handle parity |
| Intransitive Luks action | Apply the window chain rule without increasing total subdomain size |
| Small quotient or ideal domain | Enumerate a bounded image group; kernel orbits strictly reduce the actual domain |
| Constant fibers | At least halve the actual domain and retain the full kernel in the preimage |
| Intransitive ideal action | Color by orbits; descent to a large orbit introduces no exponential branching |
| Imprimitive ideal action | The block action at least halves the ideal domain |
| Non-giant, doubly transitive action | Individualize within the Wielandt bound, then reduce an action that is no longer doubly transitive |
| Non-giant, uniprimitive action | Apply Split-or-Johnson to the orbital configuration |
| Orbit partition from full certificates | Align and pull back to actual windows and a smaller block action |
| Large alternating symmetry | Use TopAction and lifting on short windows |
| Large non-giant certificate orbit | Reduce transitivity, match orbital relations, and apply Split-or-Johnson |
| Small support of full certificates | Build non-full local relations, apply joint WL, the Design Lemma, and Split-or-Johnson |

All of these branches are implemented in the default `BabaiSolver`.
`GiantActionRequired` remains only in the standalone `LuksSolver` component to
indicate that its caller has not installed a continuation. The default graph
interface uses `BabaiSolver`.

## Validation Coverage

Tests check witnesses, full isomorphism counts, and generators. Small groups and
graphs are checked against independent enumeration; the atlas script uses the
NetworkX graph atlas and VF2 as references. Coverage includes nontrivial kernels,
block and Johnson union supports, parity, choice replay, shared local-guide
colors, the large non-giant orbit of the Petersen automorphism group, and the
rook/Shrikhande pair that ordinary 2-WL cannot directly distinguish.

A few tests change shortcut thresholds to execute deeper combinatorial branches
on small instances. Small-instance aggregation tests lower the test-set threshold
but still verify the global automorphism and image-group conditions of the full
certificates they construct. These tests exercise code paths; they are not
evidence for the original theorem thresholds or asymptotic bounds.

`docs/verification-atlas.json` records the digest of the tested code so that an
old run cannot be mistaken for validation of new code. Tests check implementation
errors; the worst-case bound comes from the separate recurrence analysis.

The [2026-09-13 review](REVIEW.md) checks the arguments against the implementation
and the cited papers, and records additional independent finite checks.
