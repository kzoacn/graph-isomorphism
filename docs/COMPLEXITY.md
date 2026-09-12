# Correctness and a Conservative Quasipolynomial Bound

Let `N` be the actual string-domain size, `m` the current ideal-domain size,
and `n` the number of vertices in the original graph. The graph reduction uses
unordered vertex pairs, so `N = binom(n,2)`. The graph interface uses `0/1` colors.
The general string API requires hashable symbols with stable equality semantics:
it uses `Counter` and dictionaries. A time bound in the domain size alone assumes
polynomial-size symbols with polynomial-time hashing and equality, and a
polynomial-size input generating set. Otherwise their input size and operation
costs must also be included. The graph reduction satisfies these assumptions.

This analysis accounts for the more generous thresholds used in the code. It does
not claim the optimized exponent of the paper. The group-theoretic and
combinatorial results on which it relies are listed at the end. This is a
handwritten implementation analysis, not a machine-checked proof. The
[correctness and complexity review](REVIEW.md) records its audit and test coverage.

## 1. The Subset-Action Model

The main recursion maintains a surjection `phi: G -> P <= S_m` and an equivariant map

```
s : Omega -> { T subset Gamma : |T| = t },    1 <= t <= m/2.
```

Each `t`-subset has exactly `b` preimages, so `N = b * binom(m,t)`. The code checks
coverage, uniformity, and equivariance on generators. It does not require `P` to
be giant or `G` to act transitively on `Omega`.

The inequalities `binom(m,t) >= m` and `binom(m,t) >= 2^t` give

```
m <= N,                 t <= log2 N.
```

The kernel `ker(phi)` fixes every support subset, so its actual orbits have size
at most `b <= N/m <= N/2`. The action on proper subsets is faithful: an ideal
permutation fixing every `t`-subset also fixes every ideal point. This establishes
that subsequent model actions factor through the corresponding actual windows.

The initial Luks quotient action uses `t=1`. The graph entry point uses unordered
vertex pairs directly, taking complements when needed. Cases with `n<3` are
handled separately by preimages, retaining the kernel of the potentially
nonfaithful action on vertex pairs.

## 2. Structural Reductions Preserve the Model

For ideal color classes `C_1,...,C_r`, partition the actual domain by the vector
`(|s(v) intersect C_i|)_i`. Babai's Lemma 5.2.1 states that every window has size
at most `2N/3`, except possibly a single pure window

```
Omega(C) = {v : s(v) subset C}
```

A pure window larger than `2N/3` requires `|C| > 2m/3`. It still consists of
uniform fibers over all `t`-subsets of `C`.

There are two structural reductions on a dominant color class:

- An equipartition: the new ideal points are its blocks, and an actual point's
  support becomes the set of blocks touched by its old support.
- A Johnson embedding: the vertices of the color class are `k`-subsets of a new
  set, and an actual point's support becomes the union of those `k`-subsets.

Partition the actual domain by the new support size `u`. Each size that occurs
covers every `u`-subset with uniform fibers. The full symmetric group on the new
ideal domain induces bijections preserving the entire combinatorial point set
and its counts. This counting argument does not require the current search group
to contain that full symmetric group. The code checks the resulting property
again through `subset_action`.

The new ideal domain is at most half as large as the old one. Blocks of an
equipartition have size at least two. A Johnson embedding satisfies
`|C| = binom(m',k) >= binom(m',2)` and `m'>=5`, giving `m'<=|C|/2`.

For small `m'`, the new image group can be enumerated directly. Write `c=|C|`
and `W=b*binom(c,t)` for the pure-window size. In the partition case, the kernel
fixes every block, so each of its orbits lies in one block-count profile. For
`t<=c/2`, mixed profiles have size at most `2W/3` by Lemma 5.2.1. A profile wholly
inside one of the `m'>=2` equally sized blocks occupies at most
`W*(1/m')^t<=W/2`. If `t>c/2`, the complement bijection inside `C` gives the
same conclusion with `c-t`. In the Johnson case, the kernel fixes every old
Johnson vertex and hence every old support fiber, giving orbit size at most
`b<=W/c`. Here `0<t<c`: dominance and `t<=m/2` imply `t<3c/4`.

These small-parameter calls use the entire pure window before partitioning by
new support size. Fixing this window pointwise fixes every old `t`-subset of
`C`, hence every old vertex of `C`, and therefore the new blocks or Johnson
atoms. This proves that the auxiliary action factors through the window even
when some new union supports are the entire ideal domain. Such supports are
not passed to `subset_action`.

The small-parameter threshold here uses the pure-window size `W=|Omega(C)|`.
If `m' > 10*max(9, floor(log2 W)+3)`, the new support cannot be the entire ideal
domain. Write `c=|C|`. Dominance gives `t<3c/4`; if `t>c/2`, then
`W>=binom(c,t)>=2^(c/4)` and thus `m'<=c<=4 log2 W`, contradicting the threshold.
Consequently, the large-parameter branch has `t<=c/2`, and hence `t<=log2 W`.
The partition case has `u<=t`. In the Johnson case, the threshold also gives
`t<=log2 W<m'/10`, and `c=binom(m',k)>=binom(m',2)` implies
`m'/10<sqrt(c)`. Moreover, `k<=m'/2` gives `c>=2^k`. Consequently,
`log2 W >= t log2(c/t) >= (t/2) log2 c >= tk/2`.
Therefore `u<=tk<=2 log2 W<m'/2`.
The small-image branch covers the remaining small parameters. When needed,
nontrivial supports are complemented to normalize their size to at most half
of the ideal domain.

## 3. Each Branch Returns All Isomorphisms

The Luks window chain rule and small-quotient enumeration are exact coset
decompositions. When constant fibers are compressed, the kernel of the block
action is retained in the returned preimage.

In the giant case, isomorphisms of partitions and Johnson structures can be
constructed directly. `align_splits` finds structural isomorphisms and enforces
all selected point correspondences; `preimage_coset` lifts them to the actual
domain. The subsequent chain rule checks all actual coordinates.

The non-giant case uses orbits, blocks, and orbital configurations obtained from
the current group itself. Noncanonical choices are handled by taking their point
stabilizer and enumerating all its cosets. The code checks that this stabilizer
is contained in the structure's automorphism group, avoiding a general subgroup
intersection.

Local certificates provide global automorphisms or evidence of local
non-fullness. Aggregation uses equivalence classes of local isomorphisms across
the two inputs, so every actual global isomorphism preserves these colors. For a
large non-giant certificate orbit, one source relation is compared against every
possible target orbital relation.

In the small-support branch, every test set in the complement of the full
certificate group's support is non-full: a full test would move every point of
that test set in the projected group. The local-guide relation cannot have a
twin class of size at least `k`. Otherwise, for a test set inside that class,
all its ordered tuples would have the same color; the local comparison group
would induce every permutation of the test set, contradicting non-fullness.
Thus the Design Lemma's twin hypothesis holds when `10k<m` and the selected
domain has more than `2m/3` points.

TopAction generates the full coset from isomorphism fibers of image-group
generators. A failed test establishes only a failure of surjectivity. The
alternating automorphism projection is verified first, allowing a subsequent
failure to establish that a target branch is empty.

All possible images of noncanonical choices are covered, so the final union of
branches is exactly `Iso_G(x,y)`. When nonempty, the closure of its coset generators
stays within this isomorphism set and generates the entire coset.

## 4. Cost of Individual Subroutines

Every residual generator added by deterministic Schreier-Sims enlarges a basic
orbit, with at most `D(D-1)/2` such enlargements for an action of degree `D`.
Joint actions used for homomorphisms have degree `O(N)`. Group operations take
polynomial time in the degree and the supplied generator-list length, and
output generating sets are compressed to polynomial size. Homomorphisms, kernels,
lifts, and point stabilizers do not enumerate large groups.

F2 and higher-arity WL take `m^O(k)` time. All arities used in the main recursion
are `O(log N)`: they arise from the model's `t`, the locality parameter `k`, or
the design arity in the bipartite reduction. These steps therefore take at most
`N^O(log N)` time.

The small-right threshold in Split-or-Johnson is `64 ceil(log2 L)^2`, where `L`
is the left-side size. This case fixes every right point or fiber, using
`O(log^2 N)` original ideal points. Other steps that incur individualization cost
shrink the right side by a constant factor and fix at most `O(log N)` points each.
The full choice record therefore has length `O(log^2 N)`. Steps that decrease the
right side by only one point are canonical and introduce no additional group index.

The bipartite structure has arity `d<=6 ceil(log_R L)` with `R<L`, so
`R^d<=L^6 R^6<=L^12`. Tuple construction, F2, WL, and twin checks remain polynomial
in `L`.

The enlarged cutoff also supplies the counting hypothesis used to build this
relation. Let `L'<=L` be the retained twin-free left class, `r<=R/2` its
neighborhood size, and `d=min(r,6 ceil(log_R L'))`. The complete hypergraph is
handled separately. If `d=r`, a missing neighborhood distinguishes the tuple
colors. Otherwise put `s=d/2=3 ceil(log_R L')`. If every injective tuple had the
same color, the incidence counts would define a nonempty `d`-design. The
[Ray-Chaudhuri-Wilson inequality as stated in Section 2.6 of the exposition](https://arxiv.org/html/1710.04574#S2.SS6)
would give `L'>=binom(R,s)`. But the active branch has
`R>64 ceil(log2 L)^2`, so `s<=6 log2 L'<sqrt(R)` and

```
binom(R,s) >= (R/s)^s > R^(s/2) >= (L')^(3/2) > L'.
```

This contradiction excludes an all-vertex twin class. A proper large twin
class can cause a canonical reduction by just one right point; it adds no
individualization. The other continuing reductions after a choice shrink the
right side to at most `2R/3` or `R/2`. A choice that immediately returns a
partition is a terminal stage and needs no further size decrease.

All noncanonical matches and computation in one structure-discovery stage are
bounded by

```
Q(N) = N^O((log N)^2) = exp(O((log N)^3))
```

Naive coset traversal squares the index cost but stays within the same form of
bound. Binary WL pruning after individualization adds only polynomial work and
no additional choices.

## 5. Actual-Domain Descent for Local Certificates and Small Images

The default is `k=max(9, floor(log2 N)+3)`, and aggregation is used only when
`10k<m`. This satisfies the strict threshold of the Unaffected Stabilizer Theorem.
Each general window update enumerates at most `k!` test images and recurses on
kernel orbits. The Affected Orbit Lemma bounds every recursive window by `N/k`.

A certificate expands its window at most `N` times. There are at most `m^k` test
sets and `m^(2k)` pairs of test sets. Including all comparisons, the total number
of recursive calls is still `N^O(log N)`, with subproblems of size at most `N/k`.
Shortcuts for natural symmetric and alternating actions may inspect larger
windows, but solve them directly in polynomial time without invoking difficult
recursion.

The small-image branch satisfies one of the following conditions:

- `m=O(log N)`, allowing enumeration of at most
  `m! = exp(O(log N log log N))` image elements.
- The image-group order is at most `m^(1+ceil(log2 m))`.

The kernel fixes model fibers, so actual subproblems have size at most
`N/m<=N/2`. Small-parameter structural cases with trivial union supports use the
`2/3` window bound from Section 2. These conditions are checked at the call sites;
`m!` cannot be used arbitrarily as an enumeration bound for a general large group.

## 6. The Combined Recurrence

Fix a phase's starting actual degree `N`, and expand recursive calls until all
remaining actual subproblems have size at most `2N/3`. Actual orbit chains are
additive: their windows are disjoint and their total length does not increase.
Before entering a model they can require polynomially many steps, including
successive descents from degree `q` to `q-1`, with polynomial overhead.

Inside a model, at most one color window, and then at most one union-size
window, can exceed `2N/3`. The code processes windows in decreasing size order.
Consequently, any such continuing window is processed before the subgroup can
be changed by solving a different window. This matters for the transitivity
argument below. Other windows already end the phase.

Each continuing block or Johnson stage at least halves the ideal domain.
An ideal coloring without a dominant class ends the phase by the induced
`2/3` window bound. The step that must be analyzed together with its successor
is individualization in a non-giant
doubly transitive group. The Wielandt bound requires `O(log m)` points. The action
on the remaining large orbit is transitive but no longer doubly transitive, so
the next step is block reduction or Split-or-Johnson on an orbital configuration.
If the parameters are already small, enumeration proceeds through short kernel
orbits. Canonical ideal-orbit descent introduces no matching index. Its image
is transitive after restriction to the selected large orbit; a subsequent
complete-subset window has the same image because the ideal action factors
faithfully through that window. Thus orbit descent cannot repeatedly restart
the costly double-transitivity step without an intervening structural decrease.

On a continuing path there are therefore `O(log N)` charged stages. The calls
that start a fresh Luks problem or compress constant fibers already decrease
the actual degree by a constant factor, apart from natural symmetric or
alternating actions which are solved directly. In particular, the ideal degree
is not reset to a larger value along a continuing path within this phase.

The total multiplier from expanding this part of the model recursion is at most

```
Q(N)^O(log N) = exp(O((log N)^4)).
```

Certificate construction and unsuccessful TopAction tests may do recursive
work before a structural outcome is known. Their recursive windows already
have size at most `N/k` or a model fiber; the polynomial natural-action
shortcuts require no difficult recursion. Count these calls as additional
leaves of the phase, rather than assuming that every TopAction test terminates
the whole problem.

Including those leaves, bounded image enumeration, and the disjoint short
windows, both the phase overhead and its number of leaves are at most
`A(N)=exp(O((log N)^4))`. For the nondecreasing worst-case envelope `T`, this gives

```
T(N) <= A(N) * (1 + T(floor(2N/3)))
```

outside a fixed finite base range. Iterating through `O(log N)` phases bounds
`log T(N)` by a constant times
`sum_j (log N - j*log(3/2))^4 = O((log N)^5)`, and hence

```
T(N) <= exp(O((log N)^5)).
```

Purely additive orbit chains, even when an orbit has size `N-1`, add only
polynomial work and sums whose total subdomain length does not increase. They
do not introduce the exponential multiplier above at every level. Polynomial
costs of large integers and complete signatures are absorbed by the bound; the
analysis does not assume collision-free hashing.

Finally, `N=binom(n,2)<=n^2`, so the graph interface still has the bound
`exp(O((log n)^5))`, which is quasipolynomial.

## Sources

- [Babai's paper](https://arxiv.org/abs/1512.03547): the Luks framework,
  Lemma 5.2.1, Unaffected Stabilizer and Affected Orbit results, TopAction,
  local certificates, and the orbital-configuration route in Section 13.
- [The exposition by Helfgott et al.](https://arxiv.org/html/1710.04574): the
  Design Lemma, corrected Split-or-Johnson in Section 5.2, certificate aggregation
  in Section 6, and Appendix A.
- [The 2017 correction](https://people.cs.uchicago.edu/~laci/update.html): the
  primitive-right branch uses a single recursive instance with its size halved.

Uniform subset models explicitly track actions before and after structural
reductions. The enlarged small-right threshold simplifies the implementation.
The derivation above includes the additional costs of these implementation choices.
