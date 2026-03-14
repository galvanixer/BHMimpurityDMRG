# Pair-Distance Distribution Evaluation Philosophy

This note explains how the code evaluates the pair-distance distribution, what is
stored, and the design choices behind those decisions.

The implementation lives in:

- `src/observables/pair_distance_distribution.jl` (core evaluators)
- `src/observables/pipeline.jl` (orchestration and HDF5 writing)

The canonical HDF5 output lives under:

- `/observables/pair_distance_distribution`

## What We Measure

For one species, the observable is a pair-count profile `C(r)` and its normalized
distribution `P(r)`.

At `r=0`, the code uses on-site unordered boson pair counts:

```math
C(0)=\sum_i \frac{1}{2}\langle n_i(n_i-1)\rangle.
```

For `r>0`, it sums off-site unordered pairs:

```math
C(r)=\sum_{i<j,\;d(i,j)=r}\langle n_i n_j\rangle.
```

The normalized distribution is:

```math
P(r)=\frac{C(r)}{\sum_{r'} C(r')}.
```

For cross species (`a`,`b`), the profile is:

```math
C^{(ab)}(r)=\sum_{i,j,\;d(i,j)=r}\langle n_i^{(a)}n_j^{(b)}\rangle,
```

with:

```math
P^{(ab)}(r)=\frac{C^{(ab)}(r)}{\sum_{r'} C^{(ab)}(r')}.
```

The pipeline also stores:

- `pair_count_* = \sum_r C(r)` (normalization)
- `mean_r_* = \sum_r r P(r)`
- `var_r_* = \sum_r r^2 P(r) - mean_r_*^2`
- `std_r_* = sqrt(var_r_*)`

where `*` is `a`, `b`, or `ab`.

## Stored Datasets

When requested, the group stores metadata:

- `same_site_convention`
- `cross_species`

Species-resolved arrays/scalars are written when available:

- `r_a`, `C_a`, `P_a`, `pair_count_a`, `mean_r_a`, `var_r_a`, `std_r_a`
- `r_b`, `C_b`, `P_b`, `pair_count_b`, `mean_r_b`, `var_r_b`, `std_r_b`
- `r_ab`, `C_ab`, `P_ab`, `pair_count_ab`, `mean_r_ab`, `var_r_ab`, `std_r_ab` (if `cross_species=true`)

## Philosophy

### 1. Count pairs, do not average correlators

The primary object is `C(r)`, a pair-count measure, not a translational average of
`⟨n_i n_{i+r}⟩`.

This is deliberate. A distribution should represent "how much pair weight sits at each
distance," so we accumulate pair weight directly and normalize once at the end.

### 2. Use unordered counting for same species

For same-species observables, pairs are physically unordered.

The code reflects that by:

- using `i<j` off-site,
- and including the on-site combinatorial factor `1/2` in `C(0)`.

This keeps normalization consistent with the expected number of unordered particle pairs.

### 3. Make diagonal convention explicit

For same-species diagonal terms, users may store either:

- `⟨n_i(n_i-1)⟩` (`"factorial"`)
- or `⟨n_i^2⟩` (`"plain"`).

Because this changes `C(0)`, the chosen `same_site_convention` is propagated and stored.
Downstream analysis should never need to guess which diagonal semantics were used.

### 4. Use minimum-image distances for periodic systems

For periodic lattices, the code bins by unsigned minimum-image distance:

```math
d(i,j)=\min(|i-j|,L-|i-j|).
```

This is the natural notion of scalar separation on a ring and avoids splitting the same
physical separation into two bins (`r` and `L-r`).

### 5. Keep species-resolved structure first, cross-species optional

The pipeline treats `a` and `b` distributions as independent primary observables and
enables `ab` through `cross_species`.

This mirrors the physical workflow in impurity problems:

- `P_a(r)` and `P_b(r)` answer different questions,
- `P_ab(r)` is valuable but has distinct interpretation and normalization.

### 6. Store both raw and normalized forms

We store both `C(r)` and `P(r)`.

`P(r)` is convenient for comparison across runs, while `C(r)` preserves absolute pair
weight. Keeping both avoids forcing users to reverse-engineer one from the other.

### 7. Store compact moments for campaign-scale work

The pipeline stores `mean_r`, `var_r`, and `std_r` to provide low-dimensional summaries
for scans and ranking tasks.

The full distribution remains available when detailed interpretation is needed.

### 8. Reuse two-point primitives instead of re-measuring

This observable is built from already-computed one-point and two-point density objects.
That keeps semantics aligned with `density_density` and avoids redundant expensive
measurements.

The design goal is consistency first: all pair-based observables should agree on the
underlying correlators and conventions.

## Consistency Checks We Care About

1. Normalization:

```math
\sum_r P(r)=1
```

when `pair_count > 0`.

2. Non-negativity (within numerical tolerance):

```math
C(r)\ge 0,\quad P(r)\ge 0.
```

3. Periodic distance support:

- periodic: `r = 0,1,...,floor(L/2)`
- open: `r = 0,1,...,L-1`

4. Same-species combinatorics:

`C(0)` should track local bunching via `1/2 * <n_i(n_i-1)>` under the selected
diagonal convention.

5. Cross-species normalization:

`pair_count_ab` should match total expected cross pairs (`\sum_{i,j}\langle n_i^{(a)}n_j^{(b)}\rangle`)
up to numerical tolerance.

## Downstream View

The pair-distance distribution is intended to diagnose binding/size structure with a
clear probabilistic interpretation:

- weight concentrated near `r=0` suggests compact pairing,
- broad `P(r)` indicates extended/unbound behavior,
- peaks at finite `r>0` indicate preferred nonzero internal separation.

The guiding philosophy is:

- use a pair-count definition with explicit combinatorics,
- preserve species structure,
- enforce transparent normalization,
- and keep both rich (`C(r),P(r)`) and compact (`mean_r,var_r,std_r`) outputs.
