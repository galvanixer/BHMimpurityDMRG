# Density-Density Evaluation Philosophy

This note explains how the code evaluates density-density observables, what exactly is
stored, and the reasoning behind those choices.

The implementation lives primarily in `src/observables/density_density.jl`, with
orchestration and HDF5 writing in `src/observables/pipeline.jl`.

The canonical HDF5 output for each requested species currently lives under
`/observables/density_density`.

## What We Measure

For a fixed species, the basic two-point object is

```math
\mathrm{nnmat}_{ij} = \langle n_i n_j \rangle
```

away from the diagonal. The code currently evaluates this separately for species `a`
and species `b`, storing:

- `nn_a`, `connected_nn_a`, `variance_a`
- `nn_b`, `connected_nn_b`, `variance_b`
- `nn_ab`, `connected_nn_ab` when `cross_species=true`
- `r_a`, `transl_avg_nn_a`, `transl_avg_connected_nn_a`, `anchors_a`
- `r_b`, `transl_avg_nn_b`, `transl_avg_connected_nn_b`, `anchors_b`
- `r_ab`, `transl_avg_nn_ab`, `transl_avg_connected_nn_ab`, `anchors_ab` when `cross_species=true`

and, when requested, additional edge-trimmed translational profiles under
`/observables/density_density/profiles`.

The group also stores the metadata fields:

- `same_site_convention`
- `cross_species`
- `edge_trims`
- `default_edge_trim`

## Philosophy

### 1. Species-resolved observables come first, but mixed-species structure is explicit

The code treats `a`-boson and `b`-boson density correlations as distinct observables.

That is intentional. In this project the bath and impurity sectors generally play
different physical roles, so species-resolved observables are more informative than an
immediately combined density-density object.

The code can also store the mixed-species correlator

```math
\langle N^{(a)}_i N^{(b)}_j \rangle
```

as `nn_ab`, together with its connected counterpart and translational averages.

That is still kept explicit via the `cross_species` toggle rather than being merged into
the `a` and `b` outputs silently, because the mixed-species matrix has different symmetry
properties and a different physical interpretation.

### 2. The full matrix is the primitive object

The primary stored object is the full site-resolved matrix `nnmat`.

That matters because this code is aimed at impurity problems and other inhomogeneous
settings, where translational symmetry is often weak, broken, or simply not the main
organizing principle. In those cases the interesting questions are local:

- where correlations are strongest,
- how they change near the impurity,
- how much edges contaminate the apparent bulk behavior,
- whether the state is genuinely homogeneous or only approximately so.

If we stored only displacement-space averages, we would lose exactly that information.

So the philosophy is:

- store the primitive real-space matrix first,
- derive compressed summaries from it,
- and only average when there is a good physical reason to do so.

### 3. Connected correlations are important enough to store explicitly

For density observables, subtracting the disconnected background is often physically
meaningful:

```math
C_{ij} = \langle n_i n_j \rangle - \langle n_i \rangle \langle n_j \rangle.
```

Unlike the SPDM case, the connected form is not just an optional postprocessing trick.
It is often the quantity that isolates genuine fluctuations from trivial density
modulation.

That is why the code stores both:

- the raw matrix `nn_a` / `nn_b`
- the connected matrix `connected_nn_a` / `connected_nn_b`
- and, when requested, the mixed-species pair `nn_ab` / `connected_nn_ab`

The idea is to avoid forcing downstream analysis to choose between "only raw" and
"only connected" after information has already been discarded.

### 4. The same-site diagonal must be treated explicitly

For bosons, the diagonal `i=j` is not a trivial bookkeeping detail. There are two natural
quantities one may want on the diagonal:

```math
\langle n_i^2 \rangle
```

or

```math
\langle n_i(n_i-1) \rangle.
```

They answer different questions.

- `⟨n_i^2⟩` is the ordinary second moment.
- `⟨n_i(n_i-1)⟩` counts ordered local two-particle occupancy in the bosonic sense and is
  the more natural same-site continuation of pair-correlation language.

The code therefore makes the convention explicit via `same_site_convention`, rather than
hiding it.

For the mixed-species diagonal there is no corresponding ambiguity:

```math
\mathrm{nn}^{(ab)}_{ii} = \langle N^{(a)}_i N^{(b)}_i \rangle.
```

So the same-site convention applies to `aa` and `bb`, but not to `ab`.

### 5. The default is `"factorial"` for physical continuity of pair correlations

The default same-site convention is

```math
\mathrm{nnmat}_{ii} = \langle n_i(n_i-1) \rangle.
```

This is a deliberate physics choice.

For bosonic pair observables, the factorial form is the more natural "pair-counting"
object on-site. It avoids treating the diagonal as if it were just another copy of the
off-diagonal formula with no combinatorial distinction. In practice it is often the right
quantity when one wants local bunching, pair occupancy, or normal-ordered intuition.

But because the second moment is also useful, the code keeps `"plain"` available and
records the chosen convention in the results file so downstream readers do not have to
guess.

### 6. We store the physical local variance separately to remove ambiguity

The code now stores

```math
\mathrm{Var}(n_i) = \langle n_i^2 \rangle - \langle n_i \rangle^2
```

as `variance_a` and `variance_b`.

This is a design choice to reduce avoidable ambiguity. The diagonal of `connected_nn`
is not always the physical variance, because under the default `"factorial"` convention

```math
C_{ii} = \langle n_i(n_i-1) \rangle - \langle n_i \rangle^2,
```

which differs from `Var(n_i)` by `⟨n_i⟩`.

Rather than making every downstream script remember that correction, the pipeline stores
the physical local fluctuation explicitly and independently of the diagonal convention.

The philosophy here is:

- keep the convention-sensitive primitive matrix,
- but also store the convention-independent local fluctuation users usually mean.

### 7. Translational averages are convenience summaries, not replacements

The code also stores translational averages in displacement space:

```math
G(r) = \frac{1}{N_r}\sum_i \langle n_i n_{i+r} \rangle,
\quad
C(r) = \frac{1}{N_r}\sum_i \left[\langle n_i n_{i+r} \rangle - \langle n_i\rangle\langle n_{i+r}\rangle\right].
```

These are useful because they are compact, intuitive, and often what one plots first.
But they are not treated as the fundamental data product.

In homogeneous phases they may be close to sufficient. In impurity problems they are
only summaries of the more informative real-space matrix.

The same logic applies to the mixed-species average

```math
G^{(ab)}(r) = \frac{1}{N_r}\sum_i \langle N^{(a)}_i N^{(b)}_{i+r}\rangle,
```

which is useful, but still treated as a summary derived from the full `nn_ab` matrix.

### 8. Edge trimming is explicit because OBC can bias averages strongly

For open boundaries, translational averages can be distorted by edge effects. The code
therefore supports `edge_trim` and `edge_trims` and can store multiple profiles.

That is a philosophical choice as much as a numerical one. We do not want to pretend
that every site in an OBC chain is equally representative of the bulk. The trimming
parameters make the bulk-vs-edge tradeoff explicit and reproducible.

### 9. The fast backend is the default, but a reference backend remains available

The preferred implementation uses `ITensorMPS.correlation_matrix`, because it is the
natural efficient MPS algorithm for two-point correlators.

The legacy MPO-per-pair backend is retained as a reference path. That gives the project:

- a debugging baseline,
- a way to sanity-check small systems,
- a safeguard against future refactors that might change semantics unintentionally.

The philosophy is the same as elsewhere in the observables stack:

- use the efficient path in production,
- keep a slower path that is easy to reason about,
- and make agreement between the two a meaningful validation target.

## What We Want Users To Read Off These Data

The density-density outputs are meant to answer several different questions at once:

- the full matrix reveals spatial structure and impurity effects,
- the connected matrix isolates genuine fluctuation correlations,
- the mixed-species matrix reveals bath-impurity density locking or anticorrelation,
- the variance vectors show local compressibility-like behavior and number fluctuations,
- the translational averages provide compact radial summaries,
- the edge-trimmed profiles help separate bulk trends from boundary contamination.

No single one of these should be treated as universally "the" density-density observable.
They are a family of related views of the same underlying physics.

## Consistency Checks We Care About

The most important checks for density-density outputs are:

1. Symmetry:

```math
\mathrm{nnmat}_{ij} = \mathrm{nnmat}_{ji}
```

for these bosonic density operators.

For the mixed-species matrix, do not assume

```math
\mathrm{nn}^{(ab)}_{ij} = \mathrm{nn}^{(ab)}_{ji}.
```

That symmetry is generally absent because the row and column indices correspond to
different species.

2. Connected definition:

```math
\mathrm{connected\_nn}_{ij} = \mathrm{nnmat}_{ij} - n_i n_j.
```

3. Variance consistency:

```math
\mathrm{variance}_i = \langle n_i^2 \rangle - \langle n_i \rangle^2
```

independent of `same_site_convention`.

4. Translational averaging consistency:

the stored `transl_avg_*` data should agree with direct averaging over the stored matrix
using the recorded boundary and edge-trim conventions.

5. Mixed-species connected definition:

```math
\mathrm{connected\_nn}^{(ab)}_{ij} =
\mathrm{nn}^{(ab)}_{ij} - n^{(a)}_i n^{(b)}_j.
```

These checks are more important than the exact naming of datasets because they validate
that the stored observables still carry the intended physical meaning.

## Downstream View

The density-density group is designed to support multiple downstream analyses without
re-running DMRG:

- extracting local fluctuations,
- comparing impurity-centered and bulk correlations,
- building structure factors,
- testing edge sensitivity,
- and performing campaign-level summaries that need both full matrices and compact scalar
  indicators.

That is the guiding philosophy:

- preserve the primitive spatial information,
- make convention choices explicit,
- store connected and local-fluctuation views when they are genuinely useful,
- and let derived interpretation happen in analysis, not by throwing information away at
  measurement time.
