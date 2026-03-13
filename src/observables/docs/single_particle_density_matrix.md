# Single-Particle Density Matrix

This note explains how the code evaluates the single-particle density matrix (SPDM),
what exactly is stored, and the philosophy behind those choices.

The implementation lives in `src/observables/single_particle_density_matrix.jl`.
The canonical HDF5 output lives under:

- `/observables/single_particle_density_matrix/G1_a`
- `/observables/single_particle_density_matrix/G1_b`

for the bath (`a`) and impurity (`b`) species respectively.

## What We Measure

For each species separately, we define the one-body density matrix

```math
G^{(1)}_{ij} = \langle a_i^\dagger a_j \rangle
```

or

```math
G^{(1)}_{ij} = \langle b_i^\dagger b_j \rangle.
```

This is the standard object for diagnosing first-order coherence, natural orbitals,
condensation, and momentum-space structure.

On the diagonal,

```math
G^{(1)}_{ii} = \langle n_i \rangle,
```

so the SPDM is consistent with, and contains on its diagonal, the measured density profile.

## Philosophy

### 1. The full matrix is the primary observable

We store the full site-resolved matrix, not only a translational average and not only
its Fourier transform.

That is deliberate.

This code is built for inhomogeneous impurity problems. In those systems the interesting
physics is often local and spatially structured:

- coherence can be suppressed near the impurity but recover in the bulk,
- open boundaries can distort correlations near the edges,
- translational symmetry can be weakly broken or absent by design.

If we only stored an averaged `G^{(1)}(r)` or only stored a momentum distribution `n(k)`,
we would throw away information that is often the whole point of the calculation.

The full matrix is the most primitive object. Anything simpler can be derived from it later.

### 2. Species are kept separate

We store `G1_a` and `G1_b` separately rather than combining them.

That matches the physical model and keeps interpretation clean:

- `G1_a` probes bath coherence,
- `G1_b` probes impurity coherence,
- comparing the two is often more informative than any summed quantity.

We do not currently store mixed-species one-body correlators such as
`⟨a_i^\dagger b_j⟩`. Those are different observables with a different interpretation,
and they should be added explicitly if needed rather than being bundled into the SPDM
silently.

### 3. We store the raw SPDM, not a "connected" version

For density-density observables it is often useful to subtract disconnected pieces.
For the one-body density matrix that is not the right default object.

The standard SPDM in many-body physics is the raw correlator `⟨a_i^\dagger a_j⟩`.
That object already has a direct physical meaning:

- its eigenvectors are the natural orbitals,
- its eigenvalues are the natural occupations,
- its Fourier transform gives the momentum distribution when appropriate.

Also, in number-conserving calculations the anomalous one-point expectations
`⟨a_i⟩` and `⟨b_i⟩` are generally zero, so a connected subtraction would not add useful
information in the typical workflow here.

### 4. We keep the matrix complex

The SPDM is Hermitian for a fixed species,

```math
G^{(1)}_{ij} = \left(G^{(1)}_{ji}\right)^*,
```

but it is not guaranteed to be purely real.

That matters. Complex phases can carry physical information, especially in settings with:

- twisted boundary conditions,
- currents,
- gauge fields,
- complex variational states.

So the code stores `G1_a` and `G1_b` as complex matrices. We do not discard the imaginary
part just because many benchmark states happen to be real.

### 5. We exploit Hermiticity in the evaluator

The default backend assumes Hermiticity through `ishermitian=true`.

This is a performance choice, not a change of definition. For a fixed species,
`⟨a_i^\dagger a_j⟩` and `⟨b_i^\dagger b_j⟩` are Hermitian matrices, so filling one triangle
and reconstructing the other is legitimate and cheaper than evaluating every entry
independently.

If that assumption ever needs to be relaxed for debugging or comparison, the legacy path
can evaluate all entries explicitly.

### 6. We prefer the fast MPS correlator, but keep a reference backend

The main implementation uses `ITensorMPS.correlation_matrix` with operator pairs:

- `("Adag", "A")` for species `a`
- `("Bdag", "B")` for species `b`

This is the right production path because it is the native efficient MPS algorithm for
two-point correlators.

We also keep a legacy MPO-per-pair implementation as a reference backend. That is useful
for:

- debugging,
- validating small systems,
- checking future refactors,
- comparing results if there is ever suspicion about the fast path.

The philosophy is simple: fast by default, but never with no reference point.

### 7. We do not over-interpret translational symmetry

The SPDM can certainly be compressed to displacement-space averages,

```math
G^{(1)}(r) = \frac{1}{N_r}\sum_i \langle a_i^\dagger a_{i+r} \rangle,
```

and for homogeneous systems that can be very natural.

But it is not the default here, because this project is not centered on perfectly
homogeneous states. We want the observable layer to preserve inhomogeneity first and
average later only when the user has decided that such an average is physically justified.

## Practical Interpretation

The SPDM should be read as the coherence map of the state.

- Large off-diagonal magnitude means long-range coherence.
- Rapid decay away from the diagonal means short coherence length.
- The diagonal reproduces the density profile.
- The leading eigenvector is the dominant natural orbital.
- The leading eigenvalue measures how strongly a single orbital is occupied.

In impurity problems, a useful question is not only "how large is the condensate-like
piece?" but also "where does it live?" The full matrix answers that more directly than
an already-averaged scalar diagnostic.

## Consistency Checks We Care About

Whenever the SPDM is computed, the following checks are the most meaningful ones:

1. Diagonal consistency:

```math
\mathrm{diag}(G1_a) = n_a,\quad \mathrm{diag}(G1_b) = n_b
```

2. Hermiticity:

```math
G1 = G1^\dagger
```

up to numerical tolerance.

3. Cross-backend agreement on small systems:

the fast backend and the legacy backend should agree within numerical error.

These checks matter more than formatting details because they validate the physical
meaning of the stored matrix.

## Why This Design Is Useful Downstream

The stored SPDM is intentionally a foundation for later analysis. From `G1_a` or `G1_b`
you can derive:

- momentum distributions,
- natural orbitals,
- condensate fractions,
- coherence lengths,
- bulk-restricted or edge-restricted averages,
- impurity-centered coherence diagnostics.

Those are downstream analysis choices. The measurement layer should preserve the primary
object cleanly enough that all of those remain available.

That is the main philosophy here:

- store the primitive observable,
- preserve spatial structure,
- keep species separate,
- keep phases,
- use the fast algorithm by default,
- and make derived summaries a postprocessing decision rather than a measurement-time loss
  of information.
