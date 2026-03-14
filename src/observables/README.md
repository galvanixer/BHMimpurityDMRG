# Observables Formula Reference

This note defines the formulas implemented in `/src/observables` for density-based measurements.

## Notation

- State: $|\psi\rangle$
- Chain length: $L$
- Site indices: $i, j, k \in \{1,\dots,L\}$
- Species-resolved density operator at site $i$: $n_i$ (implemented via `opname="Na"` or `opname="Nb"`)

Expectation value:

```math
\langle O \rangle = \langle \psi | O | \psi \rangle
```

Unless noted otherwise, the formulas below are for one chosen `opname` at a time.

## 1. One-point densities

Implemented by `expect_n`, `onsite_expect`, and `measure_densities`.

```math
\langle n_i \rangle
```

`measure_densities` returns:

```math
n_a[i] = \langle n_i^{(a)} \rangle,\quad
n_b[i] = \langle n_i^{(b)} \rangle
```

## 2. Two-point density correlators

Implemented by `density_density_matrix` and `cross_density_density_matrix`.

For the design rationale, diagonal-convention policy, and storage philosophy behind this
observable family, see `docs/density_density.md`.

For $i \neq j$:

```math
\mathrm{nnmat}_{ij} = \langle n_i n_j \rangle
```

For $i=j$, the code supports two conventions:

`same_site_convention="plain"`:

```math
\mathrm{nnmat}_{ii} = \langle n_i^2 \rangle
```

`same_site_convention="factorial"`:

```math
\mathrm{nnmat}_{ii} = \langle n_i(n_i-1) \rangle = \langle n_i^2 \rangle - \langle n_i \rangle
```

One-point vector:

```math
\mathrm{nvec}_i = \langle n_i \rangle
```

### Cross-species matrix

Implemented by `cross_density_density_matrix`.

For species `a` on the left and species `b` on the right:

```math
\mathrm{nnmat}^{(ab)}_{ij} = \langle n^{(a)}_i n^{(b)}_j \rangle
```

including on the diagonal:

```math
\mathrm{nnmat}^{(ab)}_{ii} = \langle n^{(a)}_i n^{(b)}_i \rangle.
```

Unlike the same-species case, there is no separate `"plain"` versus `"factorial"`
convention for the cross-species diagonal.

### Connected two-point matrix

Implemented by `connected_density_density_matrix`.

```math
C_{ij} = \mathrm{nnmat}_{ij} - \langle n_i \rangle \langle n_j \rangle
```

Note: with factorial diagonal,

```math
C_{ii} = \langle n_i(n_i-1)\rangle - \langle n_i \rangle^2
```

which is not the usual variance.

For the mixed-species matrix, the connected version is

```math
C^{(ab)}_{ij} = \langle n^{(a)}_i n^{(b)}_j \rangle - \langle n^{(a)}_i \rangle \langle n^{(b)}_j \rangle.
```

### Local number fluctuation

Implemented by `local_density_variance` and stored as `variance_a` / `variance_b`
inside the `/observables/density_density` group.

The stored quantity is always the physical on-site variance

```math
\mathrm{Var}(n_i) = \langle n_i^2 \rangle - \langle n_i \rangle^2
```

independent of whether `nnmat[i,i]` was stored with the `"plain"` or `"factorial"`
same-site convention.

## 3. Translational averages in displacement space

Implemented by `transl_avg_density_density` and `transl_avg_density_density_pair`.

For displacement $r$, define anchors $A_r$ as sites where $j=i+r$ is valid (open BC) or wrapped (periodic BC):

```math
G(r) = \frac{1}{N_r}\sum_{i\in A_r} \langle n_i n_{i+r}\rangle
```

```math
C(r) = \frac{1}{N_r}\sum_{i\in A_r}
\left[\langle n_i n_{i+r}\rangle - \langle n_i\rangle\langle n_{i+r}\rangle\right]
```

with:

```math
N_r = |A_r|
```

returned as `anchors`.

With OBC and `edge_trim=t`, anchors are restricted to the bulk window:

```math
i \in \{1+t,\dots,L-t\},\quad j=i+r \in \{1+t,\dots,L-t\}
```

so the effective anchor count is:

```math
N_r = L-r-2t
```

and only displacements with $N_r>0$ are included.

Boundary handling (`shifted_site`):

- `periodic=false`: include only if $1 \le i+r \le L$
- `periodic=true`: wrap with `mod1(i+r, L)`

If `fold_min_image=true` with periodic BC, `r` values are folded by `min_image(r, L)` before averaging.

For the mixed-species case, the same averaging pattern is used with

```math
G^{(ab)}(r) = \frac{1}{N_r}\sum_{i\in A_r} \langle n^{(a)}_i n^{(b)}_{i+r}\rangle
```

and

```math
C^{(ab)}(r) = \frac{1}{N_r}\sum_{i\in A_r}
\left[\langle n^{(a)}_i n^{(b)}_{i+r}\rangle - \langle n^{(a)}_i\rangle\langle n^{(b)}_{i+r}\rangle\right].
```

## 4. Static structure factor

Implemented by `structure_factor_from_nn`.

Wavevectors:

```math
k_m = \frac{2\pi m}{L}
```

with centered integer modes:

```math
m=
\begin{cases}
-L/2,\dots,L/2-1, & L\ \text{even}\\
-(L-1)/2,\dots,(L-1)/2, & L\ \text{odd}
\end{cases}
```

For even $L$, $m=+L/2$ is equivalent to $m=-L/2$ (same momentum modulo $2\pi$), so only one is kept.

Construct matrix `M` as:

1. Start from `nnmat`.
2. If `factorial_diagonal=true`, convert diagonal to plain moments:

```math
M_{ii} = \mathrm{nnmat}_{ii} + \langle n_i \rangle
```

3. If `connected=true`, subtract disconnected part:

```math
M_{ij} = M_{ij} - \langle n_i \rangle \langle n_j \rangle
```

Then compute:

```math
S(k_m) = \frac{1}{L}\sum_{i=1}^{L}\sum_{j=1}^{L} e^{ik_m(i-j)} M_{ij}
```

The implementation returns $\mathrm{Re}[S(k_m)]$.

## 5. Pair-distance distribution

Implemented by `pair_distance_distribution` and `cross_pair_distance_distribution`,
stored under `/observables/pair_distance_distribution`.
For design rationale and storage philosophy, see `docs/pair_distance_distribution.md`.

For one species, define the expected unordered pair-count profile:

```math
C(r)=
\begin{cases}
\sum_i \frac{1}{2}\langle n_i(n_i-1)\rangle, & r=0, \\
\sum_{i<j,\; d(i,j)=r}\langle n_i n_j\rangle, & r>0.
\end{cases}
```

The normalized distribution is

```math
P(r)=\frac{C(r)}{\sum_{r'} C(r')}.
```

For periodic boundaries, the code uses unsigned minimum-image distance:

```math
d(i,j)=\min(|i-j|,\,L-|i-j|),
```

so `r` runs from `0` to `floor(L/2)`; for open boundaries it runs from `0` to `L-1`.

For cross species (`a`,`b`), the pair-count profile is

```math
C^{(ab)}(r)=\sum_{i,j,\; d(i,j)=r}\langle n_i^{(a)} n_j^{(b)}\rangle,
```

with normalization

```math
P^{(ab)}(r)=\frac{C^{(ab)}(r)}{\sum_{r'} C^{(ab)}(r')}.
```

The pipeline also stores moments:

```math
\langle r\rangle = \sum_r r P(r),\quad
\mathrm{Var}(r)=\sum_r r^2 P(r)-\langle r\rangle^2,\quad
\sigma_r=\sqrt{\mathrm{Var}(r)}.
```

## 6. Three-point density correlators

### Plain (non-normal-ordered)

Implemented by `expect_nnn`, `connected_nnn`.

Raw moment:

```math
\langle n_i n_j n_k \rangle
```

Connected third cumulant:

```math
C^{(3)}_{ijk} =
\langle n_i n_j n_k \rangle
- \langle n_i n_j \rangle\langle n_k \rangle
- \langle n_i n_k \rangle\langle n_j \rangle
- \langle n_j n_k \rangle\langle n_i \rangle
+ 2\langle n_i \rangle\langle n_j \rangle\langle n_k \rangle
```

### Normal-ordered moments

Implemented by `expect_nn_no`, `expect_nnn_no`.

Two-point:

```math
\langle :n_i n_j: \rangle =
\begin{cases}
\langle n_i n_j \rangle, & i\neq j \\
\langle n_i(n_i-1)\rangle, & i=j
\end{cases}
```

Three-point:

All distinct:

```math
\langle :n_i n_j n_k: \rangle = \langle n_i n_j n_k \rangle
```

One repeated index (example $i=j\neq k$):

```math
\langle :n_i^2 n_k: \rangle = \langle n_i^2 n_k \rangle - \langle n_i n_k \rangle
```

All equal:

```math
\langle :n_i^3: \rangle = \langle n_i(n_i-1)(n_i-2)\rangle
= \langle n_i^3\rangle - 3\langle n_i^2\rangle + 2\langle n_i\rangle
```

Connected normal-ordered third cumulant (implemented by `connected_nnn_no`):

```math
C^{(3),\mathrm{no}}_{ijk} =
\langle :n_i n_j n_k: \rangle
- \langle :n_i n_j: \rangle \langle n_k \rangle
- \langle :n_i n_k: \rangle \langle n_j \rangle
- \langle :n_j n_k: \rangle \langle n_i \rangle
+ 2\langle n_i \rangle\langle n_j \rangle\langle n_k \rangle
```

Cached variants (`precompute_n`, `precompute_nn`, `*_cached`) use the same formulas with reused one- and two-point moments.

## 7. Translationally averaged three-point functions

Implemented by:

- `transl_avg_nnn`
- `transl_avg_connected_nnn`
- `transl_avg_nnn_no`
- `transl_avg_connected_nnn_no`
- `transl_avg_connected_nnn_no_cached`

For displacements $(r,s)$:

```math
G^{(3)}(r,s) = \frac{1}{N_{r,s}}\sum_{i\in A_{r,s}} \langle n_i n_{i+r} n_{i+s}\rangle
```

```math
C^{(3)}(r,s) = \frac{1}{N_{r,s}}\sum_{i\in A_{r,s}} C^{(3)}_{i,i+r,i+s}
```

and similarly for normal-ordered quantities.

$N_{r,s}$ is the number of valid anchors returned by these functions.

## 8. Single-particle density matrix

Implemented by `single_particle_density_matrix`.

For the design rationale, storage policy, and interpretation philosophy behind this
observable, see `docs/single_particle_density_matrix.md`.

For each species separately, the code evaluates:

```math
G^{(1)}_{ij} = \langle a_i^\dagger a_j \rangle
```

or

```math
G^{(1)}_{ij} = \langle b_i^\dagger b_j \rangle
```

depending on the selected species.

On the diagonal, this reduces to the site density:

```math
G^{(1)}_{ii} = \langle n_i \rangle
```

For a fixed species, the matrix is Hermitian:

```math
G^{(1)}_{ij} = \left(G^{(1)}_{ji}\right)^*
```
