# Observables Formula Reference (GitHub Math)

This note defines the formulas implemented in `/src/observables` for density-based measurements.
For a plain-text-equation version compatible with non-math renderers, see `README.md`.

## Notation

- State: `|\psi\rangle`
- Chain length: `L`
- Site indices: `i, j, k \in \{1,\dots,L\}`
- Species-resolved density operator at site `i`: `n_i` (implemented via `opname="Na"` or `opname="Nb"`)

Expectation value:

```math
\langle O \rangle = \langle \psi | O | \psi \rangle
```

All formulas below are for one chosen `opname` at a time.

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

Implemented by `density_density_matrix`.

For `i \neq j`:

```math
\mathrm{nnmat}_{ij} = \langle n_i n_j \rangle
```

For `i=j`, the code supports two conventions:

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

## 3. Translational averages in displacement space

Implemented by `transl_avg_density_density`.

For displacement `r`, define anchors `A_r` as sites where `j=i+r` is valid (open BC) or wrapped (periodic BC):

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

Boundary handling (`shifted_site`):

- `periodic=false`: include only if `1 \le i+r \le L`
- `periodic=true`: wrap with `mod1(i+r, L)`

If `fold_min_image=true` with periodic BC, `r` values are folded by `min_image(r, L)` before averaging.

## 4. Static structure factor

Implemented by `structure_factor_from_nn`.

Wavevectors:

```math
k_m = \frac{2\pi m}{L},\quad m=0,\dots,L-1
```

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

The implementation returns `\mathrm{Re}[S(k_m)]`.

## 5. Three-point density correlators

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

One repeated index (example `i=j\neq k`):

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

## 6. Translationally averaged three-point functions

Implemented by:

- `transl_avg_nnn`
- `transl_avg_connected_nnn`
- `transl_avg_nnn_no`
- `transl_avg_connected_nnn_no`
- `transl_avg_connected_nnn_no_cached`

For displacements `(r,s)`:

```math
G^{(3)}(r,s) = \frac{1}{N_{r,s}}\sum_{i\in A_{r,s}} \langle n_i n_{i+r} n_{i+s}\rangle
```

```math
C^{(3)}(r,s) = \frac{1}{N_{r,s}}\sum_{i\in A_{r,s}} C^{(3)}_{i,i+r,i+s}
```

and similarly for normal-ordered quantities.

`N_{r,s}` is the number of valid anchors returned by these functions.
