# Observables Formula Reference

This note defines the formulas implemented in `/src/observables` for density-based measurements.
For a GitHub math-rendered version, see `README_math.md`.

## Notation

- State: `|psi>`
- Chain length: `L`
- Site indices: `i, j, k in {1, ..., L}`
- Species-resolved density operator at site `i`: `n_i` (implemented via `opname="Na"` or `opname="Nb"`)
- Expectation value: `<O> = <psi| O |psi>`

All formulas below are for one chosen `opname` at a time.

## 1. One-point densities

Implemented by `expect_n`, `onsite_expect`, and `measure_densities`.

```text
<n_i>
```

`measure_densities` returns:

```text
na[i] = <n_i^(a)>
nb[i] = <n_i^(b)>
```

## 2. Two-point density correlators

Implemented by `density_density_matrix`.

For `i != j`:

```text
nnmat[i,j] = <n_i n_j>
```

For `i == j`, the code supports two conventions:

- `same_site_convention="plain"`:

```text
nnmat[i,i] = <n_i^2>
```

- `same_site_convention="factorial"`:

```text
nnmat[i,i] = <n_i (n_i - 1)> = <n_i^2> - <n_i>
```

One-point vector:

```text
nvec[i] = <n_i>
```

### Connected two-point matrix

Implemented by `connected_density_density_matrix`:

```text
C[i,j] = nnmat[i,j] - <n_i><n_j>
```

Note: with factorial diagonal,

```text
C[i,i] = <n_i (n_i - 1)> - <n_i>^2
```

which is not the usual variance.

## 3. Translational averages in displacement space

Implemented by `transl_avg_density_density`.

For displacement `r`, define anchors `A_r` as sites where `j=i+r` is valid (open BC) or wrapped (periodic BC):

```text
G(r) = (1 / N_r) * sum_{i in A_r} <n_i n_{i+r}>
C(r) = (1 / N_r) * sum_{i in A_r} [<n_i n_{i+r}> - <n_i><n_{i+r}>]
```

with:

```text
N_r = |A_r|
```

returned as `anchors`.

Boundary handling (`shifted_site`):

- `periodic=false`: include only if `1 <= i+r <= L`
- `periodic=true`: wrap with `mod1(i+r, L)`

If `fold_min_image=true` with periodic BC, `r` values are folded by `min_image(r, L)` before averaging.

## 4. Static structure factor

Implemented by `structure_factor_from_nn`.

Wavevectors:

```text
k_m = 2*pi*m / L,  m = 0, ..., L-1
```

Construct matrix `M` as:

1. Start from `nnmat`.
2. If `factorial_diagonal=true`, convert diagonal to plain moments:

```text
M[i,i] = nnmat[i,i] + <n_i>
```

3. If `connected=true`, subtract disconnected part:

```text
M[i,j] = M[i,j] - <n_i><n_j>
```

Then compute:

```text
S(k_m) = (1/L) * sum_{i=1}^L sum_{j=1}^L exp(i*k_m*(i-j)) * M[i,j]
```

The implementation returns `real(S(k_m))`.

## 5. Three-point density correlators

### Plain (non-normal-ordered)

Implemented by `expect_nnn`, `connected_nnn`.

Raw moment:

```text
<n_i n_j n_k>
```

Connected third cumulant:

```text
C3[i,j,k] =
  <n_i n_j n_k>
  - <n_i n_j><n_k>
  - <n_i n_k><n_j>
  - <n_j n_k><n_i>
  + 2<n_i><n_j><n_k>
```

### Normal-ordered moments

Implemented by `expect_nn_no`, `expect_nnn_no`.

Two-point:

```text
<:n_i n_j:> = <n_i n_j>                 if i != j
<:n_i n_i:> = <n_i (n_i - 1)>           if i == j
```

Three-point:

- all distinct:

```text
<:n_i n_j n_k:> = <n_i n_j n_k>
```

- one repeated index (example `i=j!=k`):

```text
<:n_i^2 n_k:> = <n_i^2 n_k> - <n_i n_k>
```

- all equal:

```text
<:n_i^3:> = <n_i (n_i - 1) (n_i - 2)>
         = <n_i^3> - 3<n_i^2> + 2<n_i>
```

Connected normal-ordered third cumulant (implemented by `connected_nnn_no`):

```text
C3_no[i,j,k] =
  <:n_i n_j n_k:>
  - <:n_i n_j:><n_k>
  - <:n_i n_k:><n_j>
  - <:n_j n_k:><n_i>
  + 2<n_i><n_j><n_k>
```

Cached variants (`precompute_n`, `precompute_nn`, `*_cached`) use the same formulas with reused one- and two-point moments.

## 6. Translationally averaged three-point functions

Implemented by:

- `transl_avg_nnn`
- `transl_avg_connected_nnn`
- `transl_avg_nnn_no`
- `transl_avg_connected_nnn_no`
- `transl_avg_connected_nnn_no_cached`

For displacements `(r, s)`:

```text
G3(r,s) = (1 / N_{r,s}) * sum_{i in A_{r,s}} <n_i n_{i+r} n_{i+s}>
C3(r,s) = (1 / N_{r,s}) * sum_{i in A_{r,s}} C3[i, i+r, i+s]
```

and similarly for normal-ordered quantities.

`N_{r,s}` is the number of valid anchors returned by these functions.
