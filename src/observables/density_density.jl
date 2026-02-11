# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

# ----------------------------
# Two-point density correlators
# ----------------------------

"""
    density_density_matrix(psi, sites, opname::String; same_site_convention::String="factorial")

Return `(nvec, nnmat)` where:
- `nvec[i] = ⟨n_i⟩`
- `nnmat[i,j] = ⟨n_i n_j⟩` for `i != j`
- `nnmat[i,i] = ⟨n_i(n_i-1)⟩` if `same_site_convention=="factorial"`
- `nnmat[i,i] = ⟨n_i^2⟩` if `same_site_convention=="plain"`
"""
function density_density_matrix(
    psi::MPS,
    sites,
    opname::String;
    same_site_convention::String="factorial"
)
    L = length(sites)
    nvec = zeros(Float64, L)
    nnmat = zeros(Float64, L, L)
    conv = lowercase(same_site_convention)
    conv in ("factorial", "plain") ||
        throw(ArgumentError("same_site_convention must be \"factorial\" or \"plain\""))

    for i in 1:L
        nvec[i] = expect_n(psi, sites, opname, i)
        n2 = expect_nn(psi, sites, opname, i, i) # <n_i^2>
        nnmat[i, i] = conv == "factorial" ? (n2 - nvec[i]) : n2
        for j in (i + 1):L
            v = expect_nn(psi, sites, opname, i, j)
            nnmat[i, j] = v
            nnmat[j, i] = v
        end
    end
    return nvec, nnmat
end

"""
    connected_density_density_matrix(nvec, nnmat)

Return connected two-point matrix:
`C[i,j] = nnmat[i,j] - ⟨n_i⟩⟨n_j⟩`.
"""
function connected_density_density_matrix(nvec::AbstractVector, nnmat::AbstractMatrix)
    return nnmat .- nvec * transpose(nvec)
end

"""
    transl_avg_density_density(nvec, nnmat; periodic::Bool=false, max_r=nothing, fold_min_image::Bool=false)

Return translationally averaged two-point correlators over anchor site `i`:
- `rvals`: displacement values `0:rmax`
- `g`: raw average `G(r) = (1/N_r) Σ_i ⟨n_i n_{i+r}⟩`
- `c`: connected average `C(r) = (1/N_r) Σ_i [⟨n_i n_{i+r}⟩ - ⟨n_i⟩⟨n_{i+r}⟩]`
- `anchors`: number of valid anchors `N_r` for each `r`

At `r=0`, the diagonal convention follows the input `nnmat`.
If `fold_min_image=true` and `periodic=true`, displacements are folded with
`min_image(r, L)` and unique folded displacements are used.
"""
function transl_avg_density_density(
    nvec::AbstractVector,
    nnmat::AbstractMatrix;
    periodic::Bool=false,
    max_r::Union{Int,Nothing}=nothing,
    fold_min_image::Bool=false
)
    L = length(nvec)
    size(nnmat, 1) == L && size(nnmat, 2) == L ||
        throw(ArgumentError("nnmat must be LxL with L=length(nvec)"))

    rmax = max_r === nothing ? (L - 1) : min(max_r, L - 1)
    if fold_min_image && periodic
        rvals = sort(unique([min_image(r, L) for r in 0:rmax]))
    else
        rvals = collect(0:rmax)
    end
    g = zeros(Float64, length(rvals))
    c = zeros(Float64, length(rvals))
    anchors = zeros(Int, length(rvals))

    for (idx, r) in enumerate(rvals)
        acc_g = 0.0
        acc_c = 0.0
        n = 0
        for i in 1:L
            j = shifted_site(i, r, L; periodic=periodic)
            j === nothing && continue
            v = nnmat[i, j]
            acc_g += v
            acc_c += v - nvec[i] * nvec[j]
            n += 1
        end
        n == 0 && throw(ArgumentError("no valid anchors for r=$r with L=$L (periodic=$periodic)"))
        g[idx] = acc_g / n
        c[idx] = acc_c / n
        anchors[idx] = n
    end

    return rvals, g, c, anchors
end
