# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

# ----------------------------
# Two-point density correlators
# ----------------------------

@inline function _validate_same_site_convention(same_site_convention::String)
    conv = lowercase(same_site_convention)
    conv in ("factorial", "plain") ||
        throw(ArgumentError("same_site_convention must be \"factorial\" or \"plain\""))
    return conv
end

"""
    density_density_matrix_legacy(psi, sites, opname::String; same_site_convention::String="factorial")

Legacy implementation based on repeated `expect_nn` evaluations.
Returns `(nvec, nnmat)` with diagonal convention controlled by
`same_site_convention` (`"factorial"` or `"plain"`).
"""
function density_density_matrix_legacy(
    psi::MPS,
    sites,
    opname::String;
    same_site_convention::String="factorial"
)
    L = length(sites)
    nvec = zeros(Float64, L)
    nnmat = zeros(Float64, L, L)
    conv = _validate_same_site_convention(same_site_convention)

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
    density_density_matrix_correlation(psi, sites, opname::String;
                                       same_site_convention::String="factorial",
                                       ishermitian::Bool=true)

Fast implementation based on `ITensorMPS.correlation_matrix`.
"""
function density_density_matrix_correlation(
    psi::MPS,
    sites,
    opname::String;
    same_site_convention::String="factorial",
    ishermitian::Bool=true
)
    L = length(sites)
    length(psi) == L ||
        throw(ArgumentError("sites length ($L) must match MPS length ($(length(psi)))"))
    conv = _validate_same_site_convention(same_site_convention)

    site_range = 1:L
    n_raw = expect(psi, opname; sites=site_range)
    nvec = Float64.(real.(collect(n_raw)))
    nnmat = Float64.(real.(correlation_matrix(
        psi,
        opname,
        opname;
        sites=site_range,
        ishermitian=ishermitian
    )))

    if conv == "factorial"
        @inbounds for i in 1:L
            nnmat[i, i] -= nvec[i]
        end
    end

    return nvec, nnmat
end

"""
    density_density_matrix(psi, sites, opname::String;
                           same_site_convention::String="factorial",
                           backend::Union{Symbol,AbstractString}=:correlation_matrix,
                           ishermitian::Bool=true)

Compute `(nvec, nnmat)` using the selected backend:
- `backend=:correlation_matrix` (default): fast MPS correlation routine
- `backend=:legacy`: previous MPO-per-pair implementation
"""
function density_density_matrix(
    psi::MPS,
    sites,
    opname::String;
    same_site_convention::String="factorial",
    backend::Union{Symbol,AbstractString}=:correlation_matrix,
    ishermitian::Bool=true
)
    b = lowercase(String(backend))
    if b in ("correlation_matrix", "correlation", "fast")
        return density_density_matrix_correlation(
            psi,
            sites,
            opname;
            same_site_convention=same_site_convention,
            ishermitian=ishermitian
        )
    elseif b in ("legacy", "mpo", "opsum")
        return density_density_matrix_legacy(
            psi,
            sites,
            opname;
            same_site_convention=same_site_convention
        )
    end
    throw(ArgumentError(
        "backend must be one of :correlation_matrix/:legacy (also accepts correlation, fast, mpo, opsum); got \"$backend\""
    ))
end

"""
    cross_density_density_matrix_legacy(psi, sites, opname_left::String, opname_right::String)

Reference implementation of the mixed-species density-density matrix
`nnmat[i,j] = ⟨n_i^{(left)} n_j^{(right)}⟩`.
"""
function cross_density_density_matrix_legacy(
    psi::MPS,
    sites,
    opname_left::String,
    opname_right::String
)
    L = length(sites)
    length(psi) == L ||
        throw(ArgumentError("sites length ($L) must match MPS length ($(length(psi)))"))

    nnmat = zeros(Float64, L, L)
    for i in 1:L, j in 1:L
        nnmat[i, j] = real(expect_product(psi, sites, [(opname_left, i), (opname_right, j)]))
    end
    return nnmat
end

"""
    cross_density_density_matrix_correlation(psi, sites, opname_left::String, opname_right::String;
                                             ishermitian::Bool=false)

Fast implementation of the mixed-species density-density matrix based on
`ITensorMPS.correlation_matrix`.
"""
function cross_density_density_matrix_correlation(
    psi::MPS,
    sites,
    opname_left::String,
    opname_right::String;
    ishermitian::Bool=false
)
    L = length(sites)
    length(psi) == L ||
        throw(ArgumentError("sites length ($L) must match MPS length ($(length(psi)))"))

    site_range = 1:L
    return Float64.(real.(correlation_matrix(
        psi,
        opname_left,
        opname_right;
        sites=site_range,
        ishermitian=ishermitian
    )))
end

"""
    cross_density_density_matrix(psi, sites, opname_left::String, opname_right::String;
                                 backend=:correlation_matrix,
                                 ishermitian::Bool=false)

Compute the mixed-species density-density matrix `⟨n_i^{(left)} n_j^{(right)}⟩`.
"""
function cross_density_density_matrix(
    psi::MPS,
    sites,
    opname_left::String,
    opname_right::String;
    backend::Union{Symbol,AbstractString}=:correlation_matrix,
    ishermitian::Bool=false
)
    b = lowercase(String(backend))
    if b in ("correlation_matrix", "correlation", "fast")
        return cross_density_density_matrix_correlation(
            psi,
            sites,
            opname_left,
            opname_right;
            ishermitian=ishermitian
        )
    elseif b in ("legacy", "mpo", "opsum")
        return cross_density_density_matrix_legacy(
            psi,
            sites,
            opname_left,
            opname_right
        )
    end
    throw(ArgumentError(
        "backend must be one of :correlation_matrix/:legacy (also accepts correlation, fast, mpo, opsum); got \"$backend\""
    ))
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
    connected_cross_density_density_matrix(nvec_left, nvec_right, nnmat)

Return connected mixed-species two-point matrix:
`C[i,j] = nnmat[i,j] - ⟨n_i^{(left)}⟩⟨n_j^{(right)}⟩`.
"""
function connected_cross_density_density_matrix(
    nvec_left::AbstractVector,
    nvec_right::AbstractVector,
    nnmat::AbstractMatrix
)
    return nnmat .- nvec_left * transpose(nvec_right)
end

"""
    local_density_variance(nvec, nnmat; same_site_convention::String="factorial")

Return the physical local number fluctuation vector
`Var(n_i) = ⟨n_i^2⟩ - ⟨n_i⟩^2`, independent of the stored diagonal convention.
"""
function local_density_variance(
    nvec::AbstractVector,
    nnmat::AbstractMatrix;
    same_site_convention::String="factorial"
)
    L = length(nvec)
    size(nnmat, 1) == L && size(nnmat, 2) == L ||
        throw(ArgumentError("nnmat must be LxL with L=length(nvec)"))
    conv = _validate_same_site_convention(same_site_convention)

    variance = zeros(Float64, L)
    @inbounds for i in 1:L
        n2 = conv == "factorial" ? (nnmat[i, i] + nvec[i]) : nnmat[i, i]
        variance[i] = n2 - nvec[i]^2
    end
    return variance
end

"""
    transl_avg_density_density_pair(
        nvec_left,
        nvec_right,
        nnmat;
        periodic::Bool=false,
        max_r=nothing,
        fold_min_image::Bool=false,
        edge_trim::Int=0
    )

Return translationally averaged two-point correlators over anchor site `i` for a
possibly mixed pair of density operators.
"""
function transl_avg_density_density_pair(
    nvec_left::AbstractVector,
    nvec_right::AbstractVector,
    nnmat::AbstractMatrix;
    periodic::Bool=false,
    max_r::Union{Int,Nothing}=nothing,
    fold_min_image::Bool=false,
    edge_trim::Int=0
)
    L = length(nvec_left)
    length(nvec_right) == L ||
        throw(ArgumentError("nvec_right must have length L=length(nvec_left)"))
    size(nnmat, 1) == L && size(nnmat, 2) == L ||
        throw(ArgumentError("nnmat must be LxL with L=length(nvec_left)"))
    edge_trim >= 0 || throw(ArgumentError("edge_trim must be >= 0"))
    periodic && edge_trim > 0 &&
        throw(ArgumentError("edge_trim > 0 is only supported for periodic=false"))
    !periodic && (2 * edge_trim >= L) &&
        throw(ArgumentError("edge_trim too large for L=$L; require 2*edge_trim < L"))

    rmax_limit = periodic ? (L - 1) : (L - 1 - 2 * edge_trim)
    rmax = max_r === nothing ? rmax_limit : min(max_r, rmax_limit)
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
        i_min = periodic ? 1 : (1 + edge_trim)
        i_max = periodic ? L : (L - edge_trim)
        for i in i_min:i_max
            j = shifted_site(i, r, L; periodic=periodic)
            j === nothing && continue
            if !periodic && (j < i_min || j > i_max)
                continue
            end
            v = nnmat[i, j]
            acc_g += v
            acc_c += v - nvec_left[i] * nvec_right[j]
            n += 1
        end
        n == 0 && throw(ArgumentError("no valid anchors for r=$r with L=$L (periodic=$periodic)"))
        g[idx] = acc_g / n
        c[idx] = acc_c / n
        anchors[idx] = n
    end

    return rvals, g, c, anchors
end

"""
    transl_avg_density_density(
        nvec,
        nnmat;
        periodic::Bool=false,
        max_r=nothing,
        fold_min_image::Bool=false,
        edge_trim::Int=0
    )

Return translationally averaged two-point correlators over anchor site `i`:
- `rvals`: displacement values `0:rmax`
- `g`: raw average `G(r) = (1/N_r) Σ_i ⟨n_i n_{i+r}⟩`
- `c`: connected average `C(r) = (1/N_r) Σ_i [⟨n_i n_{i+r}⟩ - ⟨n_i⟩⟨n_{i+r}⟩]`
- `anchors`: number of valid anchors `N_r` for each `r`

`edge_trim` excludes `edge_trim` sites on each edge from anchor selection
for open boundaries. For OBC this corresponds to:
`i ∈ [1 + edge_trim, L - edge_trim]`, with pair partner `j=i+r` also required
to lie in the same interval. Thus for OBC the usable range is
`r <= L - 1 - 2*edge_trim`.

At `r=0`, the diagonal convention follows the input `nnmat`.
If `fold_min_image=true` and `periodic=true`, displacements are folded with
`min_image(r, L)` and unique folded displacements are used.
"""
function transl_avg_density_density(
    nvec::AbstractVector,
    nnmat::AbstractMatrix;
    periodic::Bool=false,
    max_r::Union{Int,Nothing}=nothing,
    fold_min_image::Bool=false,
    edge_trim::Int=0
)
    return transl_avg_density_density_pair(
        nvec,
        nvec,
        nnmat;
        periodic=periodic,
        max_r=max_r,
        fold_min_image=fold_min_image,
        edge_trim=edge_trim
    )
end
