# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

# ----------------------------
# Pair-distance distribution P(r)
# ----------------------------

@inline function _pair_distance_rvals(L::Int; periodic::Bool)
    rmax = periodic ? (L ÷ 2) : (L - 1)
    return collect(0:rmax)
end

"""
    pair_distance_distribution(
        nvec,
        nnmat;
        periodic::Bool=false,
        same_site_convention::String="factorial"
    )

Compute the same-species pair-distance distribution from one-point and two-point
density moments.

Returns `(r, C, P, pair_count, mean_r, var_r, std_r)` where:
- `r`: supported distances
- `C[r]`: expected number of unordered pairs at distance `r`
- `P[r]`: normalized pair-distance distribution (`C / pair_count`)
- `pair_count`: expected total number of unordered pairs
- `mean_r`, `var_r`, `std_r`: first two moments of `P(r)`

For periodic systems, distances are folded by the unsigned minimum-image rule.
"""
function pair_distance_distribution(
    nvec::AbstractVector,
    nnmat::AbstractMatrix;
    periodic::Bool=false,
    same_site_convention::String="factorial"
)
    L = length(nvec)
    size(nnmat, 1) == L && size(nnmat, 2) == L ||
        throw(ArgumentError("nnmat must be LxL with L=length(nvec)"))
    conv = _validate_same_site_convention(same_site_convention)

    rvals = _pair_distance_rvals(L; periodic=periodic)
    C = zeros(Float64, length(rvals))

    # On-site unordered boson pairs: 1/2 * <n_i (n_i - 1)>
    @inbounds for i in 1:L
        nfac = conv == "factorial" ? nnmat[i, i] : (nnmat[i, i] - nvec[i])
        C[1] += 0.5 * nfac
    end

    # Off-site unordered site pairs: sum over i < j only.
    @inbounds for i in 1:(L - 1), j in (i + 1):L
        d = j - i
        if periodic
            d = min(d, L - d)
        end
        C[d + 1] += nnmat[i, j]
    end

    pair_count = sum(C)
    if pair_count <= 0.0
        P = zeros(Float64, length(rvals))
        return rvals, C, P, pair_count, NaN, NaN, NaN
    end

    P = C ./ pair_count
    mean_r = sum(rvals .* P)
    var_r = sum((rvals .^ 2) .* P) - mean_r^2
    std_r = sqrt(max(var_r, 0.0))
    return rvals, C, P, pair_count, mean_r, var_r, std_r
end

"""
    cross_pair_distance_distribution(
        nvec_left,
        nvec_right,
        nnmat;
        periodic::Bool=false
    )

Compute the cross-species pair-distance distribution for distinguishable species
from `nnmat[i,j] = <n_i^(left) n_j^(right)>`.

Returns `(r, C, P, pair_count, mean_r, var_r, std_r)` where `pair_count` is
the expected number of cross-species pairs.
"""
function cross_pair_distance_distribution(
    nvec_left::AbstractVector,
    nvec_right::AbstractVector,
    nnmat::AbstractMatrix;
    periodic::Bool=false
)
    L = length(nvec_left)
    length(nvec_right) == L ||
        throw(ArgumentError("nvec_right must have length L=length(nvec_left)"))
    size(nnmat, 1) == L && size(nnmat, 2) == L ||
        throw(ArgumentError("nnmat must be LxL with L=length(nvec_left)"))

    rvals = _pair_distance_rvals(L; periodic=periodic)
    C = zeros(Float64, length(rvals))

    @inbounds for i in 1:L, j in 1:L
        d = abs(j - i)
        if periodic
            d = min(d, L - d)
        end
        C[d + 1] += nnmat[i, j]
    end

    pair_count = sum(C)
    if pair_count <= 0.0
        P = zeros(Float64, length(rvals))
        return rvals, C, P, pair_count, NaN, NaN, NaN
    end

    P = C ./ pair_count
    mean_r = sum(rvals .* P)
    var_r = sum((rvals .^ 2) .* P) - mean_r^2
    std_r = sqrt(max(var_r, 0.0))
    return rvals, C, P, pair_count, mean_r, var_r, std_r
end
