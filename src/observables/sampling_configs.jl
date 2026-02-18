# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

# ----------------------------
# Sampling dominant basis configurations from an MPS
# ----------------------------

using Random

"""
    top_sampled_configurations(psi::MPS; nsamples::Int=100_000, top_k::Int=20, seed::Union{Nothing,Int}=nothing)

Sample product-basis configurations `x` from the Born distribution `|psi(x)|^2` and
return the most frequently observed configurations.

The input `psi` is not mutated. Sampling is performed on an internal copy that is
orthogonalized/normalized for stable repeated calls to `sample`.

Returns a named tuple with fields:
- `configs::Matrix{Int}`: top configurations in basis-index encoding (`top_kept x L`)
- `counts::Vector{Int}`: sample counts for each top configuration
- `probs::Vector{Float64}`: estimated probabilities (`counts / nsamples`)
- `stderrs::Vector{Float64}`: binomial standard error estimate for each probability
- `nsamples::Int`
- `unique_configurations::Int`
"""
function top_sampled_configurations(
    psi::MPS;
    nsamples::Int=100_000,
    top_k::Int=20,
    seed::Union{Nothing,Int}=nothing
)
    nsamples > 0 || throw(ArgumentError("nsamples must be positive, got $nsamples"))
    top_k > 0 || throw(ArgumentError("top_k must be positive, got $top_k"))

    psi_s = copy(psi)
    orthogonalize!(psi_s, 1)
    normalize!(psi_s)
    orthogonalize!(psi_s, 1)

    rng = seed === nothing ? Random.default_rng() : MersenneTwister(seed)

    L = length(psi_s)
    KeyT = NTuple{L,Int}
    counts_map = Dict{KeyT,Int}()

    for _ in 1:nsamples
        cfg = Tuple(sample(rng, psi_s))
        counts_map[cfg] = get(counts_map, cfg, 0) + 1
    end

    ranked = collect(counts_map)
    sort!(ranked; by=p -> last(p), rev=true)

    nkeep = min(top_k, length(ranked))
    configs = Matrix{Int}(undef, nkeep, L)
    counts = Vector{Int}(undef, nkeep)
    probs = Vector{Float64}(undef, nkeep)
    stderrs = Vector{Float64}(undef, nkeep)

    for row in 1:nkeep
        cfg = first(ranked[row])
        c = last(ranked[row])
        counts[row] = c

        p = c / nsamples
        probs[row] = p
        stderrs[row] = sqrt(p * (1 - p) / nsamples)

        for i in 1:L
            configs[row, i] = cfg[i]
        end
    end

    return (
        configs=configs,
        counts=counts,
        probs=probs,
        stderrs=stderrs,
        nsamples=nsamples,
        unique_configurations=length(ranked)
    )
end

@inline function _extract_nmaxb_from_tags(s::Index)
    for t in String.(tags(s))
        if startswith(t, "nmaxb=")
            return parse(Int, split(t, "=")[2])
        end
    end
    return nothing
end

"""
    infer_nmax_b_from_sites(sites) -> Int

Infer `nmax_b` from TwoBoson site tags (`nmaxb=...`).
"""
function infer_nmax_b_from_sites(sites)
    isempty(sites) && throw(ArgumentError("sites must be non-empty"))

    nmax_b = _extract_nmaxb_from_tags(sites[1])
    nmax_b === nothing && throw(ArgumentError("could not infer nmax_b from first site tags"))

    for i in 2:length(sites)
        nb_i = _extract_nmaxb_from_tags(sites[i])
        nb_i === nmax_b ||
            throw(ArgumentError("inconsistent nmax_b in site tags at site $i: got $(nb_i), expected $(nmax_b)"))
    end
    return nmax_b
end

"""
    decode_two_boson_configs(configs::AbstractMatrix{<:Integer}, nmax_b::Int)

Decode basis indices `k` (1-based linear index of `|na,nb>`) to occupancy tables.

Returns `(na, nb)` matrices with the same shape as `configs`.
"""
function decode_two_boson_configs(configs::AbstractMatrix{<:Integer}, nmax_b::Int)
    nmax_b >= 0 || throw(ArgumentError("nmax_b must be non-negative, got $nmax_b"))

    nr, nc = size(configs)
    na = Matrix{Int}(undef, nr, nc)
    nb = Matrix{Int}(undef, nr, nc)

    stride = nmax_b + 1
    for r in 1:nr, c in 1:nc
        k = Int(configs[r, c])
        k >= 1 || throw(ArgumentError("basis index must be >= 1, got $k"))

        k0 = k - 1
        na[r, c] = k0 ÷ stride
        nb[r, c] = k0 % stride
    end
    return na, nb
end
