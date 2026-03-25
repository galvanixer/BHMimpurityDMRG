# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

# -----------------------------------
# 2) Initial product state constructor
# -----------------------------------

using Random

"""
    _require_distinct_impurity_sites(mode, L, Nb_total, nmax_b)

Validate the modes that place impurities on distinct sites rather than allowing
multiple impurities on the same site.

These modes are:
- `:centered_block`
- `:uniform_spread`
- `:random_separated`
- `:two_cluster`
- `:asymmetric_mixed`

They all require:
- `Nb_total <= L`, because each impurity starts on its own site
- `nmax_b >= 1` when `nmax_b` is provided

Example:
- `L = 32`, `Nb_total = 6` is valid
- `L = 4`, `Nb_total = 6` is invalid for these modes
"""
@inline function _require_distinct_impurity_sites(mode::Symbol, L::Int, Nb_total::Int, nmax_b)
    Nb_total <= L || throw(
        ArgumentError(
            "impurity_distribution=$mode requires Nb_total <= L because each impurity starts on a distinct site"
        )
    )
    if Nb_total > 0 && nmax_b !== nothing && nmax_b < 1
        throw(ArgumentError("impurity_distribution=$mode requires nmax_b >= 1"))
    end
    return nothing
end

"""
    _strictly_increasing_positions(targets, L)

Turn approximate target positions into a valid strictly increasing integer site list
inside `1:L`.

This helper is used by placement modes that start from "ideal" floating-point target
locations and then need a legal site list after rounding.

Example:
- targets near `[3.2, 8.1, 12.9, 18.4]` on `L = 20` become `[3, 8, 13, 18]`
- if rounding would push the last site beyond `L`, the whole pattern is shifted left
"""
function _strictly_increasing_positions(targets, L::Int)
    isempty(targets) && return Int[]

    pos = sort!(Int[round(Int, x) for x in targets])
    n = length(pos)
    # First clamp the left edge to a legal starting point. The upper clamp keeps
    # enough room for the remaining `n - 1` impurities.
    pos[1] = clamp(pos[1], 1, max(1, L - n + 1))
    for i in 2:n
        # Enforce distinct occupied sites after rounding.
        pos[i] = max(pos[i], pos[i - 1] + 1)
    end

    overflow = pos[end] - L
    if overflow > 0
        # If the right edge spills past the chain, shift the whole pattern left and
        # re-enforce strict ordering from right to left.
        pos .-= overflow
        for i in (n - 1):-1:1
            pos[i] = min(pos[i], pos[i + 1] - 1)
        end
    end

    pos[1] >= 1 || throw(ArgumentError("failed to place distinct impurity sites in a chain of length $L"))
    return pos
end

"""
    _nb_from_positions(L, positions)

Convert a list of occupied impurity sites into an `nb` occupation vector of length `L`.

Example:
- `positions = [4, 5, 6, 15, 24, 30]`
- result has `nb[4] = nb[5] = ... = nb[30] = 1` and zero elsewhere
"""
@inline function _nb_from_positions(L::Int, positions::AbstractVector{<:Integer})
    nb = zeros(Int, L)
    for i in positions
        nb[i] += 1
    end
    return nb
end

"""
    _centered_pileup_occupations(L, Nb_total)

Place all impurities on the single center site.

This is the strongest possible same-site clustered initial condition.

Examples:
- `L = 32`, `Nb_total = 6` gives `nb[16] = 6`
- `L = 11`, `Nb_total = 3` gives `nb[6] = 3`
"""
function _centered_pileup_occupations(L::Int, Nb_total::Int)
    nb = zeros(Int, L)
    Nb_total == 0 && return nb
    c = (L + 1) ÷ 2
    nb[c] = Nb_total
    return nb
end

"""
    _random_capped_occupations(L, Nb_total; nmax_b=nothing, seed=nothing)

Randomly place impurities on the chain, optionally respecting a per-site cap `nmax_b`.

This mode allows repeated occupation of the same site. It is the stochastic counterpart
to `:centered_pileup`, not a separated-site mode.

Examples:
- with `nmax_b = nothing`, repeated hits on the same site are allowed without limit
- with `nmax_b = 2`, no site can start with more than two impurities
- with a fixed `seed`, the placement is reproducible
"""
function _random_capped_occupations(L::Int, Nb_total::Int; nmax_b::Union{Int,Nothing}=nothing, seed::Union{Int,Nothing}=nothing)
    nb = zeros(Int, L)
    Nb_total == 0 && return nb

    if seed !== nothing
        Random.seed!(seed)
    end

    remaining = Nb_total
    if nmax_b === nothing
        # No cap: drop each impurity on an independently sampled site.
        for _ in 1:remaining
            nb[rand(1:L)] += 1
        end
    else
        # Respect per-site cap nmax_b while preserving the existing rejection-sampling behavior.
        while remaining > 0
            i = rand(1:L)
            if nb[i] < nmax_b
                nb[i] += 1
                remaining -= 1
            end
        end
    end
    return nb
end

"""
    _centered_block_sites(L, Nb_total)

Return a contiguous centered block of distinct impurity sites.

Example:
- `L = 32`, `Nb_total = 6` gives `[14, 15, 16, 17, 18, 19]`
- `L = 12`, `Nb_total = 3` gives `[5, 6, 7]`
"""
function _centered_block_sites(L::Int, Nb_total::Int)
    Nb_total == 0 && return Int[]
    start = (L - Nb_total) ÷ 2 + 1
    return collect(start:(start + Nb_total - 1))
end

"""
    _uniform_spread_sites(L, Nb_total)

Return a deterministic set of distinct impurity sites with near-uniform spacing.

The construction is meant to give a clean "spread out" competitor to the centered
clustered states.

Example:
- `L = 32`, `Nb_total = 6` gives `[3, 8, 13, 18, 23, 28]`
"""
function _uniform_spread_sites(L::Int, Nb_total::Int)
    Nb_total == 0 && return Int[]

    step = max(1, fld(L, Nb_total))
    targets = if step == 1
        # In crowded cases use fractional targets across the whole system, then repair
        # them into legal distinct sites.
        [round(Int, j * (L + 1) / (Nb_total + 1)) for j in 1:Nb_total]
    else
        # In the sparse case use a simple arithmetic progression centered away from the
        # very edge of the chain.
        first = fld(L, 2 * Nb_total) + 1
        [first + step * (j - 1) for j in 1:Nb_total]
    end
    return _strictly_increasing_positions(targets, L)
end

"""
    _random_separated_sites(L, Nb_total; seed=nothing)

Return a random set of distinct impurity sites with an enforced minimum spacing.

The spacing target is derived from `fld(L, Nb_total)`, so this mode produces random
placements that are still biased away from immediate clustering.

Example:
- `L = 32`, `Nb_total = 6`, `seed = 123` gives a reproducible separated pattern
- unlike `:random_capped`, every impurity starts on a distinct site
"""
function _random_separated_sites(L::Int, Nb_total::Int; seed::Union{Int,Nothing}=nothing)
    Nb_total == 0 && return Int[]

    min_gap = max(1, fld(L, Nb_total))
    reduced_length = L - (Nb_total - 1) * (min_gap - 1)
    reduced_length >= Nb_total || throw(
        ArgumentError(
            "impurity_distribution=random_separated cannot fit Nb_total=$Nb_total distinct impurities on L=$L sites"
        )
    )

    rng = seed === nothing ? Random.default_rng() : Random.MersenneTwister(seed)
    reduced_sites = randperm(rng, reduced_length)[1:Nb_total]
    sort!(reduced_sites)
    return [reduced_sites[i] + (i - 1) * (min_gap - 1) for i in 1:Nb_total]
end

"""
    _two_cluster_sites(L, Nb_total)

Split impurities into two compact separated clusters.

Examples:
- `L = 32`, `Nb_total = 6` gives `[9, 10, 11, 22, 23, 24]`
- `L = 32`, `Nb_total = 5` gives two clusters of sizes `2` and `3`
"""
function _two_cluster_sites(L::Int, Nb_total::Int)
    Nb_total <= 1 && return _centered_block_sites(L, Nb_total)

    left_count = fld(Nb_total, 2)
    right_count = Nb_total - left_count
    empty_sites = L - Nb_total
    side_gap = fld(empty_sites, 3)
    middle_gap = empty_sites - 2 * side_gap

    left_start = side_gap + 1
    right_start = left_start + left_count + middle_gap
    positions = vcat(
        collect(left_start:(left_start + left_count - 1)),
        collect(right_start:(right_start + right_count - 1))
    )
    return positions
end

"""
    _asymmetric_mixed_sites(L, Nb_total)

Build a deterministic, intentionally irregular impurity pattern.

The construction starts with a small left-side cluster and then places the remaining
impurities as unevenly spaced singles biased toward the right half of the chain.

Example:
- `L = 32`, `Nb_total = 6` gives `[4, 5, 6, 15, 24, 30]`
"""
function _asymmetric_mixed_sites(L::Int, Nb_total::Int)
    Nb_total == 0 && return Int[]

    cluster_count = min(Nb_total, max(1, cld(Nb_total, 2)))
    cluster_start = min(max(1, fld(L, 8)), max(1, L - Nb_total + 1))
    cluster_sites = collect(cluster_start:(cluster_start + cluster_count - 1))

    tail_count = Nb_total - cluster_count
    tail_count == 0 && return cluster_sites

    tail_start = min(L, cluster_sites[end] + max(2, fld(L, 4)) + 1)
    tail_stop = max(tail_start, L - max(2, fld(L, 16)))
    tail_targets = if tail_count == 1
        [tail_stop]
    else
        span = tail_stop - tail_start
        [tail_start + span * ((j - 1) / (tail_count - 1))^0.7 for j in 1:tail_count]
    end

    return _strictly_increasing_positions(vcat(cluster_sites, tail_targets), L)
end

"""
    _distinct_impurity_positions(mode, L, Nb_total; nmax_b=nothing, seed=nothing)

Dispatch helper for the distinct-site impurity modes.

This returns site positions such as `[14, 15, 16, 17, 18, 19]`, not the full `nb`
occupation vector.
"""
function _distinct_impurity_positions(
    mode::Symbol,
    L::Int,
    Nb_total::Int;
    nmax_b::Union{Int,Nothing}=nothing,
    seed::Union{Int,Nothing}=nothing
)
    _require_distinct_impurity_sites(mode, L, Nb_total, nmax_b)

    if mode == :centered_block
        return _centered_block_sites(L, Nb_total)
    elseif mode == :uniform_spread
        return _uniform_spread_sites(L, Nb_total)
    elseif mode == :random_separated
        return _random_separated_sites(L, Nb_total; seed=seed)
    elseif mode == :two_cluster
        return _two_cluster_sites(L, Nb_total)
    elseif mode == :asymmetric_mixed
        return _asymmetric_mixed_sites(L, Nb_total)
    end
    error("internal error: unsupported distinct impurity mode $mode")
end

"""
    _impurity_occupations(mode, L, Nb_total; nmax_b=nothing, seed=nothing)

Top-level impurity initializer for all supported modes.

This is the single place where the public `impurity_distribution` mode names are
interpreted. It returns the full `nb` occupation vector used by
`initial_configuration`.
"""
function _impurity_occupations(
    mode::Symbol,
    L::Int,
    Nb_total::Int;
    nmax_b::Union{Int,Nothing}=nothing,
    seed::Union{Int,Nothing}=nothing
)
    if mode == :centered_pileup
        return _centered_pileup_occupations(L, Nb_total)
    elseif mode == :random_capped
        return _random_capped_occupations(L, Nb_total; nmax_b=nmax_b, seed=seed)
    elseif mode in (:centered_block, :uniform_spread, :random_separated, :two_cluster, :asymmetric_mixed)
        positions = _distinct_impurity_positions(mode, L, Nb_total; nmax_b=nmax_b, seed=seed)
        return _nb_from_positions(L, positions)
    end

    throw(
        ArgumentError(
            "impurity_distribution must be one of :centered_pileup, :random_capped, :centered_block, :uniform_spread, :random_separated, :two_cluster, or :asymmetric_mixed"
        )
    )
end

"""
    initial_configuration(L::Int; Na_total::Int, Nb_total::Int,
                          impurity_distribution::Symbol=:centered_pileup,
                          nmax_a::Union{Int,Nothing}=nothing,
                          nmax_b::Union{Int,Nothing}=nothing,
                          seed::Union{Int,Nothing}=nothing)

Generate an initial configuration for a lattice system with `L` sites, distributing two
types of particles: bath particles (`Na_total`) and impurity particles (`Nb_total`).

`impurity_distribution` can be:
- `:centered_pileup` (default): place all impurities at center site.
- `:random_capped`: randomly distribute impurities across sites, respecting `nmax_b` if provided.
- `:centered_block`: occupy a contiguous centered block of `Nb_total` distinct sites.
- `:uniform_spread`: place impurities on distinct sites with near-uniform spacing.
- `:random_separated`: random distinct sites with a built-in minimum spacing.
- `:two_cluster`: split impurities into two compact separated clusters.
- `:asymmetric_mixed`: use a deterministic irregular pattern with a left cluster and right-biased singles.

Bath particles are initialized separately from impurities:
- species `a` is distributed approximately uniformly across the chain
- species `b` is assigned using the selected impurity mode

Examples for `L = 32`, `Nb_total = 6`:
- `:centered_pileup` -> site `16` gets all 6 impurities
- `:centered_block` -> `[14, 15, 16, 17, 18, 19]`
- `:uniform_spread` -> `[3, 8, 13, 18, 23, 28]`
- `:two_cluster` -> `[9, 10, 11, 22, 23, 24]`
- `:asymmetric_mixed` -> `[4, 5, 6, 15, 24, 30]`

Examples for the return value:
- if `Na_total = 4`, `Nb_total = 2`, `L = 4`, the output is a vector of
  `(na, nb)` tuples, one per site
- a possible result is `[(1, 0), (1, 2), (1, 0), (1, 0)]` for `:centered_pileup`
"""
function initial_configuration(L::Int; Na_total::Int, Nb_total::Int,
    impurity_distribution::Symbol=:centered_pileup,
    nmax_a::Union{Int,Nothing}=nothing,
    nmax_b::Union{Int,Nothing}=nothing,
    seed::Union{Int,Nothing}=nothing)
    na = fill(0, L)

    # Species-a particles are not treated as special "impurities" here. They are
    # simply spread as evenly as possible over the chain before the species-b pattern
    # is overlaid.
    base = Na_total ÷ L
    rem = Na_total % L
    for i in 1:L
        na[i] = base + (i <= rem ? 1 : 0)
    end
    if nmax_a !== nothing && maximum(na) > nmax_a
        throw(ArgumentError("Na_total too large for nmax_a=$nmax_a (max site occupancy=$(maximum(na)))"))
    end

    # Species-b occupations come entirely from the chosen impurity mode.
    nb = _impurity_occupations(
        impurity_distribution,
        L,
        Nb_total;
        nmax_b=nmax_b,
        seed=seed
    )

    return collect(zip(na, nb))
end

"""
    product_state_mps(sites, conf::Vector{Tuple{Int,Int}}; nmax_b::Int)

Convert a configuration [(na, nb), ...] into a product MPS.

Example:
- `conf = [(1, 0), (1, 0), (1, 2), (1, 0)]` means site 3 starts with two impurities
- each tuple is converted into the corresponding local basis state `|na, nb>`
"""
function product_state_mps(sites, conf::Vector{Tuple{Int,Int}}; nmax_b::Int)
    L = length(sites)
    psi = MPS(sites)

    for i in eachindex(sites)
        s = sites[i]
        na, nb = conf[i]
        # `_state_index` maps the local occupation pair `(na, nb)` into the linear
        # basis ordering used by the custom `TwoBoson` site type.
        k = _state_index(na, nb, nmax_b)
        A = ITensor(s)
        A[s => k] = 1.0
        psi[i] = A
    end
    return psi
end
