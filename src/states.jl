# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

# -----------------------------------
# 2) Initial product state constructor
# -----------------------------------

using Random

"""
    initial_configuration(L::Int; Na_total::Int, Nb_total::Int,
                          impurity_distribution::Symbol=:center,
                          nmax_a::Union{Int,Nothing}=nothing,
                          nmax_b::Union{Int,Nothing}=nothing,
                          seed::Union{Int,Nothing}=nothing)

Generate an initial configuration for a lattice system with `L` sites, distributing two
types of particles: bath particles (`Na_total`) and impurity particles (`Nb_total`).

`impurity_distribution` can be:
- `:center` (default): place all impurities at center site.
- `:random`: randomly distribute impurities across sites, respecting `nmax_b` if provided.
"""
function initial_configuration(L::Int; Na_total::Int, Nb_total::Int,
    impurity_distribution::Symbol=:center,
    nmax_a::Union{Int,Nothing}=nothing,
    nmax_b::Union{Int,Nothing}=nothing,
    seed::Union{Int,Nothing}=nothing)
    na = fill(0, L)
    nb = fill(0, L)

    # distribute bath approximately uniformly
    base = Na_total ÷ L
    rem = Na_total % L
    for i in 1:L
        na[i] = base + (i <= rem ? 1 : 0)
    end
    if nmax_a !== nothing && maximum(na) > nmax_a
        throw(ArgumentError("Na_total too large for nmax_a=$nmax_a (max site occupancy=$(maximum(na)))"))
    end

    if impurity_distribution == :center
        # put impurities in the center
        c = (L + 1) ÷ 2
        nb[c] = Nb_total
    elseif impurity_distribution == :random
        if seed !== nothing
            Random.seed!(seed)
        end
        remaining = Nb_total
        if nmax_b === nothing
            # no cap, just drop each impurity on a random site
            for _ in 1:remaining
                nb[rand(1:L)] += 1
            end
        else
            # respect per-site cap nmax_b
            while remaining > 0
                i = rand(1:L)
                if nb[i] < nmax_b
                    nb[i] += 1
                    remaining -= 1
                end
            end
        end
    else
        throw(ArgumentError("impurity_distribution must be :center or :random"))
    end

    return collect(zip(na, nb))
end

"""
    product_state_mps(sites, conf::Vector{Tuple{Int,Int}}; nmax_b::Int)

Convert a configuration [(na, nb), ...] into a product MPS.
"""
function product_state_mps(sites, conf::Vector{Tuple{Int,Int}}; nmax_b::Int)
    L = length(sites)
    psi = MPS(sites)

    for i in eachindex(sites)
        s = sites[i]
        na, nb = conf[i]
        k = _state_index(na, nb, nmax_b)
        A = ITensor(s)
        A[s => k] = 1.0
        psi[i] = A
    end
    return psi
end
