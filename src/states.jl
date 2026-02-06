# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

# -----------------------------------
# 2) Initial product state constructor
# -----------------------------------

"""
    initial_configuration(L::Int; Na_total::Int, Nb_total::Int)

Generate an initial configuration for a lattice system with `L` sites, distributing two
types of particles: bath particles (`Na_total`) and impurity particles (`Nb_total`).
"""
function initial_configuration(L::Int; Na_total::Int, Nb_total::Int)
    na = fill(0, L)
    nb = fill(0, L)

    # distribute bath approximately uniformly
    base = Na_total ÷ L
    rem = Na_total % L
    for i in 1:L
        na[i] = base + (i <= rem ? 1 : 0)
    end

    # put impurities in the center (or spread if you want)
    c = (L + 1) ÷ 2
    nb[c] = Nb_total

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
