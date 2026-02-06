# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

# ----------------------------
# Simple densities
# ----------------------------

function onsite_expect(psi::MPS, sites, opname::String, i::Int)
    orthogonalize!(psi, i)

    s = sites[i]
    Ai = psi[i]

    Oi = op(opname, s)

    # Build bra tensor with correct priming to match Oi
    braAi = dag(prime(Ai, s))
    ketAi = Ai

    return (braAi * Oi * ketAi)[]
end

function measure_densities(psi, sites)
    L = length(sites)
    na = zeros(L)
    nb = zeros(L)
    for i in 1:L
        na[i] = onsite_expect(psi, sites, "Na", i)
        nb[i] = onsite_expect(psi, sites, "Nb", i)
    end
    return na, nb
end
