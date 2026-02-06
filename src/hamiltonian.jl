# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

# ----------------------------
# 3) Build your Hamiltonian MPO
# ----------------------------

"""
    build_hamiltonian(sites; t_a::Real, t_b::Real, U_a::Real, U_b::Real=0.0,
                      U_ab::Real, mu_a::Real, mu_b::Real, periodic::Bool=false)

Construct the Hamiltonian MPO for a two-species Bose-Hubbard model with impurity terms:

H = -t_a Σ (adag_i a_{i+1} + h.c.)
    -t_b Σ (bdag_i b_{i+1} + h.c.)
    + (U_a/2) Σ Na_i (Na_i - 1)
    + (U_b/2) Σ Nb_i (Nb_i - 1)
    + U_ab Σ Na_i Nb_i
    - mu_a Σ Na_i - mu_b Σ Nb_i

Notes:
- `U_b` defaults to 0.0 so you can omit B-B interactions if desired.
- For backward compatibility, a second method accepts old keyword names
  (`U`, `UIB`, `muB`, `muI`) and forwards to this one.
"""
function build_hamiltonian(sites; t_a::Real, t_b::Real, U_a::Real, U_b::Real=0.0,
    U_ab::Real, mu_a::Real, mu_b::Real, periodic::Bool=false)
    L = length(sites)
    os = OpSum()

    # Hopping terms (Hermitian)
    for i in 1:(L - 1)
        os += -t_a, "Adag", i, "A", i + 1
        os += -t_a, "Adag", i + 1, "A", i
        os += -t_b, "Bdag", i, "B", i + 1
        os += -t_b, "Bdag", i + 1, "B", i
    end

    # PBC wrap-around bond: L <-> 1
    if periodic && L > 2
        os += -t_a, "Adag", L, "A", 1
        os += -t_a, "Adag", 1, "A", L
        os += -t_b, "Bdag", L, "B", 1
        os += -t_b, "Bdag", 1, "B", L
    end

    # Onsite terms
    for i in 1:L
        # (U_a/2) * Na*(Na-1) = (U_a/2) * (Na*Na - Na)
        os += 0.5 * U_a, "Na", i, "Na", i
        os += -0.5 * U_a, "Na", i

        # (U_b/2) * Nb*(Nb-1) = (U_b/2) * (Nb*Nb - Nb)
        os += 0.5 * U_b, "Nb", i, "Nb", i
        os += -0.5 * U_b, "Nb", i

        # U_ab * Na * Nb
        os += U_ab, "Na", i, "Nb", i

        # -mu * Na and -mu * Nb
        os += -mu_a, "Na", i
        os += -mu_b, "Nb", i
    end

    return MPO(os, sites)
end