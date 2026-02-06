# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

# -----------------------------------
# 6) Binding energy estimators
# -----------------------------------

"""
    effective_energy(energy::Real, Nb::Real, mu_b::Real)

Return the energy with the impurity chemical potential term removed.
For H = K + V - mu_a*Na - mu_b*Nb, this is:
    Ẽ = <H> + mu_b*<Nb>
"""
function effective_energy(energy::Real, Nb::Real, mu_b::Real)
    return energy
    # return energy + mu_b * Nb
end

"""
    sector_energy(; Nb_total::Int, mu_b::Real=0.0, keep_state::Bool=false, kwargs...)

Run DMRG in a fixed-impurity sector and return the effective energy Ẽ.
If keep_state=false, only returns scalars.
"""
function sector_energy(; Nb_total::Int, nmax_b_sector::Int=0, mu_b::Real=0.0, keep_state::Bool=false, kwargs...)
    energy, psi, sites, H = run_dmrg(; Nb_total=Nb_total, nmax_b=nmax_b_sector, mu_b=mu_b, kwargs...)
    na, nb = measure_densities(psi, sites)
    # Print average nb and na density
    println("Average Nb density = ", mean(nb))
    println("Average Na density = ", mean(na))

    Nb = sum(nb)
    Etilde = effective_energy(energy, Nb, mu_b)

    if keep_state
        return (; energy, Etilde, Nb, psi, sites, H)
    end
    return (; energy, Etilde, Nb)
end

"""
    binding_energies(; mu_b::Real=0.0, kwargs...)

Compute binding energy combinations from sectors Nb = 0,1,2,3.

Returns a NamedTuple with:
- E0..E3 (effective energies)
- Ebind2: E2 - 2E1 + E0
- Ebind3111: E3 - 3E1 + 2E0
- Ebind321: E3 - E2 - E1 + E0
"""
function binding_energies(; mu_b::Real=0.0, kwargs...)
    r0 = sector_energy(; Nb_total=0, nmax_b_sector=1, mu_b=mu_b, kwargs...)
    r1 = sector_energy(; Nb_total=1, nmax_b_sector=1, mu_b=mu_b, kwargs...)
    r2 = sector_energy(; Nb_total=2, nmax_b_sector=2, mu_b=mu_b, kwargs...)
    r3 = sector_energy(; Nb_total=3, nmax_b_sector=3, mu_b=mu_b, kwargs...)

    E0 = r0.Etilde
    E1 = r1.Etilde
    E2 = r2.Etilde
    E3 = r3.Etilde

    Ebind2 = E2 - 2 * E1 + E0
    Ebind3111 = E3 - 3 * E1 + 2 * E0
    Ebind321 = E3 - E2 - E1 + E0

    return (; E0, E1, E2, E3, Ebind2, Ebind3111, Ebind321)
end
