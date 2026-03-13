# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

# ----------------------------
# Single-particle density matrix
# ----------------------------

@inline function _spdm_ops(species::String)
    s = lowercase(species)
    if s == "a"
        return "Adag", "A"
    elseif s == "b"
        return "Bdag", "B"
    end
    throw(ArgumentError("species must be one of: \"a\", \"b\""))
end

"""
    single_particle_density_matrix_legacy(psi, sites, species::String; ishermitian::Bool=true)

Reference implementation of the single-particle density matrix
`G1[i,j] = ⟨bᵢ† bⱼ⟩` (or the corresponding `a`-species analogue),
evaluated using explicit MPO expectations.
"""
function single_particle_density_matrix_legacy(
    psi::MPS,
    sites,
    species::String;
    ishermitian::Bool=true
)
    L = length(sites)
    length(psi) == L ||
        throw(ArgumentError("sites length ($L) must match MPS length ($(length(psi)))"))

    create_op, annihilate_op = _spdm_ops(species)
    G1 = zeros(ComplexF64, L, L)

    for i in 1:L
        G1[i, i] = expect_product(psi, sites, [(create_op, i), (annihilate_op, i)])
        if ishermitian
            for j in (i + 1):L
                v = expect_product(psi, sites, [(create_op, i), (annihilate_op, j)])
                G1[i, j] = v
                G1[j, i] = conj(v)
            end
        else
            for j in (i + 1):L
                G1[i, j] = expect_product(psi, sites, [(create_op, i), (annihilate_op, j)])
            end
            for j in 1:(i - 1)
                G1[i, j] = expect_product(psi, sites, [(create_op, i), (annihilate_op, j)])
            end
        end
    end

    return G1
end

"""
    single_particle_density_matrix_correlation(psi, sites, species::String; ishermitian::Bool=true)

Fast implementation of the single-particle density matrix based on
`ITensorMPS.correlation_matrix`.
"""
function single_particle_density_matrix_correlation(
    psi::MPS,
    sites,
    species::String;
    ishermitian::Bool=true
)
    L = length(sites)
    length(psi) == L ||
        throw(ArgumentError("sites length ($L) must match MPS length ($(length(psi)))"))

    create_op, annihilate_op = _spdm_ops(species)
    site_range = 1:L
    return ComplexF64.(correlation_matrix(
        psi,
        create_op,
        annihilate_op;
        sites=site_range,
        ishermitian=ishermitian
    ))
end

"""
    single_particle_density_matrix(psi, sites, species::String;
                                   backend=:correlation_matrix,
                                   ishermitian::Bool=true)

Compute the species-resolved single-particle density matrix using the selected backend:
- `backend=:correlation_matrix` (default): fast MPS correlation routine
- `backend=:legacy`: MPO-per-pair reference implementation
"""
function single_particle_density_matrix(
    psi::MPS,
    sites,
    species::String;
    backend::Union{Symbol,AbstractString}=:correlation_matrix,
    ishermitian::Bool=true
)
    b = lowercase(String(backend))
    if b in ("correlation_matrix", "correlation", "fast")
        return single_particle_density_matrix_correlation(
            psi,
            sites,
            species;
            ishermitian=ishermitian
        )
    elseif b in ("legacy", "mpo", "opsum")
        return single_particle_density_matrix_legacy(
            psi,
            sites,
            species;
            ishermitian=ishermitian
        )
    end
    throw(ArgumentError(
        "backend must be one of :correlation_matrix/:legacy (also accepts correlation, fast, mpo, opsum); got \"$backend\""
    ))
end
