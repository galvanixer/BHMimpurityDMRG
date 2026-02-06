# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

# Common measurement utilities used by observables.

"""
    expect_product(psi::MPS, sites, ops::Vector{Tuple{String,Int}})
Compute the expectation value of a product of on-site operators specified by `ops`.

# Arguments
- `psi::MPS`: The MPS state on which to measure.
- `sites`: The site indices corresponding to the MPS.
- `ops::Vector{Tuple{String,Int}}`: A vector of tuples where each tuple is `(opname, i)` specifying the operator name (e.g. "Na", "Nb") and the site index `i` on which it acts.

# Returns
- The expectation value ⟨O⟩ where O is the product of the specified operators.

# Details
- Constructs an `MPO` representing the product of the specified operators and applies it to `psi` to compute the expectation value.
"""
function expect_product(psi::MPS, sites, ops::Vector{Tuple{String,Int}})
    isempty(ops) && throw(ArgumentError("ops must be non-empty"))
    args = Any[]
    for (opname, i) in ops
        push!(args, opname)
        push!(args, i)
    end
    os = OpSum()
    os += 1.0, args...
    O = MPO(os, sites)

    return inner(psi, Apply(O, psi))
end


"""
    expect_n(psi::MPS, sites, opname::String, i::Int)

Compute the expectation value of a single operator at site `i`.

# Arguments
- `psi::MPS`: The many-body wavefunction
- `sites`: The site indices
- `opname::String`: The name of the operator to measure
- `i::Int`: The site index at which to measure the operator

# Returns
The expectation value of the operator `opname` at site `i`.

# See Also
- [`expect_product`](@ref): For computing expectation values of products of operators
"""
expect_n(psi::MPS, sites, opname::String, i::Int) =
    expect_product(psi, sites, [(opname, i)])

"""
    expect_nn(psi::MPS, sites, opname::String, i::Int, j::Int)

Calculate the expectation value of a product of two operators acting on sites i and j.

# Arguments
- `psi::MPS`: The matrix product state (wavefunction)
- `sites`: The site indices for the MPS
- `opname::String`: The name of the operator to apply (e.g., "Sz", "N")
- `i::Int`: The first site index
- `j::Int`: The second site index

# Returns
- The expectation value of the product operator at sites i and j

# See Also
- [`expect_product`](@ref): For computing expectation values of products of operators
"""
expect_nn(psi::MPS, sites, opname::String, i::Int, j::Int) =
    expect_product(psi, sites, [(opname, i), (opname, j)])


"""
    expect_nnn(psi::MPS, sites, opname::String, i::Int, j::Int, k::Int)

Calculate the expectation value of a product of three operators at positions i, j, and k.

# Arguments
- `psi::MPS`: The matrix product state
- `sites`: The sites object containing site information
- `opname::String`: The name of the operator to apply at each position
- `i::Int`: The first site position
- `j::Int`: The second site position
- `k::Int`: The third site position

# Returns
The expectation value `⟨ψ | O_i O_j O_k | ψ⟩` where O is the operator specified by `opname`.

# See Also
- [`expect_product`](@ref): For computing expectation values of products of operators
"""
expect_nnn(psi::MPS, sites, opname::String, i::Int, j::Int, k::Int) =
    expect_product(psi, sites, [(opname, i), (opname, j), (opname, k)])
