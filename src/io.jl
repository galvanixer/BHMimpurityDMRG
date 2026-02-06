# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

"""
    save_state(path::AbstractString, psi::MPS; energy=nothing, sites=siteinds(psi))

Save the ground state `psi` (and optionally `energy` and `sites`) to an HDF5 file.

Requires `HDF5.jl` to be installed:
    import Pkg; Pkg.add("HDF5")
"""
function save_state(path::AbstractString, psi::MPS; energy=nothing, sites=siteinds(psi))
    @eval using HDF5
    HDF5.h5open(path, "w") do f
        write(f, "psi", psi)
        write(f, "sites", sites)
        if energy !== nothing
            write(f, "energy", energy)
        end
    end
    return nothing
end

"""
    load_state(path::AbstractString)

Load a saved MPS ground state from an HDF5 file.

Returns a NamedTuple `(psi, sites, energy)` where `energy` may be `nothing`
if it was not stored.

Requires `HDF5.jl` to be installed.
"""
function load_state(path::AbstractString)
    @eval using HDF5
    HDF5.h5open(path, "r") do f
        psi = read(f, "psi", MPS)
        sites = read(f, "sites", Vector{Index})
        energy = haskey(f, "energy") ? read(f, "energy") : nothing
        return (; psi, sites, energy)
    end
end
