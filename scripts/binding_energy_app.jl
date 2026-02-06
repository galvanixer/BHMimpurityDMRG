# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

using Statistics

include(joinpath(@__DIR__, "..", "src", "BHMimpurityDMRG.jl"))
using .BHMimpurityDMRG

function main()
    params_path = get(ENV, "PARAMS", joinpath(@__DIR__, "..", "parameters.yaml"))
    cfg = isfile(params_path) ? load_params(params_path) : Dict{String,Any}()
    dmrg_cfg = haskey(cfg, "dmrg") ? dict_to_namedtuple(cfg["dmrg"]) : (;)
    bind_cfg = haskey(cfg, "binding_energy") ? dict_to_namedtuple(cfg["binding_energy"]) : (;)

    # Example run for binding energies (Nb = 0..3)
    res = binding_energies(; dmrg_cfg..., bind_cfg...)

    println("E0 = ", res.E0, " E0/L = ", res.E0 / 12)
    println("E1 = ", res.E1, " E1/L = ", res.E1 / 12)
    println("E2 = ", res.E2, " E2/L = ", res.E2 / 12)
    println("E3 = ", res.E3, " E3/L = ", res.E3 / 12)
    println("Ebind2 (dimer) = ", res.Ebind2)
    println("Ebind3111 (trimer vs 3 polarons) = ", res.Ebind3111)
    println("Ebind321 (trimer vs dimer+monomer) = ", res.Ebind321)

    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
