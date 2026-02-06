# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

include(joinpath(@__DIR__, "..", "src", "BHMimpurityDMRG.jl"))
using .BHMimpurityDMRG

function main()
    params_path = get(ENV, "PARAMS", joinpath(@__DIR__, "..", "parameters.yaml"))
    cfg = isfile(params_path) ? load_params(params_path) : Dict{String,Any}()
    dmrg_cfg = haskey(cfg, "dmrg") ? dict_to_namedtuple(cfg["dmrg"]) : (;)

    # Example run: single impurity (Nb_total=1), bath ~ unit filling (Na_total=L)
    energy, psi, sites, H = run_dmrg(; dmrg_cfg...)

    println("Ground-state energy = ", energy)

    na, nb = measure_densities(psi, sites)
    println("Na (first 10 sites) = ", na[1:min(end, 10)])
    println("Nb (first 10 sites) = ", nb[1:min(end, 10)])

    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
