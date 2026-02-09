# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

using Statistics

include(joinpath(@__DIR__, "..", "src", "BHMimpurityDMRG.jl"))
using .BHMimpurityDMRG
using HDF5
using Logging

function main()
    params_path = get(ENV, "PARAMS", joinpath(@__DIR__, "..", "configs", "parameters.yaml"))
    cfg = isfile(params_path) ? load_params(params_path) : Dict{String,Any}()
    dmrg_cfg = merge_sections(cfg, ["lattice", "local_hilbert", "initial_state", "hamiltonian", "dmrg"])
    bind_cfg = haskey(cfg, "binding_energy") ? dict_to_namedtuple(cfg["binding_energy"]) : (;)

    io_cfg = get(cfg, "io", Dict{String,Any}())
    results_path = get(io_cfg, "results_path", "results.h5")

    log_path = get(io_cfg, "log_path", "run.log")
    logio = open(log_path, "a")
    logger = SimpleLogger(logio, Logging.Info)

    with_logger(logger) do
        @info "Starting binding_energy" params_path=params_path results_path=results_path
        # Example run for binding energies (Nb = 0..3)
        res = binding_energies(; dmrg_cfg..., bind_cfg...)

        mode = isfile(results_path) ? "r+" : "w"
        HDF5.h5open(results_path, mode) do f
            g_meta = ensure_group(f, "meta")
            write_meta!(g_meta; params_path=params_path)
            g_obs = ensure_group(f, "observables")
            g_be = ensure_group(g_obs, "binding_energy")
            write_or_replace(g_be, "E0", res.E0)
            write_or_replace(g_be, "E1", res.E1)
            write_or_replace(g_be, "E2", res.E2)
            write_or_replace(g_be, "E3", res.E3)
            write_or_replace(g_be, "Ebind2", res.Ebind2)
            write_or_replace(g_be, "Ebind3111", res.Ebind3111)
            write_or_replace(g_be, "Ebind321", res.Ebind321)
        end
        @info "Wrote results" results_path=results_path
    end

    close(logio)

    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
