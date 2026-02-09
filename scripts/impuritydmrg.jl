# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

include(joinpath(@__DIR__, "..", "src", "BHMimpurityDMRG.jl"))
using .BHMimpurityDMRG
using HDF5
using Logging

function main()
    params_path = get(ENV, "PARAMS", joinpath(@__DIR__, "..", "configs", "parameters.yaml"))
    cfg = isfile(params_path) ? load_params(params_path) : Dict{String,Any}()
    dmrg_cfg = merge_sections(cfg, ["lattice", "local_hilbert", "initial_state", "hamiltonian", "dmrg"])
    io_cfg = get(cfg, "io", Dict{String,Any}())
    results_path = get(io_cfg, "results_path", "results.h5")
    log_path = get(io_cfg, "log_path", "run.log")

    logio = open(log_path, "a")
    logger = SimpleLogger(logio, Logging.Info)

    with_logger(logger) do
        @info "Starting impuritydmrg" params_path=params_path results_path=results_path
        energy, psi, sites, H = run_dmrg(; dmrg_cfg...)

        na, nb = measure_densities(psi, sites)
        mode = isfile(results_path) ? "r+" : "w"
        HDF5.h5open(results_path, mode) do f
            g_meta = ensure_group(f, "meta")
            write_meta!(g_meta; params_path=params_path)
            g_obs = ensure_group(f, "observables")
            g_den = ensure_group(g_obs, "densities")
            write_or_replace(g_den, "energy", energy)
            write_or_replace(g_den, "na", na)
            write_or_replace(g_den, "nb", nb)
        end
        @info "Wrote results" results_path=results_path
    end

    close(logio)

    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
