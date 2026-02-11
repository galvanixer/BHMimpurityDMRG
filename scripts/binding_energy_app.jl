# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

using Statistics

include(joinpath(@__DIR__, "..", "src", "BHMimpurityDMRG.jl"))
using .BHMimpurityDMRG
using HDF5
using Logging

function main()
    params_path = length(ARGS) >= 1 ? ARGS[1] :
                  get(ENV, "PARAMS", joinpath(@__DIR__, "..", "configs", "binding_energy.yaml"))

    # Load parameters with fallback to empty dict if file not found or loading fails
    cfg = isfile(params_path) ? load_params(params_path) : Dict{String,Any}()
    dmrg_cfg = merge_sections(cfg, ["lattice", "local_hilbert", "initial_state", "hamiltonian", "dmrg"])
    if hasproperty(dmrg_cfg, :Nb_total)
        dmrg_cfg = Base.structdiff(dmrg_cfg, (; Nb_total=nothing))
    end
    bind_cfg = get(cfg, "binding_energy", Dict{String,Any}())

    io_cfg = get(cfg, "io", Dict{String,Any}())
    results_path = get(io_cfg, "results_path", "results.h5")

    logcfg = setup_logger(io_cfg; default_log_path="run.log")
    logger = logcfg.logger

    with_logger(logger) do
        @info "Starting binding_energy estimation" params_path=params_path results_path=results_path
        # Binding energy over specified sectors
        sectors = get(bind_cfg, "sectors", nothing)
        save_states = get(bind_cfg, "save_states", false)
        state_tpl = get(bind_cfg, "state_path_template", "dmrg_state_Nb{Nb_total}.h5")
        if sectors === nothing
            res = binding_energies(; checkpoint_params_path=params_path, dmrg_cfg...)
        else
            results = Dict{Int, Any}()
            for s in sectors
                nb = if s isa Integer
                    Int(s)
                elseif s isa AbstractDict
                    nb_val = haskey(s, "Nb_total") ? s["Nb_total"] :
                             haskey(s, :Nb_total) ? s[:Nb_total] :
                             error("binding_energy.sectors entries must include Nb_total")
                    Int(nb_val)
                else
                    error("binding_energy.sectors entries must be integers or dicts with Nb_total")
                end
                nmax_b = if s isa AbstractDict
                    haskey(s, "nmax_b") ? s["nmax_b"] :
                    haskey(s, :nmax_b) ? s[:nmax_b] :
                    get(get(cfg, "local_hilbert", Dict{String,Any}()), "nmax_b", nothing)
                else
                    get(get(cfg, "local_hilbert", Dict{String,Any}()), "nmax_b", nothing)
                end
                @info "Running DMRG for sector" Nb_total=nb nmax_b=nmax_b save_state=save_states
                r = sector_energy(; Nb_total_sector=nb, nmax_b_sector=nmax_b, keep_state=save_states, checkpoint_params_path=params_path, dmrg_cfg...)
                
                # Optionally save state for this sector
                if save_states && haskey(r, :psi)
                    state_path = replace(state_tpl, "{Nb_total}" => string(nb))
                    na, nb_dens = measure_densities(r.psi, r.sites)
                    _, init_na, init_nb = dmrg_initial_configuration(;
                        Nb_total=nb,
                        nmax_b=nmax_b,
                        dmrg_cfg...
                    )
                    save_state(
                        state_path,
                        r.psi;
                        energy=r.energy,
                        params_path=params_path,
                        na=na,
                        nb=nb_dens,
                        init_na=init_na,
                        init_nb=init_nb
                    )
                    @info "Saved state" state_path=state_path Nb_total=nb
                end
                results[nb] = r
            end
            E0 = results[0].Etilde
            E1 = results[1].Etilde
            E2 = results[2].Etilde
            E3 = results[3].Etilde
            res = (; E0, E1, E2, E3,
                Ebind2 = E2 - 2 * E1 + E0,
                Ebind3111 = E3 - 3 * E1 + 2 * E0,
                Ebind321 = E3 - E2 - E1 + E0)
        end

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

    close(logcfg.logio)

    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
