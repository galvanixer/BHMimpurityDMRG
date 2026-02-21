# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

include(joinpath(@__DIR__, "..", "src", "BHMimpurityDMRG.jl"))
using .BHMimpurityDMRG
using HDF5
using Logging
using SHA

function run_dmrg_to_log(log_path::AbstractString; kwargs...)
    open(log_path, "a") do io
        return redirect_stdout(io) do
            redirect_stderr(io) do
                return run_dmrg(; kwargs...)
            end
        end
    end
end

function main()
    params_path = length(ARGS) >= 1 ? ARGS[1] :
                  get(ENV, "PARAMS", joinpath(@__DIR__, "..", "configs", "parameters.yaml"))
    observables_path_arg = length(ARGS) >= 2 ? ARGS[2] : nothing

    # Load parameters with fallback to empty dict if file not found or loading fails
    cfg_base = isfile(params_path) ? load_params(params_path) : Dict{String,Any}()
    merged_cfg = with_observables_config(cfg_base; observables_path=observables_path_arg)
    cfg = merged_cfg.cfg
    observables_path = merged_cfg.observables_path
    observables_loaded = merged_cfg.observables_loaded
    params_has_observables = haskey(cfg_base, "observables") || haskey(cfg_base, :observables)

    dmrg_cfg = merge_sections(cfg_base, ["lattice", "local_hilbert", "initial_state", "hamiltonian", "dmrg"])
    io_cfg = get(cfg_base, "io", Dict{String,Any}())
    results_path = get(io_cfg, "results_path", "results.h5")
    state_path = get(io_cfg, "state_save_path", "dmrg_state.h5")
    save_state_flag = get(io_cfg, "save_state", true)
    logcfg = setup_logger(io_cfg; default_log_path="run.log")
    logger = logcfg.logger

    with_logger(logger) do
        @info "Solving for ground state using DMRG" params_path=params_path observables_path=observables_path observables_loaded=observables_loaded results_path=results_path
        if !observables_loaded && params_has_observables
            @info "Observables file not found; falling back to observables in parameters config" observables_path=observables_path
        elseif !observables_loaded
            @warn "Observables file not found and no observables section in parameters config; defaults will be used" observables_path=observables_path
        end
        params_text = isfile(params_path) ? read(params_path, String) : nothing
        current_hash = params_text === nothing ? nothing : bytes2hex(SHA.sha256(params_text))
        _, init_na, init_nb = dmrg_initial_configuration(; dmrg_cfg...)

        na = nothing
        nb = nothing
        if isfile(state_path)
            st = load_state(state_path)
            if current_hash !== nothing && st.params_sha256 !== nothing &&
               st.params_sha256 == current_hash && st.energy !== nothing
                psi = st.psi
                sites = st.sites
                energy = st.energy
                na = st.na
                nb = st.nb
                @info "Using cached state (hash matched)" state_path=state_path
            else
                @info "Cached state missing/mismatched hash; running DMRG" state_path=state_path
                energy, psi, sites, _ = run_dmrg_to_log(
                    logcfg.log_path;
                    checkpoint_params_path=params_path,
                    dmrg_cfg...
                )
                if save_state_flag
                    na_tmp, nb_tmp = measure_densities(psi, sites)
                    save_state(
                        state_path,
                        psi;
                        energy=energy,
                        params_path=params_path,
                        na=na_tmp,
                        nb=nb_tmp,
                        init_na=init_na,
                        init_nb=init_nb
                    )
                    na = na_tmp
                    nb = nb_tmp
                    @info "Saved state" state_path=state_path
                end
            end
        else
            @info "No cached state found; running DMRG" state_path=state_path
            energy, psi, sites, _ = run_dmrg_to_log(
                logcfg.log_path;
                checkpoint_params_path=params_path,
                dmrg_cfg...
            )
            if save_state_flag
                na_tmp, nb_tmp = measure_densities(psi, sites)
                save_state(
                    state_path,
                    psi;
                    energy=energy,
                    params_path=params_path,
                    na=na_tmp,
                    nb=nb_tmp,
                    init_na=init_na,
                    init_nb=init_nb
                )
                na = na_tmp
                nb = nb_tmp
                @info "Saved state" state_path=state_path
            end
        end

        if na === nothing || nb === nothing
            na, nb = measure_densities(psi, sites)
        end
        obs = compute_observables(
            psi,
            sites;
            energy=energy,
            na=na,
            nb=nb,
            cfg=cfg,
            compute_density_density=true,
            compute_structure_factor=true,
            compute_triple_corr=false
        )

        mode = isfile(results_path) ? "r+" : "w"
        HDF5.h5open(results_path, mode) do f
            g_meta = ensure_group(f, "meta")
            write_meta!(g_meta;
                params_path=params_path,
                state_path=state_path,
                state_params_sha256=current_hash
            )
            write_results_schema!(g_meta)
            if observables_loaded
                write_or_replace(g_meta, "observables_path", abspath(observables_path))
                write_or_replace(g_meta, "observables_sha256", bytes2hex(SHA.sha256(read(observables_path, String))))
            end
            write_observables_hdf5!(f, obs)
        end
        @info "Wrote results" results_path=results_path
    end

    close(logcfg.logio)

    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
