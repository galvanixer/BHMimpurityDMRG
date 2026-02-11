# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

include(joinpath(@__DIR__, "..", "src", "BHMimpurityDMRG.jl"))
using .BHMimpurityDMRG
using HDF5
using SHA
using Logging
using Dates
function main()
    params_path = length(ARGS) >= 1 ? ARGS[1] :
                  get(ENV, "PARAMS", joinpath(@__DIR__, "..", "configs", "parameters.yaml"))

    # Load parameters with fallback to empty dict if file not found or loading fails
    cfg = isfile(params_path) ? load_params(params_path) : Dict{String,Any}()
    dmrg_cfg = merge_sections(cfg, ["lattice", "local_hilbert", "initial_state", "hamiltonian", "dmrg"])
    io_cfg = get(cfg, "io", Dict{String,Any}())
    state_path = get(io_cfg, "state_save_path", "dmrg_state.h5")
    save_state_flag = get(io_cfg, "save_state", true)
    
    logcfg = setup_logger(io_cfg; default_log_path="run.log")
    logger = logcfg.logger

    with_logger(logger) do
        logstarttime = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
        @info "--------Logging starts here -------- ($logstarttime)--------"
        @info "Starting triple_corr" params_path = params_path state_path = state_path
        params_text = isfile(params_path) ? read(params_path, String) : nothing
        current_hash = params_text === nothing ? nothing : bytes2hex(SHA.sha256(params_text))
        _, init_na, init_nb = dmrg_initial_configuration(; dmrg_cfg...)
        t_total = time_ns()
        t_state = time_ns()
        na = nothing
        nb = nothing
        if isfile(state_path)
            st = load_state(state_path)
            if current_hash !== nothing && st.params_sha256 !== nothing && st.params_sha256 == current_hash
                psi = st.psi
                sites = st.sites
                energy = st.energy
                na = st.na
                nb = st.nb
                @info "Loaded cached state" state_path = state_path
            else
                @info "State hash mismatch; recomputing"
                energy, psi, sites, H = run_dmrg(; checkpoint_params_path=params_path, dmrg_cfg...)
                if save_state_flag
                    na, nb = measure_densities(psi, sites)
                    save_state(
                        state_path,
                        psi;
                        energy=energy,
                        params_path=params_path,
                        na=na,
                        nb=nb,
                        init_na=init_na,
                        init_nb=init_nb
                    )
                    @info "Saved state" state_path = state_path
                end
            end
        else
            @info "No cached state; running DMRG"
            energy, psi, sites, H = run_dmrg(; checkpoint_params_path=params_path, dmrg_cfg...)
            if save_state_flag
                na, nb = measure_densities(psi, sites)
                save_state(
                    state_path,
                    psi;
                    energy=energy,
                    params_path=params_path,
                    na=na,
                    nb=nb,
                    init_na=init_na,
                    init_nb=init_nb
                )
                @info "Saved state" state_path = state_path
            end
        end
        @info "State ready" seconds = (time_ns() - t_state) / 1e9
        flush(logcfg.logio) # Ensure all logs are written before proceeding

        if na === nothing || nb === nothing
            na, nb = measure_densities(psi, sites)
        end
        Na, Nb = total_numbers(na, nb)

        # Measure translationally averaged connected 3-point correlator:
        # C^(3)(r, s) = (1/L) Σ_i ⟨n_i n_{i+r} n_{i+s}⟩_c
        #
        # Since `run_dmrg` currently builds a periodic Hamiltonian, use periodic=true here.
        obs_cfg = get(cfg, "observables", Dict{String,Any}())
        tc_cfg = get(obs_cfg, "triple_corr", Dict{String,Any}())
        results_path = get(io_cfg, "results_path", "results.h5")
        periodic = get(get(cfg, "lattice", Dict{String,Any}()), "periodic", true)
        species = lowercase(String(get(tc_cfg, "species", "both")))
        opnames = species == "a" ? ["Na"] :
                  species == "b" ? ["Nb"] :
                  species == "both" ? ["Na", "Nb"] :
                  error("observables.triple_corr.species must be one of: \"a\", \"b\", \"both\"")
        precompute = get(tc_cfg, "precompute", true)

        nvec_a = nothing
        nnmat_a = nothing
        nvec_b = nothing
        nnmat_b = nothing
        if precompute
            t_pre = time_ns()
            if "Na" in opnames
                nvec_a = precompute_n(psi, sites, "Na")
                nnmat_a = precompute_nn(psi, sites, "Na")
            end
            if "Nb" in opnames
                nvec_b = precompute_n(psi, sites, "Nb")
                nnmat_b = precompute_nn(psi, sites, "Nb")
            end
            @info "Precompute done" seconds = (time_ns() - t_pre) / 1e9
        end

        # Example separations (r, s). Add more as needed.
        pairs_spec = get(tc_cfg, "pairs", nothing)
        all_pairs = get(tc_cfg, "all_pairs", false) ||
                    (pairs_spec isa AbstractString && lowercase(pairs_spec) == "all")
        rmax = get(tc_cfg, "rmax", nothing)
        smax = get(tc_cfg, "smax", nothing)

        pairs = Vector{Tuple{Int,Int}}()
        if all_pairs || rmax !== nothing || smax !== nothing
            L = length(sites)
            rmax_eff = rmax === nothing ? (L - 1) : min(Int(rmax), L - 1)
            smax_eff = smax === nothing ? (L - 1) : min(Int(smax), L - 1)
            for r in 0:rmax_eff, s in 0:smax_eff
                push!(pairs, (r, s))
            end
        else
            raw_pairs = pairs_spec === nothing ?
                        [[0, 0], [0, 1], [0, 2], [1, 1], [1, 2], [1, 3], [2, 3]] :
                        pairs_spec
            pairs = [(Int(p[1]), Int(p[2])) for p in raw_pairs]
        end

        npairs = length(pairs)
        anchors = zeros(Int, npairs)
        results_a = ("Na" in opnames) ? zeros(Float64, npairs) : nothing
        results_b = ("Nb" in opnames) ? zeros(Float64, npairs) : nothing

        t_corr = time_ns()
        for (idx, (r, s)) in enumerate(pairs)
            if results_a !== nothing
                if precompute && nvec_a !== nothing && nnmat_a !== nothing
                    C, N = transl_avg_connected_nnn_no_cached(psi, sites, "Na", r, s;
                        periodic=periodic, nvec=nvec_a, nnmat=nnmat_a)
                else
                    C, N = transl_avg_connected_nnn_no(psi, sites, "Na", r, s; periodic=periodic)
                end
                results_a[idx] = C
                anchors[idx] = N
            end
            if results_b !== nothing
                if precompute && nvec_b !== nothing && nnmat_b !== nothing
                    C, N = transl_avg_connected_nnn_no_cached(psi, sites, "Nb", r, s;
                        periodic=periodic, nvec=nvec_b, nnmat=nnmat_b)
                else
                    C, N = transl_avg_connected_nnn_no(psi, sites, "Nb", r, s; periodic=periodic)
                end
                results_b[idx] = C
                anchors[idx] = N
            end
        end
        @info "Correlator loop done" pairs = npairs seconds = (time_ns() - t_corr) / 1e9

        # Write results to HDF5
        t_write = time_ns()
        mode = isfile(results_path) ? "r+" : "w"
        HDF5.h5open(results_path, mode) do f
            g_meta = ensure_group(f, "meta")
            write_meta!(g_meta;
                params_path=params_path,
                state_path=state_path,
                state_params_sha256=current_hash
            )
            g_obs = ensure_group(f, "observables")
            g_energy = ensure_group(g_obs, "energy")
            g_den = ensure_group(g_obs, "densities")
            g_tot = ensure_group(g_obs, "totals")
            g_tc = ensure_group(g_obs, "triple_corr")
            write_or_replace(g_energy, "E0", energy)
            write_or_replace(g_den, "na", na)
            write_or_replace(g_den, "nb", nb)
            write_or_replace(g_tot, "Na", Na)
            write_or_replace(g_tot, "Nb", Nb)
            write_or_replace(g_tc, "pairs", hcat([p[1] for p in pairs], [p[2] for p in pairs]))
            write_or_replace(g_tc, "anchors", anchors)
            if results_a !== nothing
                write_or_replace(g_tc, "C_no_a", results_a)
            end
            if results_b !== nothing
                write_or_replace(g_tc, "C_no_b", results_b)
            end
            if energy !== nothing
                write_or_replace(g_tc, "energy", energy)
            end
        end
        @info "Wrote results" results_path = results_path seconds = (time_ns() - t_write) / 1e9
        @info "Total time" seconds = (time_ns() - t_total) / 1e9
    end

    close(logcfg.logio)

    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
