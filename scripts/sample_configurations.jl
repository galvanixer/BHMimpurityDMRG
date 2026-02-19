# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

include(joinpath(@__DIR__, "..", "src", "BHMimpurityDMRG.jl"))
using .BHMimpurityDMRG
using HDF5
using Logging
using SHA

function print_help(io::IO=stdout)
    script = basename(@__FILE__)
    println(io, "Usage:")
    println(io, "  julia --project=. scripts/$script [params_path] [observables_path] [state_path] [results_path]")
    println(io, "")
    println(io, "Positional arguments:")
    println(io, "  params_path       Optional model/DMRG YAML path. Defaults to ./parameters.yaml")
    println(io, "  observables_path  Optional observables YAML path. Defaults to ./observables.yaml")
    println(io, "  state_path        Optional state HDF5 path. Defaults to io.state_save_path or dmrg_state.h5")
    println(io, "  results_path      Optional results HDF5 path. Defaults to io.results_path or results.h5")
    println(io, "")
    println(io, "Config section (in observables YAML):")
    println(io, "  observables.sampled_configs.nsamples")
    println(io, "  observables.sampled_configs.top_k")
    println(io, "  observables.sampled_configs.seed")
    println(io, "  observables.sampled_configs.write_decoded_occupations")
    println(io, "")
    println(io, "Examples:")
    println(io, "  julia --project=. scripts/$script")
    println(io, "  ./bin/sample_configurations")
    println(io, "  julia --project=. scripts/$script parameters.yaml observables.yaml dmrg_state.h5 results.h5")
end

function parse_bool(x, default::Bool=false)
    x === nothing && return default
    x isa Bool && return x
    s = lowercase(strip(String(x)))
    if s in ("1", "true", "yes", "y", "on")
        return true
    elseif s in ("0", "false", "no", "n", "off")
        return false
    end
    return default
end

function parse_optional_int(x)
    x === nothing && return nothing
    x isa String && isempty(strip(x)) && return nothing
    return Int(x)
end

function main()
    if any(a -> a in ("-h", "--help"), ARGS)
        print_help()
        return nothing
    end

    if length(ARGS) > 4
        print_help(stderr)
        error("Expected at most 4 positional arguments, got $(length(ARGS)).")
    end

    params_path = length(ARGS) >= 1 ? ARGS[1] : "parameters.yaml"
    observables_path_arg = length(ARGS) >= 2 ? ARGS[2] : "observables.yaml"

    cfg_base = isfile(params_path) ? load_params(params_path) : Dict{String,Any}()
    merged_cfg = with_observables_config(cfg_base; observables_path=observables_path_arg)
    cfg = merged_cfg.cfg
    observables_path = merged_cfg.observables_path
    observables_loaded = merged_cfg.observables_loaded
    params_has_observables = haskey(cfg_base, "observables") || haskey(cfg_base, :observables)

    io_cfg = get(cfg_base, "io", Dict{String,Any}())
    obs_cfg = get(cfg, "observables", Dict{String,Any}())
    samp_cfg = get(obs_cfg, "sampled_configs", Dict{String,Any}())

    state_path = length(ARGS) >= 3 ? ARGS[3] : get(io_cfg, "state_save_path", "dmrg_state.h5")
    results_path = length(ARGS) >= 4 ? ARGS[4] : get(io_cfg, "results_path", "results.h5")

    nsamples = Int(get(samp_cfg, "nsamples", 100_000))
    top_k = Int(get(samp_cfg, "top_k", 20))
    seed = parse_optional_int(get(samp_cfg, "seed", nothing))
    write_decoded = parse_bool(get(samp_cfg, "write_decoded_occupations", true), true)

    logcfg = setup_logger(io_cfg; default_log_path="sampled_configs.log")
    logger = logcfg.logger

    with_logger(logger) do
        @info "Sampling top configurations from saved DMRG state" params_path=params_path observables_path=observables_path observables_loaded=observables_loaded state_path=state_path results_path=results_path nsamples=nsamples top_k=top_k seed=seed
        if !observables_loaded && params_has_observables
            @info "Observables file not found; falling back to observables in parameters config" observables_path=observables_path
        elseif !observables_loaded
            @warn "Observables file not found and no observables section in parameters config; defaults will be used" observables_path=observables_path
        end

        isfile(state_path) || error("State file not found: $state_path. Run scripts/solve_ground_state.jl first.")
        st = load_state(state_path)

        sampled = top_sampled_configurations(st.psi; nsamples=nsamples, top_k=top_k, seed=seed)

        na_top = nothing
        nb_top = nothing
        nmax_b = nothing
        if write_decoded
            try
                nmax_b = infer_nmax_b_from_sites(st.sites)
                na_top, nb_top = decode_two_boson_configs(sampled.configs, nmax_b)
            catch err
                @warn "Could not decode basis indices to (na, nb); writing basis-index output only" error=string(err)
            end
        end

        params_text = isfile(params_path) ? read(params_path, String) : nothing
        current_hash = params_text === nothing ? nothing : bytes2hex(SHA.sha256(params_text))

        mode = isfile(results_path) ? "r+" : "w"
        HDF5.h5open(results_path, mode) do f
            g_meta = ensure_group(f, "meta")
            write_meta!(g_meta;
                params_path=params_path,
                state_path=state_path,
                state_params_sha256=current_hash
            )
            if observables_loaded
                write_or_replace(g_meta, "observables_path", abspath(observables_path))
                write_or_replace(g_meta, "observables_sha256", bytes2hex(SHA.sha256(read(observables_path, String))))
            end

            g_obs = ensure_group(f, "observables")
            g_samp = ensure_group(g_obs, "sampled_configs")
            write_or_replace(g_samp, "basis_indices_top", sampled.configs)
            write_or_replace(g_samp, "counts_top", sampled.counts)
            write_or_replace(g_samp, "probabilities_top", sampled.probs)
            write_or_replace(g_samp, "stderr_top", sampled.stderrs)
            write_or_replace(g_samp, "nsamples", Int(sampled.nsamples))
            write_or_replace(g_samp, "unique_configurations", Int(sampled.unique_configurations))
            write_or_replace(g_samp, "top_k_requested", Int(top_k))
            if seed !== nothing
                write_or_replace(g_samp, "seed", Int(seed))
            end
            if na_top !== nothing && nb_top !== nothing && nmax_b !== nothing
                write_or_replace(g_samp, "nmax_b", Int(nmax_b))
                write_or_replace(g_samp, "na_top", na_top)
                write_or_replace(g_samp, "nb_top", nb_top)
            end
        end

        nshow = min(size(sampled.configs, 1), 5)
        for r in 1:nshow
            if na_top === nothing || nb_top === nothing
                @info "Top sampled config (basis indices)" rank=r probability=sampled.probs[r] stderr=sampled.stderrs[r] count=sampled.counts[r] config=collect(sampled.configs[r, :])
            else
                occ = [(na_top[r, c], nb_top[r, c]) for c in 1:size(na_top, 2)]
                @info "Top sampled config (decoded)" rank=r probability=sampled.probs[r] stderr=sampled.stderrs[r] count=sampled.counts[r] occupations=occ
            end
        end

        @info "Wrote sampled-configuration results" results_path=results_path top_kept=size(sampled.configs, 1) unique_configurations=sampled.unique_configurations
    end

    close(logcfg.logio)
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
