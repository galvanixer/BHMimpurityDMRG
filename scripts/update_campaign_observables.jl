# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

include(joinpath(@__DIR__, "..", "src", "BHMimpurityDMRG.jl"))
using .BHMimpurityDMRG
using HDF5
using YAML
using SHA
using ITensors
using ITensorMPS
using Dates

const VALID_MISSING_RESULTS_POLICIES = Set(["skip", "from-checkpoint"])
const CHECKPOINT_RESULTS_FILENAME = "results_from_checkpoint.h5"

function print_help(io::IO=stdout)
    script = basename(@__FILE__)
    println(io, "Usage:")
    println(io, "  julia --project=. scripts/$script [options] <campaign_root> <observables_yaml>")
    println(io, "")
    println(io, "Options:")
    println(io, "  --missing-results <skip|from-checkpoint>  What to do when configured results file is missing (default: skip)")
    println(io, "  --quiet                                    Reduce progress output")
    println(io, "  -h, --help                                 Show this help")
    println(io, "")
    println(io, "Behavior:")
    println(io, "  - Updates /observables in each run's configured results file (io.results_path; default: results.h5).")
    println(io, "  - If the configured results file is missing:")
    println(io, "      skip            -> do nothing for that run")
    println(io, "      from-checkpoint -> generate $CHECKPOINT_RESULTS_FILENAME from dmrg_state_checkpoint.h5 with same observables YAML")
    println(io, "")
    println(io, "Examples:")
    println(io, "  ./bin/update_campaign_observables runs/my_campaign configs/observables.yaml")
    println(io, "  ./bin/update_campaign_observables --missing-results from-checkpoint runs/my_campaign configs/observables.yaml")
end

@inline function parse_bool(x, default::Bool=false)
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

function normalize_yaml(x)
    if x isa AbstractDict
        return Dict{String,Any}(String(k) => normalize_yaml(v) for (k, v) in x)
    elseif x isa AbstractVector
        return [normalize_yaml(v) for v in x]
    else
        return x
    end
end

@inline function clean_error(err)
    return replace(sprint(showerror, err), '\n' => ' ')
end

@inline function resolved_path(base_dir::AbstractString, p)
    p === nothing && return nothing
    s = String(p)
    return isabspath(s) ? s : abspath(joinpath(base_dir, s))
end

function parse_args(args::Vector{String})
    show_help = false
    verbose = true
    missing_results_policy = "skip"

    positional = String[]
    i = 1
    while i <= length(args)
        a = args[i]
        if a in ("-h", "--help")
            show_help = true
            i += 1
        elseif a == "--quiet"
            verbose = false
            i += 1
        elseif a == "--missing-results"
            i < length(args) || error("--missing-results requires a value")
            missing_results_policy = lowercase(strip(args[i + 1]))
            i += 2
        elseif startswith(a, "--")
            error("Unknown option: $a")
        else
            push!(positional, a)
            i += 1
        end
    end

    missing_results_policy in VALID_MISSING_RESULTS_POLICIES ||
        error("Invalid --missing-results: $missing_results_policy (expected one of: $(join(sort(collect(VALID_MISSING_RESULTS_POLICIES)), ", ")))")

    campaign_root = nothing
    observables_path = nothing
    if length(positional) >= 1
        campaign_root = positional[1]
    end
    if length(positional) >= 2
        observables_path = positional[2]
    end
    if length(positional) > 2
        error("Expected 2 positional arguments, got $(length(positional)).")
    end

    return (
        show_help=show_help,
        verbose=verbose,
        missing_results_policy=missing_results_policy,
        campaign_root=campaign_root,
        observables_path=observables_path
    )
end

@inline function looks_like_run_dir_name(name::AbstractString)
    return occursin(r"^run_[0-9]+$", name)
end

function parse_runs_csv(path::AbstractString)
    lines = readlines(path)
    isempty(lines) && return NamedTuple[]
    header = split(lines[1], ",")
    idx_run_id = findfirst(==("run_id"), header)
    idx_run_dir = findfirst(==("run_dir"), header)
    idx_params = findfirst(==("params_path"), header)

    (idx_run_id === nothing || idx_run_dir === nothing) && return NamedTuple[]

    rows = NamedTuple[]
    for ln in lines[2:end]
        s = strip(ln)
        isempty(s) && continue
        cols = split(s, ",")
        maxidx = max(idx_run_id, idx_run_dir, idx_params === nothing ? 0 : idx_params)
        length(cols) < maxidx && continue
        run_id = strip(cols[idx_run_id])
        run_dir = strip(cols[idx_run_dir])
        params_path = idx_params === nothing ? "" : strip(cols[idx_params])
        isempty(run_id) && continue
        isempty(run_dir) && continue
        push!(rows, (
            run_id=run_id,
            run_dir=run_dir,
            params_path=params_path
        ))
    end
    return rows
end

function discover_runs(campaign_root::AbstractString)
    campaign_root_abs = abspath(campaign_root)
    isdir(campaign_root_abs) || error("campaign_root is not a directory: $campaign_root_abs")

    runs_csv = joinpath(campaign_root_abs, "runs.csv")
    if isfile(runs_csv)
        parsed = parse_runs_csv(runs_csv)
        if !isempty(parsed)
            rows = NamedTuple[]
            for r in parsed
                run_dir_abs = resolved_path(campaign_root_abs, r.run_dir)
                params_abs = isempty(r.params_path) ? joinpath(run_dir_abs, "parameters.yaml") : resolved_path(campaign_root_abs, r.params_path)
                push!(rows, (run_id=r.run_id, run_dir=run_dir_abs, params_path=params_abs))
            end
            sort!(rows, by=x -> x.run_id)
            return rows
        end
    end

    rows = NamedTuple[]
    for ent in sort(readdir(campaign_root_abs))
        looks_like_run_dir_name(ent) || continue
        run_dir = joinpath(campaign_root_abs, ent)
        isdir(run_dir) || continue
        push!(rows, (run_id=ent, run_dir=run_dir, params_path=joinpath(run_dir, "parameters.yaml")))
    end
    return rows
end

function compute_energy_stats_from_state(st, cfg::AbstractDict)
    need_hamiltonian = st.energy === nothing || st.energy_variance === nothing
    H = need_hamiltonian ? build_hamiltonian_from_config(st.sites, cfg) : nothing
    energy = st.energy !== nothing ? Float64(real(st.energy)) : expect_operator(st.psi, H)
    energy_variance = st.energy_variance !== nothing ?
                      Float64(real(st.energy_variance)) :
                      operator_variance(st.psi, H; expectation=energy)
    return energy, energy_variance
end

function compute_initial_densities_from_state(st, cfg::AbstractDict)
    if st.init_na !== nothing && st.init_nb !== nothing
        return st.init_na, st.init_nb
    end
    init_cfg = merge_sections(cfg, ["lattice", "local_hilbert", "initial_state"])
    _, init_na, init_nb = dmrg_initial_configuration(; init_cfg...)
    return init_na, init_nb
end

function build_observables_from_state(st, cfg::AbstractDict)
    obs_cfg = get(cfg, "observables", Dict{String,Any}())
    lattice_cfg = get(cfg, "lattice", Dict{String,Any}())
    periodic = parse_bool(get(lattice_cfg, "periodic", true), true)

    psi = st.psi
    sites = st.sites
    energy, energy_variance = compute_energy_stats_from_state(st, cfg)
    init_na, init_nb = compute_initial_densities_from_state(st, cfg)

    na, nb = if st.na !== nothing && st.nb !== nothing
        st.na, st.nb
    else
        measure_densities(psi, sites)
    end

    dd_requested = haskey(obs_cfg, "density_density")
    sf_requested = haskey(obs_cfg, "structure_factor")
    pdd_requested = haskey(obs_cfg, "pair_distance_distribution")
    spdm_requested = haskey(obs_cfg, "single_particle_density_matrix")
    tc_requested = haskey(obs_cfg, "triple_corr")

    obs = compute_observables(
        psi,
        sites;
        energy=energy,
        energy_variance=energy_variance,
        init_na=init_na,
        init_nb=init_nb,
        na=na,
        nb=nb,
        cfg=cfg,
        periodic=periodic,
        compute_density_density=dd_requested,
        compute_structure_factor=sf_requested,
        compute_pair_distance_distribution=pdd_requested,
        compute_single_particle_density_matrix=spdm_requested,
        compute_triple_corr=tc_requested,
        progress=false
    )
    return obs
end

function update_results_observables!(
    results_path::AbstractString,
    state_source_path::AbstractString,
    params_path::AbstractString,
    observables_path::AbstractString,
    cfg::AbstractDict;
    mode::AbstractString
)
    st = load_state(state_source_path)
    obs = build_observables_from_state(st, cfg)

    params_text = isfile(params_path) ? read(params_path, String) :
                  (st.params_yaml === nothing ? nothing : String(st.params_yaml))
    current_hash = params_text === nothing ? nothing : bytes2hex(SHA.sha256(params_text))

    HDF5.h5open(results_path, mode) do f
        g_meta = ensure_group(f, "meta")
        write_meta!(
            g_meta;
            params_path=isfile(params_path) ? params_path : nothing,
            params_text=params_text,
            state_path=state_source_path,
            state_params_sha256=current_hash
        )
        write_results_schema!(g_meta)
        write_or_replace(g_meta, "observables_path", abspath(observables_path))
        write_or_replace(g_meta, "observables_sha256", bytes2hex(SHA.sha256(read(observables_path, String))))
        if st.checkpoint_sweep !== nothing
            write_or_replace(g_meta, "checkpoint_sweep", Int(st.checkpoint_sweep))
        end

        write_observables_hdf5!(f, obs)
    end
    return nothing
end

function process_run(
    run_id::AbstractString,
    run_dir::AbstractString,
    params_path::AbstractString,
    observables_path::AbstractString,
    missing_results_policy::AbstractString
)
    cfg_base = isfile(params_path) ? load_params(params_path) : Dict{String,Any}()
    merged_cfg = with_observables_config(cfg_base; observables_path=observables_path)
    cfg = merged_cfg.cfg

    io_cfg = get(cfg_base, "io", Dict{String,Any}())
    dmrg_cfg = get(cfg_base, "dmrg", Dict{String,Any}())

    results_path = resolved_path(run_dir, get(io_cfg, "results_path", "results.h5"))
    checkpoint_results_path = joinpath(dirname(results_path), CHECKPOINT_RESULTS_FILENAME)
    state_path = resolved_path(run_dir, get(io_cfg, "state_save_path", "dmrg_state.h5"))
    checkpoint_path = resolved_path(run_dir, get(dmrg_cfg, "checkpoint_path", "dmrg_state_checkpoint.h5"))

    results_exists = isfile(results_path)
    if results_exists
        if isfile(state_path)
            state_source = state_path
        elseif isfile(checkpoint_path)
            state_source = checkpoint_path
        else
            return (run_id=String(run_id), action="skipped", reason="no_state_or_checkpoint_for_existing_results", results_path=results_path)
        end

        update_results_observables!(results_path, state_source, params_path, observables_path, cfg; mode="r+")
        return (run_id=String(run_id), action="updated", reason="updated_existing_results", results_path=results_path)
    end

    if missing_results_policy == "skip"
        return (run_id=String(run_id), action="skipped", reason="results_missing_policy_skip", results_path=results_path)
    end

    # missing_results_policy == "from-checkpoint"
    if !isfile(checkpoint_path)
        return (run_id=String(run_id), action="skipped", reason="results_missing_no_checkpoint", results_path=checkpoint_results_path)
    end

    mkpath(dirname(checkpoint_results_path))
    update_results_observables!(checkpoint_results_path, checkpoint_path, params_path, observables_path, cfg; mode="w")
    return (run_id=String(run_id), action="generated", reason="generated_from_checkpoint", results_path=checkpoint_results_path)
end

function main(args=ARGS)
    opts = parse_args(args)
    if opts.show_help
        print_help()
        return nothing
    end
    opts.campaign_root === nothing && error("campaign_root is required. Use --help for usage.")
    opts.observables_path === nothing && error("observables_yaml is required. Use --help for usage.")

    campaign_root_abs = abspath(opts.campaign_root)
    observables_path_abs = abspath(opts.observables_path)
    isdir(campaign_root_abs) || error("campaign_root is not a directory: $campaign_root_abs")
    isfile(observables_path_abs) || error("observables_yaml file not found: $observables_path_abs")

    runs = discover_runs(campaign_root_abs)
    isempty(runs) && error("No run directories found under: $campaign_root_abs")

    n_updated = 0
    n_generated = 0
    n_skipped = 0
    n_failed = 0

    opts.verbose && println("Updating observables across campaign runs")
    opts.verbose && println("campaign_root      : $campaign_root_abs")
    opts.verbose && println("observables_yaml   : $observables_path_abs")
    opts.verbose && println("missing_results    : $(opts.missing_results_policy)")
    opts.verbose && println("runs_discovered    : $(length(runs))")

    for (idx, r) in enumerate(runs)
        opts.verbose && println("[$idx/$(length(runs))] $(r.run_id)")
        try
            out = process_run(r.run_id, r.run_dir, r.params_path, observables_path_abs, opts.missing_results_policy)
            if out.action == "updated"
                n_updated += 1
            elseif out.action == "generated"
                n_generated += 1
            elseif out.action == "skipped"
                n_skipped += 1
            end
            opts.verbose && println("  -> $(out.action): $(out.reason)")
        catch err
            n_failed += 1
            println(stderr, "  -> failed: run=$(r.run_id) error=$(clean_error(err))")
        end
    end

    println("Finished observables update at $(Dates.format(now(), dateformat"yyyy-mm-dd HH:MM:SS"))")
    println("Updated existing results : $n_updated")
    println("Generated missing results: $n_generated")
    println("Skipped                 : $n_skipped")
    println("Failed                  : $n_failed")
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
