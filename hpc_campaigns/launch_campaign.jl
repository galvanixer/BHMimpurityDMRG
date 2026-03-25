#!/usr/bin/env julia
# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

using Printf
using Dates
using Random
using OrderedCollections: OrderedDict
include(joinpath(@__DIR__, "yaml_helpers.jl"))
using .YAMLHelpers: as_dict
using .YAMLHelpers: load_ordered_yaml
using .YAMLHelpers: set_nested!
using .YAMLHelpers: sweep_assignments
using .YAMLHelpers: write_yaml_ordered

function shell_quote(s::AbstractString)
    return "'" * replace(String(s), "'" => "'\"'\"'") * "'"
end

function apply_overrides!(cfg::AbstractDict, overrides::AbstractDict)
    for (k, v) in overrides
        set_nested!(cfg, String(k), v)
    end
    return cfg
end

function as_sweep_vector(value)
    return value isa AbstractVector ? collect(value) : [value]
end

function linked_sweep_assignments(linked_sweep::AbstractDict)
    linked = OrderedDict{String,Any}(String(k) => v for (k, v) in linked_sweep)
    keys_sorted = sort(collect(keys(linked)))
    if isempty(keys_sorted)
        return [OrderedDict{String,Any}()]
    end

    values = [as_sweep_vector(linked[k]) for k in keys_sorted]
    lengths = [length(v) for v in values]
    n = first(lengths)
    n > 0 || error("linked_sweep values must be non-empty")
    all(==(n), lengths) || error("All linked_sweep arrays must have the same length, got lengths: " * join(lengths, ", "))

    assigns = OrderedDict{String,Any}[]
    for i in 1:n
        a = OrderedDict{String,Any}()
        for (k, v) in zip(keys_sorted, values)
            a[k] = v[i]
        end
        push!(assigns, a)
    end
    return assigns
end

function combine_sweep_assignments(cart_assigns::AbstractVector, linked_assigns::AbstractVector)
    out = OrderedDict{String,Any}[]
    for cart in cart_assigns, linked in linked_assigns
        merged = OrderedDict{String,Any}()
        for source in (cart, linked)
            for (k, v) in source
                key = String(k)
                if haskey(merged, key) && merged[key] != v
                    error("Conflicting assignments for key '$key' between sweep and linked_sweep")
                end
                merged[key] = v
            end
        end
        push!(out, merged)
    end
    return out
end

function parse_int_setting(value, key::AbstractString)
    if value isa Integer
        return Int(value)
    elseif value isa AbstractFloat
        isfinite(value) || error("$key must be finite, got: $value")
        isinteger(value) || error("$key must be an integer value, got: $value")
        return Int(round(value))
    elseif value isa AbstractString
        s = strip(String(value))
        isempty(s) && error("$key must be non-empty")
        try
            return parse(Int, s)
        catch
            error("$key must be an integer, got: $value")
        end
    else
        error("$key must be an integer, got $(typeof(value))")
    end
end

function parse_positive_int_setting(value, key::AbstractString)
    v = parse_int_setting(value, key)
    v > 0 || error("$key must be a positive integer, got: $v")
    return v
end

function resolve_seed_count(seed_cfg::AbstractDict)
    count_raw = get(seed_cfg, "count", "prompt")
    if count_raw isa AbstractString && lowercase(strip(String(count_raw))) == "prompt"
        isatty(stdin) || error("seed_generation.count=prompt requires interactive stdin; set seed_generation.count to an integer for non-interactive runs.")
        while true
            print("Enter number of seeds to generate: ")
            flush(stdout)
            line = strip(readline(stdin))
            isempty(line) && continue
            try
                count = parse(Int, line)
                if count > 0
                    return count
                end
            catch
            end
            println("Please enter a positive integer.")
        end
    end
    return parse_positive_int_setting(count_raw, "seed_generation.count")
end

function generate_auto_seeds(seed_cfg::AbstractDict)
    mode_raw = get(seed_cfg, "mode", "random")
    mode = lowercase(strip(String(mode_raw)))
    count = resolve_seed_count(seed_cfg)

    if mode == "sequential"
        start = parse_positive_int_setting(get(seed_cfg, "start", 14_061_990), "seed_generation.start")
        step = parse_positive_int_setting(get(seed_cfg, "step", 1), "seed_generation.step")
        max_seed = start + (count - 1) * step
        max_seed >= start || error("seed_generation produced overflow; reduce count/start/step.")
        seeds = [start + (i - 1) * step for i in 1:count]
        return seeds, mode
    elseif mode == "random"
        min_seed = parse_positive_int_setting(get(seed_cfg, "min", 1), "seed_generation.min")
        max_seed = parse_positive_int_setting(get(seed_cfg, "max", 2_147_483_647), "seed_generation.max")
        min_seed <= max_seed || error("seed_generation.min must be <= seed_generation.max")
        range_size = max_seed - min_seed + 1
        count <= range_size || error("seed_generation.count=$count exceeds available unique random seeds in range [$min_seed, $max_seed]")

        rng = if haskey(seed_cfg, "generator_seed")
            Random.MersenneTwister(parse_int_setting(seed_cfg["generator_seed"], "seed_generation.generator_seed"))
        else
            Random.default_rng()
        end

        seeds = Int[]
        seen = Set{Int}()
        while length(seeds) < count
            candidate = rand(rng, min_seed:max_seed)
            if !(candidate in seen)
                push!(seen, candidate)
                push!(seeds, candidate)
            end
        end
        return seeds, mode
    else
        error("Unsupported seed_generation.mode '$mode_raw'. Expected 'sequential' or 'random'.")
    end
end

function maybe_expand_auto_seeds!(root_cfg::AbstractDict, sweep::AbstractDict)
    haskey(sweep, "initial_state.seed") || return
    seed_spec = sweep["initial_state.seed"]
    if !(seed_spec isa AbstractString && uppercase(strip(String(seed_spec))) == "AUTO")
        return
    end

    seed_cfg = haskey(root_cfg, "seed_generation") ? as_dict(root_cfg["seed_generation"]) : OrderedDict{String,Any}()
    seeds, mode = generate_auto_seeds(seed_cfg)
    sweep["initial_state.seed"] = seeds

    n_show = min(length(seeds), 10)
    preview = join(seeds[1:n_show], ", ")
    if length(seeds) > n_show
        preview *= ", ..."
    end

    println("Auto seed generation:")
    println("  mode:  ", mode)
    println("  count: ", length(seeds))
    println("  seeds: ", preview)
end

@inline function dict_get(d::AbstractDict, key::String, default=nothing)
    if haskey(d, key)
        return d[key]
    end
    ks = Symbol(key)
    if haskey(d, ks)
        return d[ks]
    end
    return default
end

@inline function nested_get(d::AbstractDict, path::AbstractVector{<:AbstractString}, default=nothing)
    cur = d
    for p in path
        if cur isa AbstractDict && haskey(cur, p)
            cur = cur[p]
        else
            return default
        end
    end
    return cur
end

@inline function dotted_or_nested_get(d::AbstractDict, dotted_key::String, default=nothing)
    if haskey(d, dotted_key)
        return d[dotted_key]
    end
    return nested_get(d, split(dotted_key, "."), default)
end

@inline function normalize_impurity_distribution(value)
    value === nothing && return nothing
    s = lowercase(strip(string(value)))
    isempty(s) && return nothing
    return Symbol(s)
end

@inline function impurity_distribution_uses_seed(impdist)
    impdist === nothing && return false
    return impdist in (:random_capped, :random_separated)
end

function effective_impurity_distribution(base_cfg::AbstractDict, layers::AbstractDict...)
    value = nested_get(base_cfg, ["initial_state", "impurity_distribution"], nothing)
    for layer in layers
        candidate = dotted_or_nested_get(layer, "initial_state.impurity_distribution", nothing)
        if candidate !== nothing
            value = candidate
        end
    end
    impdist = normalize_impurity_distribution(value)
    return impdist === nothing ? :centered_pileup : impdist
end

function all_impurity_distributions_seed_independent(
    base_cfg::AbstractDict,
    overrides::AbstractDict,
    sweep::AbstractDict,
    linked_sweep::AbstractDict
)
    found_explicit_distribution = false
    for source in (sweep, linked_sweep)
        if haskey(source, "initial_state.impurity_distribution")
            found_explicit_distribution = true
            for value in as_sweep_vector(source["initial_state.impurity_distribution"])
                impurity_distribution_uses_seed(normalize_impurity_distribution(value)) && return false
            end
        end
    end

    if found_explicit_distribution
        return true
    end
    return !impurity_distribution_uses_seed(effective_impurity_distribution(base_cfg, overrides))
end

function prune_irrelevant_seed_sweeps!(base_cfg::AbstractDict, overrides::AbstractDict, sweep::AbstractDict, linked_sweep::AbstractDict)
    has_seed_axis = haskey(sweep, "initial_state.seed") || haskey(linked_sweep, "initial_state.seed")
    has_seed_axis || return false

    all_impurity_distributions_seed_independent(base_cfg, overrides, sweep, linked_sweep) || return false

    if haskey(sweep, "initial_state.seed")
        delete!(sweep, "initial_state.seed")
    end
    if haskey(linked_sweep, "initial_state.seed")
        delete!(linked_sweep, "initial_state.seed")
    end
    return true
end

function collapse_seed_independent_assignments(base_cfg::AbstractDict, overrides::AbstractDict, assigns::AbstractVector)
    out = OrderedDict{String,Any}[]
    seen_seed_independent = Set{String}()
    dropped = 0

    for assign in assigns
        assign_norm = copy(assign)
        if !impurity_distribution_uses_seed(effective_impurity_distribution(base_cfg, overrides, assign_norm)) &&
           haskey(assign_norm, "initial_state.seed")
            delete!(assign_norm, "initial_state.seed")
            sig = repr(assign_norm)
            if sig in seen_seed_independent
                dropped += 1
                continue
            end
            push!(seen_seed_independent, sig)
        end
        push!(out, assign_norm)
    end

    return out, dropped
end

function maybe_expand_auto_meta_date!(cfg::AbstractDict; timestamp::AbstractString)
    meta = dict_get(cfg, "meta", nothing)
    meta isa AbstractDict || return false

    date_value = dict_get(meta, "date", nothing)
    if date_value isa AbstractString && uppercase(strip(String(date_value))) == "AUTO"
        set_nested!(cfg, "meta.date", String(timestamp))
        return true
    end
    return false
end

function ddmmmyyyy_tag(dt::Dates.TimeType=now())
    month_abbr = ("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")
    d = Dates.Date(dt)
    return lpad(string(Dates.day(d)), 2, '0') * month_abbr[Dates.month(d)] * string(Dates.year(d))
end

function next_versioned_campaign_name(output_root::AbstractString, stem::AbstractString)
    prefix = stem * "_v"
    max_v = 0
    if isdir(output_root)
        for entry in readdir(output_root)
            startswith(entry, prefix) || continue
            length(entry) > length(prefix) || continue
            suffix = entry[(length(prefix) + 1):end]
            v = tryparse(Int, suffix)
            if v !== nothing
                max_v = max(max_v, v)
            end
        end
    end
    next_v = max_v + 1
    return "$(stem)_v$(next_v)"
end

function main()
    length(ARGS) >= 1 ||
        error("Usage: julia hpc_campaigns/launch_campaign.jl <campaign.yaml> [output_root_override]")
    campaign_file = abspath(ARGS[1])
    isfile(campaign_file) || error("Campaign file not found: $campaign_file")

    c = load_ordered_yaml(campaign_file)
    campaign = haskey(c, "campaign") ? as_dict(c["campaign"]) : OrderedDict{String,Any}()

    campaign_base_name = String(get(campaign, "name", "campaign"))
    output_root_raw = if length(ARGS) >= 2
        String(ARGS[2])
    elseif haskey(campaign, "output_root_abs")
        String(campaign["output_root_abs"])
    else
        String(get(campaign, "output_root", "runs"))
    end
    output_root = isabspath(output_root_raw) ? output_root_raw :
                  abspath(joinpath(dirname(campaign_file), output_root_raw))
    campaign_date_tag = ddmmmyyyy_tag(now())
    campaign_stem = "$(campaign_base_name)_$(campaign_date_tag)"
    campaign_name = next_versioned_campaign_name(output_root, campaign_stem)
    app_script = String(get(campaign, "app_script", "scripts/solve_ground_state.jl"))
    params_filename = get(campaign, "params_filename", "parameters.yaml")

    base_config_path = haskey(c, "base_config") ? String(c["base_config"]) :
                       error("Missing required key: base_config")
    repo_root = abspath(joinpath(@__DIR__, ".."))
    base_config_abs = if isabspath(base_config_path)
        base_config_path
    else
        by_campaign = abspath(joinpath(dirname(campaign_file), base_config_path))
        isfile(by_campaign) ? by_campaign : abspath(joinpath(repo_root, base_config_path))
    end
    isfile(base_config_abs) || error("Base config not found: $base_config_abs")
    base_cfg = load_ordered_yaml(base_config_abs)
    app_script_abs = if isabspath(app_script)
        app_script
    else
        by_campaign = abspath(joinpath(dirname(campaign_file), app_script))
        isfile(by_campaign) ? by_campaign : abspath(joinpath(repo_root, app_script))
    end
    isfile(app_script_abs) || error("App script not found: $app_script_abs")

    overrides = haskey(c, "overrides") ? as_dict(c["overrides"]) : OrderedDict{String,Any}()
    sweep = haskey(c, "sweep") ? as_dict(c["sweep"]) : OrderedDict{String,Any}()
    linked_sweep = haskey(c, "linked_sweep") ? as_dict(c["linked_sweep"]) : OrderedDict{String,Any}()
    ignored_seed_sweep_for_seed_independent = prune_irrelevant_seed_sweeps!(base_cfg, overrides, sweep, linked_sweep)
    maybe_expand_auto_seeds!(c, sweep)
    cart_assigns = sweep_assignments(sweep)
    linked_assigns = linked_sweep_assignments(linked_sweep)
    assigns = combine_sweep_assignments(cart_assigns, linked_assigns)
    assigns, dropped_seed_independent_duplicates = collapse_seed_independent_assignments(base_cfg, overrides, assigns)

    campaign_dir = joinpath(output_root, campaign_name)
    mkpath(campaign_dir)
    campaign_yaml_copy = joinpath(campaign_dir, basename(campaign_file))
    cp(campaign_file, campaign_yaml_copy; force=true)

    index_path = joinpath(campaign_dir, "runs.csv")
    jobfile_path = joinpath(campaign_dir, "jobfile")
    run_dirs_path = joinpath(campaign_dir, "run_dirs.txt")
    open(index_path, "w") do io
        println(io, "run_id,run_dir,params_path,app_script,status")
    end
    open(jobfile_path, "w") do _ end
    open(run_dirs_path, "w") do _ end

    campaign_timestamp = Dates.format(now(), dateformat"yyyy-mm-dd HH:MM:SS")
    auto_date_expanded_count = 0

    for (i, assign) in enumerate(assigns)
        run_id = @sprintf("run_%04d", i)
        run_dir = joinpath(campaign_dir, run_id)
        mkpath(run_dir)

        cfg = deepcopy(base_cfg)
        apply_overrides!(cfg, overrides)
        apply_overrides!(cfg, assign)
        auto_date_expanded_count += maybe_expand_auto_meta_date!(cfg; timestamp=campaign_timestamp) ? 1 : 0

        # Force output files to stay in each run directory.
        set_nested!(cfg, "io.state_save_path", abspath(joinpath(run_dir, "dmrg_state.h5")))
        set_nested!(cfg, "io.results_path", abspath(joinpath(run_dir, "results.h5")))
        set_nested!(cfg, "io.log_path", abspath(joinpath(run_dir, "run.log")))
        set_nested!(cfg, "dmrg.checkpoint_path", abspath(joinpath(run_dir, "dmrg_state_checkpoint.h5")))

        # Keep run id in metadata for traceability.
        set_nested!(cfg, "meta.run_name", run_id)

        params_path = abspath(joinpath(run_dir, params_filename))
        write_yaml_ordered(params_path, cfg; reference=base_cfg, inline_keys=["maxdim"])

        open(index_path, "a") do io
            println(io, join([run_id, run_dir, params_path, app_script_abs, "PENDING"], ","))
        end
        open(jobfile_path, "a") do io
            println(io, "julia --project=$(shell_quote(repo_root)) $(shell_quote(app_script_abs)) $(shell_quote(params_path))")
        end
        open(run_dirs_path, "a") do io
            println(io, run_dir)
        end
    end

    println("Campaign generated:")
    println("  base_name: ", campaign_base_name)
    println("  date_tag: ", campaign_date_tag)
    println("  name: ", campaign_name)
    println("  runs: ", length(assigns))
    println("  dir:  ", campaign_dir)
    println("  campaign_yaml: ", campaign_yaml_copy)
    println("  runs_csv: ", index_path)
    println("  jobfile: ", jobfile_path)
    if ignored_seed_sweep_for_seed_independent
        println("  note: initial_state.seed ignored because impurity_distribution is fixed to a deterministic mode")
    elseif dropped_seed_independent_duplicates > 0
        println("  note: removed $dropped_seed_independent_duplicates deterministic-mode seed duplicate run(s)")
    end
    if auto_date_expanded_count > 0
        println("  meta.date AUTO expanded in $auto_date_expanded_count run(s) as: ", campaign_timestamp)
    end
    println("Next step:")
    println("  bash hpc_campaigns/slurm/submit_multilauncher.sh ", campaign_dir)
end

main()
