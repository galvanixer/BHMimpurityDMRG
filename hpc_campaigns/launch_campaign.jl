#!/usr/bin/env julia
# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

using YAML
using Printf
using Dates

function as_dict(x)
    x isa AbstractDict || throw(ArgumentError("expected dictionary, got $(typeof(x))"))
    return Dict{String,Any}(String(k) => v for (k, v) in x)
end

function set_nested!(d::AbstractDict, dotted_key::AbstractString, value)
    parts = split(String(dotted_key), ".")
    isempty(parts) && throw(ArgumentError("empty key"))
    cur = d
    for p in parts[1:(end - 1)]
        if !haskey(cur, p) || !(cur[p] isa AbstractDict)
            cur[p] = Dict{String,Any}()
        end
        cur = cur[p]
    end
    cur[parts[end]] = value
    return d
end

function as_vector(v)
    return v isa AbstractVector ? collect(v) : [v]
end

function sweep_assignments(sweep::AbstractDict)
    sweepd = Dict{String,Any}(String(k) => v for (k, v) in sweep)
    keys_sorted = sort(collect(keys(sweepd)))
    if isempty(keys_sorted)
        return [Dict{String,Any}()]
    end
    values = [as_vector(sweepd[k]) for k in keys_sorted]
    assigns = Dict{String,Any}[]
    for tuple_vals in Iterators.product(values...)
        a = Dict{String,Any}()
        for (k, v) in zip(keys_sorted, tuple_vals)
            a[k] = v
        end
        push!(assigns, a)
    end
    return assigns
end

function apply_overrides!(cfg::AbstractDict, overrides::AbstractDict)
    for (k, v) in overrides
        set_nested!(cfg, String(k), v)
    end
    return cfg
end

function main()
    length(ARGS) >= 1 ||
        error("Usage: julia hpc_campaigns/launch_campaign.jl <campaign.yaml> [output_root_override]")
    campaign_file = abspath(ARGS[1])
    isfile(campaign_file) || error("Campaign file not found: $campaign_file")

    c = as_dict(YAML.load_file(campaign_file))
    campaign = haskey(c, "campaign") ? as_dict(c["campaign"]) : Dict{String,Any}()

    campaign_name = get(campaign, "name", "campaign_" * Dates.format(now(), "yyyymmdd_HHMMSS"))
    output_root_raw = if length(ARGS) >= 2
        String(ARGS[2])
    elseif haskey(campaign, "output_root_abs")
        String(campaign["output_root_abs"])
    else
        String(get(campaign, "output_root", "runs"))
    end
    output_root = isabspath(output_root_raw) ? output_root_raw :
                  abspath(joinpath(dirname(campaign_file), output_root_raw))
    app_script = get(campaign, "app_script", "scripts/solve_ground_state.jl")
    params_filename = get(campaign, "params_filename", "parameters.yaml")

    base_config_path = haskey(c, "base_config") ? String(c["base_config"]) :
                       error("Missing required key: base_config")
    base_config_abs = abspath(base_config_path)
    isfile(base_config_abs) || error("Base config not found: $base_config_abs")
    base_cfg = as_dict(YAML.load_file(base_config_abs))

    overrides = haskey(c, "overrides") ? as_dict(c["overrides"]) : Dict{String,Any}()
    sweep = haskey(c, "sweep") ? as_dict(c["sweep"]) : Dict{String,Any}()
    assigns = sweep_assignments(sweep)

    campaign_dir = joinpath(output_root, campaign_name)
    mkpath(campaign_dir)

    index_path = joinpath(campaign_dir, "index.csv")
    run_dirs_path = joinpath(campaign_dir, "run_dirs.txt")
    open(index_path, "w") do io
        println(io, "run_id,run_dir,params_path,app_script,status")
    end
    open(run_dirs_path, "w") do _ end

    for (i, assign) in enumerate(assigns)
        run_id = @sprintf("run_%04d", i)
        run_dir = joinpath(campaign_dir, run_id)
        mkpath(run_dir)

        cfg = deepcopy(base_cfg)
        apply_overrides!(cfg, overrides)
        apply_overrides!(cfg, assign)

        # Force output files to stay in each run directory.
        set_nested!(cfg, "io.state_save_path", abspath(joinpath(run_dir, "dmrg_state.h5")))
        set_nested!(cfg, "io.results_path", abspath(joinpath(run_dir, "results.h5")))
        set_nested!(cfg, "io.log_path", abspath(joinpath(run_dir, "run.log")))
        set_nested!(cfg, "dmrg.checkpoint_path", abspath(joinpath(run_dir, "dmrg_state_checkpoint.h5")))

        # Keep run id in metadata for traceability.
        set_nested!(cfg, "meta.run_name", run_id)

        params_path = abspath(joinpath(run_dir, params_filename))
        YAML.write_file(params_path, cfg)

        open(index_path, "a") do io
            println(io, join([run_id, run_dir, params_path, app_script, "PENDING"], ","))
        end
        open(run_dirs_path, "a") do io
            println(io, run_dir)
        end
    end

    println("Campaign generated:")
    println("  name: ", campaign_name)
    println("  runs: ", length(assigns))
    println("  dir:  ", campaign_dir)
    println("  index:", index_path)
    println("Next step:")
    println("  bash hpc_campaigns/slurm/submit_array.sh ", campaign_dir)
end

main()
