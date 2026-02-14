#!/usr/bin/env julia
# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

using Printf
using Dates
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

function main()
    length(ARGS) >= 1 ||
        error("Usage: julia hpc_campaigns/launch_campaign.jl <campaign.yaml> [output_root_override]")
    campaign_file = abspath(ARGS[1])
    isfile(campaign_file) || error("Campaign file not found: $campaign_file")

    c = load_ordered_yaml(campaign_file)
    campaign = haskey(c, "campaign") ? as_dict(c["campaign"]) : OrderedDict{String,Any}()

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
    assigns = sweep_assignments(sweep)

    campaign_dir = joinpath(output_root, campaign_name)
    mkpath(campaign_dir)

    index_path = joinpath(campaign_dir, "runs.csv")
    jobfile_path = joinpath(campaign_dir, "jobfile")
    run_dirs_path = joinpath(campaign_dir, "run_dirs.txt")
    open(index_path, "w") do io
        println(io, "run_id,run_dir,params_path,app_script,status")
    end
    open(jobfile_path, "w") do _ end
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
    println("  name: ", campaign_name)
    println("  runs: ", length(assigns))
    println("  dir:  ", campaign_dir)
    println("  runs_csv: ", index_path)
    println("  jobfile: ", jobfile_path)
    println("Next step:")
    println("  bash hpc_campaigns/slurm/submit_multilauncher.sh ", campaign_dir)
end

main()
