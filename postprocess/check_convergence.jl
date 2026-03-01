# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

using CSV
using DataFrames

if !isdefined(@__MODULE__, :_POSTPROCESS_CONVERGENCE_CORE_INCLUDED)
    const _POSTPROCESS_CONVERGENCE_CORE_INCLUDED = true
    include(joinpath(@__DIR__, "convergence_core.jl"))
end

function print_help(io::IO=stdout)
    script = basename(@__FILE__)
    println(io, "Usage:")
    println(io, "  julia --startup-file=no --project=postprocess postprocess/$script [options] <campaign_root_or_runs_root> [output_csv]")
    println(io, "")
    println(io, "Options:")
    println(io, "  --quiet                Reduce progress output")
    println(io, "  -h, --help             Show this help")
    println(io, "")
    println(io, "Notes:")
    println(io, "  - If <campaign_root_or_runs_root> points to one campaign directory, analyze that campaign.")
    println(io, "  - Otherwise, analyze all immediate subdirectories that look like launch_campaign outputs.")
    println(io, "  - Uses convergence logic from postprocess/convergence_core.jl.")
    println(io, "")
    println(io, "Output columns include:")
    println(io, "  campaign_name, run_id, seed_initial_state, hamiltonian params, convergence_status,")
    println(io, "  run_status, last_stored_sweep, diagnostics/log signals, and evidence.")
    println(io, "")
    println(io, "Examples:")
    println(io, "  julia --startup-file=no --project=postprocess postprocess/$script runs/deep_mi_scan_22feb2026_v1")
    println(io, "  ./bin/check_convergence runs")
end

function parse_args(args::Vector{String})
    show_help = false
    verbose = true
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
        elseif startswith(a, "--")
            error("Unknown option: $a")
        else
            push!(positional, a)
            i += 1
        end
    end

    root = nothing
    output_csv = nothing
    if !isempty(positional)
        root = positional[1]
    end
    if length(positional) >= 2
        output_csv = positional[2]
    end
    if length(positional) > 2
        error("Expected at most 2 positional arguments, got $(length(positional))")
    end

    return (
        show_help=show_help,
        verbose=verbose,
        root=root,
        output_csv=output_csv
    )
end

function default_output_path(root_abs::AbstractString, campaign_dirs::Vector{String})
    if length(campaign_dirs) == 1
        campaign_name = basename(normpath(campaign_dirs[1]))
        return joinpath(campaign_dirs[1], "Convergence_$(campaign_name).csv")
    end
    return joinpath(root_abs, "Convergence_all_campaigns.csv")
end

function summarize(df::DataFrame; io::IO=stdout)
    n_total = nrow(df)
    n_converged = sum(df.convergence_status .== "converged")
    n_not_converged = sum(df.convergence_status .== "not_converged")
    n_unknown = sum(df.convergence_status .== "unknown")
    n_failed = sum(df.run_status .== "failed")

    println(io, "Runs analyzed        : $n_total")
    println(io, "Converged            : $n_converged")
    println(io, "Not converged        : $n_not_converged")
    println(io, "Unknown convergence  : $n_unknown")
    println(io, "Run failures         : $n_failed")
end

function main(args=ARGS)
    opts = parse_args(args)
    if opts.show_help
        print_help()
        return nothing
    end
    opts.root === nothing && error("campaign_root_or_runs_root is required. Use --help for usage.")

    root_abs = abspath(opts.root)
    campaign_dirs = discover_campaign_dirs(root_abs)
    rows = NamedTuple[]

    for campaign_dir in campaign_dirs
        campaign_name = basename(normpath(campaign_dir))
        run_dirs = discover_run_dirs(campaign_dir)
        if isempty(run_dirs)
            @warn "No run directories found in campaign" campaign_dir=campaign_dir
            continue
        end
        for run_dir in run_dirs
            push!(rows, assess_run(campaign_name, run_dir))
        end
    end

    isempty(rows) && error("No runs found to analyze under: $root_abs")

    df = DataFrame(rows)
    sort!(df, [:campaign_name, :run_id])

    out_path = opts.output_csv === nothing ? default_output_path(root_abs, campaign_dirs) : abspath(opts.output_csv)
    mkpath(dirname(out_path))
    CSV.write(out_path, df)

    if opts.verbose
        println("Campaign directories : $(length(campaign_dirs))")
        summarize(df)
        println("Output CSV           : $out_path")
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
