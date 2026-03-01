# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

using CSV
using DataFrames
using HDF5

const EARLY_STOP_MARKER = "Early stopping DMRG"
const WROTE_RESULTS_MARKER = "Wrote results"
const LOG_CHECKPOINT_SWEEP_REGEX = r"Wrote DMRG checkpoint at sweep\s+([0-9]+)"
const ERROR_PATTERNS = [
    r"(?m)^ERROR:",
    r"(?m)^Stacktrace:",
    r"(?i)segmentation fault",
    r"(?i)\bkilled\b"
]

@inline function as_int_or_nothing(x)
    if x === nothing
        return nothing
    elseif x isa Integer
        return Int(x)
    elseif x isa AbstractFloat
        return isfinite(x) ? Int(round(x)) : nothing
    elseif x isa AbstractString
        try
            return parse(Int, strip(x))
        catch
            return nothing
        end
    elseif x isa AbstractArray && length(x) == 1
        return as_int_or_nothing(first(x))
    end
    return nothing
end

@inline function as_string_or_nothing(x)
    x === nothing && return nothing
    try
        return String(x)
    catch
        return string(x)
    end
end

function print_help(io::IO=stdout)
    script = basename(@__FILE__)
    println(io, "Usage:")
    println(io, "  julia --startup-file=no --project=postprocess postprocess/$script [options] <campaign_root_or_runs_root> [output_csv]")
    println(io, "")
    println(io, "Options:")
    println(io, "  --absolute-paths       Write absolute paths in CSV (default: relative to campaign)")
    println(io, "  --quiet                Reduce progress output")
    println(io, "  -h, --help             Show this help")
    println(io, "")
    println(io, "Notes:")
    println(io, "  - If <campaign_root_or_runs_root> points to one campaign directory, analyze that campaign.")
    println(io, "  - Otherwise, analyze all immediate subdirectories that look like launch_campaign outputs.")
    println(io, "  - Convergence source priority: results.h5 diagnostics first, run.log fallback second.")
    println(io, "  - last_stored_sweep source priority: results diagnostics, then checkpoint.h5, then run.log.")
    println(io, "")
    println(io, "Output columns include:")
    println(io, "  campaign_name, run_id, run_dir, convergence_status, run_status, last_stored_sweep, and evidence.")
    println(io, "")
    println(io, "Examples:")
    println(io, "  julia --startup-file=no --project=postprocess postprocess/$script runs/deep_mi_scan_22feb2026_v1")
    println(io, "  ./bin/check_convergence runs")
end

function parse_args(args::Vector{String})
    show_help = false
    verbose = true
    absolute_paths = false
    positional = String[]

    i = 1
    while i <= length(args)
        a = args[i]
        if a in ("-h", "--help")
            show_help = true
            i += 1
        elseif a == "--absolute-paths"
            absolute_paths = true
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
        absolute_paths=absolute_paths,
        root=root,
        output_csv=output_csv
    )
end

@inline function clean_error(err)
    return replace(sprint(showerror, err), '\n' => ' ')
end

@inline function looks_like_run_dir_name(name::AbstractString)
    return occursin(r"^run_[0-9]+$", name)
end

function is_campaign_dir(path::AbstractString)
    isdir(path) || return false
    isfile(joinpath(path, "run_dirs.txt")) && return true
    isfile(joinpath(path, "runs.csv")) && return true
    for ent in readdir(path)
        if looks_like_run_dir_name(ent) && isdir(joinpath(path, ent))
            return true
        end
    end
    return false
end

function discover_campaign_dirs(root::AbstractString)
    root_abs = abspath(root)
    isdir(root_abs) || error("Root is not a directory: $root_abs")

    if is_campaign_dir(root_abs)
        return [root_abs]
    end

    dirs = String[]
    for ent in sort(readdir(root_abs))
        p = joinpath(root_abs, ent)
        isdir(p) || continue
        is_campaign_dir(p) || continue
        push!(dirs, p)
    end

    isempty(dirs) && error("No campaign directories found under: $root_abs")
    return dirs
end

function parse_runs_csv(path::AbstractString)
    out = String[]
    try
        for row in CSV.File(path)
            if hasproperty(row, :run_dir)
                rd = strip(String(getproperty(row, :run_dir)))
                isempty(rd) || push!(out, rd)
            end
        end
    catch err
        @warn "Failed to parse runs.csv, will fallback to directory scan" path=path error=clean_error(err)
    end
    return out
end

function discover_run_dirs(campaign_dir::AbstractString)
    run_dirs = String[]

    run_dirs_txt = joinpath(campaign_dir, "run_dirs.txt")
    if isfile(run_dirs_txt)
        for ln in eachline(run_dirs_txt)
            p = strip(ln)
            isempty(p) && continue
            push!(run_dirs, abspath(p))
        end
    end

    if isempty(run_dirs)
        runs_csv = joinpath(campaign_dir, "runs.csv")
        if isfile(runs_csv)
            append!(run_dirs, abspath.(parse_runs_csv(runs_csv)))
        end
    end

    if isempty(run_dirs)
        for ent in sort(readdir(campaign_dir))
            looks_like_run_dir_name(ent) || continue
            p = joinpath(campaign_dir, ent)
            isdir(p) || continue
            push!(run_dirs, abspath(p))
        end
    end

    unique!(run_dirs)
    sort!(run_dirs)
    return run_dirs
end

function read_diagnostics_converged(results_path::AbstractString)
    if !isfile(results_path)
        return (
            results_present=false,
            diagnostics_present=false,
            converged=nothing,
            last_stored_sweep=nothing,
            read_error=nothing
        )
    end

    try
        return HDF5.h5open(results_path, "r") do f
            has_diag = haskey(f, "diagnostics") && haskey(f["diagnostics"], "dmrg")
            if !has_diag
                return (
                    results_present=true,
                    diagnostics_present=false,
                    converged=nothing,
                    last_stored_sweep=nothing,
                    read_error=nothing
                )
            end
            g_dmrg = f["diagnostics"]["dmrg"]

            sweeps_completed = if haskey(g_dmrg, "sweeps_completed")
                as_int_or_nothing(read(g_dmrg["sweeps_completed"]))
            elseif haskey(g_dmrg, "sweep_trace")
                try
                    size(g_dmrg["sweep_trace"], 1)
                catch
                    nothing
                end
            else
                nothing
            end

            checkpoint_sweep_start = haskey(g_dmrg, "checkpoint_sweep_start") ?
                as_int_or_nothing(read(g_dmrg["checkpoint_sweep_start"])) : 0
            checkpoint_sweep_start === nothing && (checkpoint_sweep_start = 0)

            resume_mode = haskey(g_dmrg, "resume_mode") ?
                lowercase(strip(String(as_string_or_nothing(read(g_dmrg["resume_mode"]))))) : ""

            last_stored_sweep = if sweeps_completed === nothing
                checkpoint_sweep_start > 0 ? checkpoint_sweep_start : nothing
            elseif resume_mode == "remaining"
                checkpoint_sweep_start + sweeps_completed
            elseif sweeps_completed == 0 && checkpoint_sweep_start > 0
                checkpoint_sweep_start
            else
                sweeps_completed
            end

            if !haskey(g_dmrg, "converged")
                return (
                    results_present=true,
                    diagnostics_present=true,
                    converged=nothing,
                    last_stored_sweep=last_stored_sweep,
                    read_error=nothing
                )
            end
            raw = read(g_dmrg["converged"])
            val = if raw isa Bool
                raw
            elseif raw isa Integer
                raw != 0
            elseif raw isa AbstractString
                lowercase(strip(raw)) in ("1", "true", "yes", "y", "on")
            elseif raw isa AbstractArray && length(raw) == 1
                x = first(raw)
                x isa Bool ? x : (x isa Integer ? x != 0 : nothing)
            else
                nothing
            end
            return (
                results_present=true,
                diagnostics_present=true,
                converged=val,
                last_stored_sweep=last_stored_sweep,
                read_error=val === nothing ? "unparseable_diagnostics_converged" : nothing
            )
        end
    catch err
        return (
            results_present=true,
            diagnostics_present=false,
            converged=nothing,
            last_stored_sweep=nothing,
            read_error=clean_error(err)
        )
    end
end

function read_checkpoint_sweep(checkpoint_path::AbstractString)
    if !isfile(checkpoint_path)
        return (
            checkpoint_present=false,
            checkpoint_sweep=nothing,
            read_error=nothing
        )
    end

    try
        return HDF5.h5open(checkpoint_path, "r") do f
            sweep = nothing
            if haskey(f, "meta") && haskey(f["meta"], "checkpoint_sweep")
                sweep = as_int_or_nothing(read(f["meta"]["checkpoint_sweep"]))
            elseif haskey(f, "checkpoint_sweep")
                # Backward compatibility for older checkpoint layouts.
                sweep = as_int_or_nothing(read(f["checkpoint_sweep"]))
            end
            return (
                checkpoint_present=true,
                checkpoint_sweep=sweep,
                read_error=sweep === nothing ? "checkpoint_sweep_missing_or_unparseable" : nothing
            )
        end
    catch err
        return (
            checkpoint_present=true,
            checkpoint_sweep=nothing,
            read_error=clean_error(err)
        )
    end
end

function read_log_signals(log_path::AbstractString)
    if !isfile(log_path)
        return (
            log_present=false,
            early_stop=false,
            wrote_results=false,
            has_error=false,
            checkpoint_sweep=nothing,
            read_error=nothing
        )
    end

    try
        txt = read(log_path, String)
        early_stop = occursin(EARLY_STOP_MARKER, txt)
        wrote_results = occursin(WROTE_RESULTS_MARKER, txt)
        has_error = any(p -> occursin(p, txt), ERROR_PATTERNS)
        checkpoint_sweep = nothing
        for m in eachmatch(LOG_CHECKPOINT_SWEEP_REGEX, txt)
            sw = as_int_or_nothing(m.captures[1])
            sw === nothing && continue
            checkpoint_sweep = checkpoint_sweep === nothing ? sw : max(checkpoint_sweep, sw)
        end
        return (
            log_present=true,
            early_stop=early_stop,
            wrote_results=wrote_results,
            has_error=has_error,
            checkpoint_sweep=checkpoint_sweep,
            read_error=nothing
        )
    catch err
        return (
            log_present=true,
            early_stop=false,
            wrote_results=false,
            has_error=false,
            checkpoint_sweep=nothing,
            read_error=clean_error(err)
        )
    end
end

function assess_run(
    campaign_name::AbstractString,
    campaign_dir::AbstractString,
    run_dir::AbstractString;
    absolute_paths::Bool=false
)
    campaign_dir_abs = abspath(campaign_dir)
    run_dir_abs = abspath(run_dir)
    run_id = basename(normpath(run_dir))
    results_path = joinpath(run_dir_abs, "results.h5")
    checkpoint_path = joinpath(run_dir_abs, "dmrg_state_checkpoint.h5")
    log_path = joinpath(run_dir_abs, "run.log")

    diag = read_diagnostics_converged(results_path)
    checkpoint = read_checkpoint_sweep(checkpoint_path)
    log = read_log_signals(log_path)

    convergence_status = "unknown"
    used_log_fallback = false
    if diag.converged !== nothing
        convergence_status = diag.converged ? "converged" : "not_converged"
    elseif log.log_present && log.early_stop && !log.has_error
        convergence_status = "converged"
        used_log_fallback = true
    end

    run_status = "unknown"
    if diag.converged !== nothing
        run_status = "ok"
    elseif log.has_error && !log.wrote_results
        run_status = "failed"
    elseif diag.results_present || log.wrote_results
        run_status = "ok"
    elseif !diag.results_present && !log.log_present
        run_status = "missing"
    end

    last_stored_sweep = nothing
    if diag.last_stored_sweep !== nothing
        last_stored_sweep = diag.last_stored_sweep
    elseif checkpoint.checkpoint_sweep !== nothing
        last_stored_sweep = checkpoint.checkpoint_sweep
    elseif log.checkpoint_sweep !== nothing
        last_stored_sweep = log.checkpoint_sweep
    end

    evidence = String[]
    if diag.converged !== nothing
        push!(evidence, "results.h5:/diagnostics/dmrg/converged=$(diag.converged)")
    elseif diag.results_present && diag.diagnostics_present
        push!(evidence, "results.h5 has diagnostics but converged flag unavailable")
    elseif diag.results_present
        push!(evidence, "results.h5 missing diagnostics/dmrg")
    else
        push!(evidence, "results.h5 missing")
    end

    if log.log_present
        log.early_stop && push!(evidence, "run.log contains '$EARLY_STOP_MARKER'")
        log.wrote_results && push!(evidence, "run.log contains '$WROTE_RESULTS_MARKER'")
        log.checkpoint_sweep !== nothing && push!(evidence, "run.log checkpoint sweep=$(log.checkpoint_sweep)")
        log.has_error && push!(evidence, "run.log contains error signature")
    else
        push!(evidence, "run.log missing")
    end

    if diag.last_stored_sweep !== nothing
        push!(evidence, "last_stored_sweep from results.h5 diagnostics=$(diag.last_stored_sweep)")
    elseif checkpoint.checkpoint_sweep !== nothing
        push!(evidence, "last_stored_sweep from checkpoint.h5 meta/checkpoint_sweep=$(checkpoint.checkpoint_sweep)")
    end

    diag.read_error !== nothing && push!(evidence, "results_read_error=$(diag.read_error)")
    checkpoint.read_error !== nothing && push!(evidence, "checkpoint_read_error=$(checkpoint.read_error)")
    log.read_error !== nothing && push!(evidence, "log_read_error=$(log.read_error)")

    run_dir_out = absolute_paths ? run_dir_abs : relpath(run_dir_abs, campaign_dir_abs)

    return (
        campaign_name=String(campaign_name),
        run_id=String(run_id),
        run_dir=String(run_dir_out),
        convergence_status=convergence_status,
        run_status=run_status,
        last_stored_sweep=last_stored_sweep === nothing ? missing : last_stored_sweep,
        converged_from_diagnostics=diag.converged === nothing ? missing : diag.converged,
        used_log_fallback=used_log_fallback,
        results_present=diag.results_present,
        diagnostics_present=diag.diagnostics_present,
        log_present=log.log_present,
        log_has_error=log.has_error,
        log_early_stop=log.early_stop,
        log_wrote_results=log.wrote_results,
        evidence=join(evidence, "; ")
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
            push!(rows, assess_run(campaign_name, campaign_dir, run_dir; absolute_paths=opts.absolute_paths))
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
