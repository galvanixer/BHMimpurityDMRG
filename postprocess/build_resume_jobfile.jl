# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

using CSV
using DataFrames
using Dates

const DEFAULT_HALFWAY_RUN_STATUSES = Set(["failed", "unknown", "missing"])

function print_help(io::IO=stdout)
    script = basename(@__FILE__)
    println(io, "Usage:")
    println(io, "  julia --startup-file=no --project=postprocess postprocess/$script [options] <campaign_root> [convergence_csv] [output_dir]")
    println(io, "")
    println(io, "Options:")
    println(io, "  --mode <halfway|unconverged>   Selection mode (default: halfway)")
    println(io, "  --quiet                        Reduce progress output")
    println(io, "  -h, --help                     Show this help")
    println(io, "")
    println(io, "Mode details:")
    println(io, "  halfway     : convergence_status != converged, run_status in {failed,unknown,missing}, and last_stored_sweep present")
    println(io, "  unconverged : convergence_status != converged")
    println(io, "")
    println(io, "Defaults:")
    println(io, "  convergence_csv: <campaign_root>/Convergence_<campaign_name>.csv")
    println(io, "  output_dir     : <campaign_root>/resume_<yyyymmdd_HHMMSS>")
    println(io, "")
    println(io, "Outputs in output_dir:")
    println(io, "  jobfile               Resume-only commands")
    println(io, "  resume_manifest.csv   Selection + mapping audit")
    println(io, "")
    println(io, "Examples:")
    println(io, "  ./bin/build_resume_jobfile runs/my_campaign")
    println(io, "  ./bin/build_resume_jobfile --mode unconverged runs/my_campaign")
end

@inline function as_str(x)
    x === missing && return ""
    x === nothing && return ""
    return String(x)
end

@inline function lower_str(x)
    return lowercase(strip(as_str(x)))
end

@inline function has_nonmissing_value(x)
    if x === nothing || x === missing
        return false
    elseif x isa AbstractString
        return !isempty(strip(x))
    end
    return true
end

function parse_args(args::Vector{String})
    show_help = false
    verbose = true
    mode = "halfway"

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
        elseif a == "--mode"
            i < length(args) || error("--mode requires a value")
            mode = lowercase(strip(args[i + 1]))
            i += 2
        elseif startswith(a, "--")
            error("Unknown option: $a")
        else
            push!(positional, a)
            i += 1
        end
    end

    mode in ("halfway", "unconverged") || error("Unsupported mode: $mode (expected: halfway|unconverged)")

    campaign_root = nothing
    convergence_csv = nothing
    output_dir = nothing
    if !isempty(positional)
        campaign_root = positional[1]
    end
    if length(positional) >= 2
        convergence_csv = positional[2]
    end
    if length(positional) >= 3
        output_dir = positional[3]
    end
    if length(positional) > 3
        error("Expected at most 3 positional arguments, got $(length(positional))")
    end

    return (
        show_help=show_help,
        verbose=verbose,
        mode=mode,
        campaign_root=campaign_root,
        convergence_csv=convergence_csv,
        output_dir=output_dir
    )
end

@inline function required_columns(df::DataFrame, cols::Vector{String}, label::AbstractString)
    missing_cols = [c for c in cols if !(c in names(df))]
    isempty(missing_cols) || error("$label is missing required column(s): $(join(missing_cols, ", "))")
end

function select_resume_candidates(conv::DataFrame, mode::AbstractString)
    required_columns(conv, ["run_id", "convergence_status"], "convergence CSV")
    if mode == "halfway"
        required_columns(conv, ["run_status", "last_stored_sweep"], "convergence CSV")
    end

    selected = falses(nrow(conv))
    reasons = Vector{String}(undef, nrow(conv))

    for i in 1:nrow(conv)
        rid = strip(as_str(conv[i, :run_id]))
        cstat = lower_str(conv[i, :convergence_status])
        if isempty(rid)
            selected[i] = false
            reasons[i] = "skip_empty_run_id"
            continue
        end

        if mode == "unconverged"
            if cstat != "converged"
                selected[i] = true
                reasons[i] = "selected_unconverged"
            else
                selected[i] = false
                reasons[i] = "skip_converged"
            end
            continue
        end

        rstat = lower_str(conv[i, :run_status])
        sweep_ok = has_nonmissing_value(conv[i, :last_stored_sweep])
        not_converged = cstat != "converged"
        status_ok = rstat in DEFAULT_HALFWAY_RUN_STATUSES

        if not_converged && status_ok && sweep_ok
            selected[i] = true
            reasons[i] = "selected_halfway"
        else
            selected[i] = false
            reasons[i] = "skip_halfway_filter"
        end
    end

    out = copy(conv)
    out[!, :selected] = selected
    out[!, :selection_reason] = reasons
    return out
end

function build_runid_to_job_line_map(runs_df::DataFrame, job_lines::Vector{String})
    required_columns(runs_df, ["run_id"], "runs.csv")

    runid_to_line = Dict{String,Int}()
    for (i, row) in enumerate(eachrow(runs_df))
        rid = strip(as_str(row[:run_id]))
        isempty(rid) && continue
        if !haskey(runid_to_line, rid)
            runid_to_line[rid] = i
        end
    end

    return runid_to_line
end

function empty_manifest_df()
    return DataFrame(
        run_id=String[],
        selected=Bool[],
        included=Bool[],
        reason=String[],
        jobfile_line=Vector{Union{Missing,Int}}(),
        convergence_status=String[],
        run_status=String[],
        last_stored_sweep=Any[]
    )
end

function default_convergence_csv(campaign_root_abs::AbstractString)
    campaign_name = basename(normpath(campaign_root_abs))
    return joinpath(campaign_root_abs, "Convergence_$(campaign_name).csv")
end

function default_output_dir(campaign_root_abs::AbstractString)
    stamp = Dates.format(now(), dateformat"yyyymmdd_HHMMSS")
    return joinpath(campaign_root_abs, "resume_$stamp")
end

function main(args=ARGS)
    opts = parse_args(args)
    if opts.show_help
        print_help()
        return nothing
    end

    opts.campaign_root === nothing && error("campaign_root is required. Use --help for usage.")

    campaign_root_abs = abspath(opts.campaign_root)
    isdir(campaign_root_abs) || error("campaign_root is not a directory: $campaign_root_abs")

    convergence_csv = opts.convergence_csv === nothing ? default_convergence_csv(campaign_root_abs) : abspath(opts.convergence_csv)
    runs_csv = joinpath(campaign_root_abs, "runs.csv")
    src_jobfile = joinpath(campaign_root_abs, "jobfile")
    output_dir = opts.output_dir === nothing ? default_output_dir(campaign_root_abs) : abspath(opts.output_dir)

    isfile(convergence_csv) || error("Convergence CSV not found: $convergence_csv")
    isfile(runs_csv) || error("runs.csv not found: $runs_csv")
    isfile(src_jobfile) || error("jobfile not found: $src_jobfile")

    conv_df = CSV.read(convergence_csv, DataFrame)
    selected_df = select_resume_candidates(conv_df, opts.mode)
    rows_selected = selected_df[selected_df.selected .== true, :]

    job_lines = readlines(src_jobfile)
    runs_df = CSV.read(runs_csv, DataFrame)
    runid_to_line = build_runid_to_job_line_map(runs_df, job_lines)

    manifest_rows = NamedTuple[]
    out_cmds = String[]

    for row in eachrow(rows_selected)
        rid = strip(as_str(row[:run_id]))
        cstat = as_str(row[:convergence_status])
        rstat = hasproperty(row, :run_status) ? as_str(row[:run_status]) : ""
        sweep = hasproperty(row, :last_stored_sweep) ? row[:last_stored_sweep] : missing

        line_no = get(runid_to_line, rid, 0)
        if line_no <= 0
            push!(manifest_rows, (
                run_id=rid,
                selected=true,
                included=false,
                reason="run_id_not_found_in_runs_csv",
                jobfile_line=missing,
                convergence_status=cstat,
                run_status=rstat,
                last_stored_sweep=sweep
            ))
            continue
        end
        if line_no > length(job_lines)
            push!(manifest_rows, (
                run_id=rid,
                selected=true,
                included=false,
                reason="jobfile_line_out_of_range",
                jobfile_line=line_no,
                convergence_status=cstat,
                run_status=rstat,
                last_stored_sweep=sweep
            ))
            continue
        end

        cmd = strip(job_lines[line_no])
        if isempty(cmd)
            push!(manifest_rows, (
                run_id=rid,
                selected=true,
                included=false,
                reason="empty_jobfile_line",
                jobfile_line=line_no,
                convergence_status=cstat,
                run_status=rstat,
                last_stored_sweep=sweep
            ))
            continue
        end

        push!(out_cmds, cmd)
        push!(manifest_rows, (
            run_id=rid,
            selected=true,
            included=true,
            reason="included",
            jobfile_line=line_no,
            convergence_status=cstat,
            run_status=rstat,
            last_stored_sweep=sweep
        ))
    end

    mkpath(output_dir)
    out_jobfile = joinpath(output_dir, "jobfile")
    open(out_jobfile, "w") do io
        for cmd in out_cmds
            println(io, cmd)
        end
    end

    manifest_df = isempty(manifest_rows) ? empty_manifest_df() : DataFrame(manifest_rows)
    if nrow(manifest_df) > 0
        sort!(manifest_df, [:included, :jobfile_line], rev=[true, false])
    end
    manifest_path = joinpath(output_dir, "resume_manifest.csv")
    CSV.write(manifest_path, manifest_df)

    if opts.verbose
        n_candidates = nrow(rows_selected)
        n_included = isempty(manifest_df) ? 0 : sum(manifest_df.included .== true)
        n_excluded = n_candidates - n_included
        println("Campaign root        : $campaign_root_abs")
        println("Convergence CSV      : $convergence_csv")
        println("Mode                 : $(opts.mode)")
        println("Candidates selected  : $n_candidates")
        println("Commands included    : $n_included")
        println("Candidates excluded  : $n_excluded")
        println("Output directory     : $output_dir")
        println("Resume jobfile       : $out_jobfile")
        println("Resume manifest      : $manifest_path")
        println("")
        println("Submit with:")
        println("  bash hpc_campaigns/slurm/submit_multilauncher.sh $output_dir")
    end

    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
