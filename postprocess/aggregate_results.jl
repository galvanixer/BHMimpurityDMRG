# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

using Arrow
using CSV
using DataFrames
using Dates
using HDF5
using JLD2
using YAML

include(joinpath(@__DIR__, "parsers.jl"))
if !isdefined(@__MODULE__, :_POSTPROCESS_CONVERGENCE_CORE_INCLUDED)
    const _POSTPROCESS_CONVERGENCE_CORE_INCLUDED = true
    include(joinpath(@__DIR__, "convergence_core.jl"))
end

const PREFERRED_SUMMARY_COLUMN_ORDER = [
    "campaign_name",
    "run_id",
    "status",
    "convergence_status",
    "last_stored_sweep",
    "results_present",
    "results_source",
    "evidence",
    "cfg_meta_run_name",
    "cfg_lattice_L",
    "cfg_lattice_periodic",
    "cfg_initial_Na_total",
    "cfg_initial_Nb_total",
    "cfg_initial_impurity_distribution",
    "cfg_initial_seed",
    "cfg_local_nmax_a",
    "cfg_local_nmax_b",
    "cfg_t_a",
    "cfg_t_b",
    "cfg_U_a",
    "cfg_U_b",
    "cfg_U_ab",
    "cfg_mu_a",
    "cfg_mu_b",
    "cfg_dmrg_nsweeps",
    "cfg_dmrg_maxdim_max",
    "E0",
    "Na",
    "Nb",
    "na_mean",
    "nb_mean",
    "L_from_density",
    "has_observables",
    "has_density_density",
    "has_structure_factor",
    "has_single_particle_density_matrix",
    "has_triple_corr",
    "has_sampled_configs",
    "schema_id",
    "schema_version",
    "issues"
]

function add_convergence_row_fields!(row::Dict{String,Any}, conv)
    row["convergence_status"] = conv.convergence_status
    row["last_stored_sweep"] = conv.last_stored_sweep
    row["results_present"] = conv.results_present
    row["evidence"] = conv.evidence
    return row
end

function print_help(io::IO=stdout)
    script = basename(@__FILE__)
    println(io, "Usage:")
    println(io, "  julia --startup-file=no --project=postprocess postprocess/$script [options] <campaign_root> [output_jld2]")
    println(io, "")
    println(io, "Options:")
    println(io, "  --profile <name>       Extraction profile name from extract_profile.yaml (default: summary_only)")
    println(io, "  --profile-file <path>  YAML file containing profiles (default: postprocess/extract_profile.yaml)")
    println(io, "  --pattern <regex>      Regex for result files within campaign_root (default: ^results.*\\\\.(h5|hdf5)\$)")
    println(io, "  --quiet                Reduce progress output")
    println(io, "  -h, --help             Show this help")
    println(io, "")
    println(io, "Outputs:")
    println(io, "  <output_jld2>                    Aggregated nested data")
    println(io, "  Summary_<campaign_name>.arrow    Tabular summary export")
    println(io, "  Summary_<campaign_name>.csv      Tabular summary export")
    println(io, "")
    println(io, "Examples:")
    println(io, "  julia --startup-file=no --project=postprocess postprocess/$script runs/my_campaign")
    println(io, "  julia --startup-file=no --project=postprocess postprocess/$script --profile full runs/my_campaign runs/my_campaign/Results_my_campaign_full.jld2")
    println(io, "  ./bin/aggregate_results --profile summary_only runs/my_campaign")
end

function parse_args(args)
    campaign_root = nothing
    output_path = nothing
    profile_name = "summary_only"
    profile_file = joinpath(@__DIR__, "extract_profile.yaml")
    pattern_str = raw"^results.*\.(h5|hdf5)$"
    verbose = true
    show_help = false

    positional = String[]
    i = 1
    while i <= length(args)
        a = args[i]
        if a in ("-h", "--help")
            show_help = true
            i += 1
        elseif a == "--profile"
            i < length(args) || error("--profile requires a value")
            profile_name = args[i + 1]
            i += 2
        elseif a == "--profile-file"
            i < length(args) || error("--profile-file requires a value")
            profile_file = args[i + 1]
            i += 2
        elseif a == "--pattern"
            i < length(args) || error("--pattern requires a value")
            pattern_str = args[i + 1]
            i += 2
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

    if !isempty(positional)
        campaign_root = positional[1]
    end
    if length(positional) >= 2
        output_path = positional[2]
    end
    if length(positional) > 2
        error("Expected at most 2 positional arguments, got $(length(positional))")
    end

    return (
        campaign_root=campaign_root,
        output_path=output_path,
        profile_name=profile_name,
        profile_file=profile_file,
        pattern_str=pattern_str,
        verbose=verbose,
        show_help=show_help
    )
end

function load_profile(profile_file::AbstractString, profile_name::AbstractString)
    isfile(profile_file) || error("Profile file not found: $profile_file")
    raw = normalize_yaml(YAML.load_file(profile_file))
    profiles = get(raw, "profiles", Dict{String,Any}())
    profiles isa AbstractDict || error("Invalid profile file: missing 'profiles' dictionary")
    haskey(profiles, profile_name) || error("Profile '$profile_name' not found in $profile_file")
    profile = normalize_yaml(profiles[profile_name])
    return profile, collect(keys(profiles))
end

function discover_results_files(campaign_root::AbstractString, pattern::Regex)
    files = String[]
    for (root, _, fns) in walkdir(campaign_root)
        for fn in fns
            occursin(pattern, fn) || continue
            push!(files, joinpath(root, fn))
        end
    end
    sort!(files)
    return files
end

function group_results_files_by_run_dir(results_files::Vector{String})
    grouped = Dict{String,Vector{String}}()
    for path in results_files
        run_dir = dirname(abspath(path))
        if !haskey(grouped, run_dir)
            grouped[run_dir] = String[]
        end
        push!(grouped[run_dir], abspath(path))
    end
    for paths in values(grouped)
        sort!(paths)
    end
    return grouped
end

function pick_preferred_checkpoint_results(run_dir::AbstractString, run_files::Vector{String})
    preferred_names = ("results_from_checkpoint.h5", "results_checkpoint.h5")
    for name in preferred_names
        p = joinpath(run_dir, name)
        if p in run_files
            return p
        end
    end

    for p in run_files
        if occursin("checkpoint", lowercase(basename(p)))
            return p
        end
    end
    return nothing
end

function select_results_file_for_run(run_dir::AbstractString, run_files::Vector{String}, conv)
    main_results = joinpath(run_dir, "results.h5")
    has_main = main_results in run_files
    checkpoint_results = pick_preferred_checkpoint_results(run_dir, run_files)
    is_converged = getproperty(conv, :convergence_status) == "converged"

    # Policy: converged runs use canonical results.h5; otherwise prefer checkpoint-derived results.
    if is_converged && has_main
        return (path=main_results, source="results.h5")
    end
    if checkpoint_results !== nothing
        return (path=checkpoint_results, source=basename(checkpoint_results))
    end
    if has_main
        return (path=main_results, source="results.h5")
    end
    return (path=run_files[1], source=basename(run_files[1]))
end

function derive_candidate_run_id(campaign_root::AbstractString, results_path::AbstractString, idx::Int)
    run_dir = basename(dirname(results_path))
    if occursin(r"^run_[0-9]+$", run_dir)
        return run_dir
    end
    rel_dir = relpath(dirname(results_path), campaign_root)
    base = splitext(basename(results_path))[1]
    if rel_dir == "."
        return isempty(base) ? "run_$(lpad(idx, 4, '0'))" : base
    end
    rel_clean = replace(rel_dir, r"[\\/]+" => "__")
    if isempty(base)
        return rel_clean
    end
    return rel_clean * "__" * base
end

function build_unique_run_ids(campaign_root::AbstractString, results_paths::Vector{String})
    used = Dict{String,Int}()
    ids = String[]
    for (idx, path) in enumerate(results_paths)
        base = derive_candidate_run_id(campaign_root, path, idx)
        n = get(used, base, 0) + 1
        used[base] = n
        push!(ids, n == 1 ? base : "$(base)_$n")
    end
    return ids
end

function to_serializable(x)
    if x isa AbstractDict
        out = Dict{String,Any}()
        for (k, v) in x
            out[String(k)] = to_serializable(v)
        end
        return out
    elseif x isa NamedTuple
        out = Dict{String,Any}()
        for (k, v) in pairs(x)
            out[String(k)] = to_serializable(v)
        end
        return out
    elseif x isa AbstractVector
        if eltype(x) <: Number || eltype(x) <: AbstractString || eltype(x) <: Bool
            return x
        end
        return [to_serializable(v) for v in x]
    elseif x isa AbstractArray
        if eltype(x) <: Number || eltype(x) <: Bool || eltype(x) <: AbstractString
            return x
        end
        return map(to_serializable, x)
    elseif x === nothing || x isa Number || x isa Bool || x isa AbstractString
        return x
    else
        return string(x)
    end
end

@inline function normalize_summary_cell(x)
    if x === nothing
        return missing
    elseif x isa Missing
        return missing
    elseif x isa Bool
        return x
    elseif x isa Integer
        return Int64(x)
    elseif x isa AbstractFloat
        return Float64(x)
    elseif x isa AbstractString
        return String(x)
    elseif x isa AbstractVector
        return join(string.(x), "; ")
    elseif x isa AbstractDict || x isa NamedTuple
        return sprint(show, x)
    else
        return string(x)
    end
end

function typed_summary_column(values::Vector{Any})
    nonmissing = [v for v in values if v !== missing]
    if isempty(nonmissing)
        return fill(missing, length(values))
    end

    has_bool = any(v -> v isa Bool, nonmissing)
    has_int = any(v -> (v isa Integer) && !(v isa Bool), nonmissing)
    has_float = any(v -> v isa Float64, nonmissing)
    has_string = any(v -> v isa String, nonmissing)
    has_numeric = has_int || has_float
    kind_count = count(identity, (has_bool, has_numeric, has_string))

    if kind_count == 1 && has_bool
        col = Vector{Union{Missing,Bool}}(undef, length(values))
        for i in eachindex(values)
            v = values[i]
            col[i] = (v === missing) ? missing : Bool(v)
        end
        return col
    elseif kind_count == 1 && has_int
        col = Vector{Union{Missing,Int64}}(undef, length(values))
        for i in eachindex(values)
            v = values[i]
            col[i] = (v === missing) ? missing : Int64(v)
        end
        return col
    elseif kind_count == 1 && has_float
        col = Vector{Union{Missing,Float64}}(undef, length(values))
        for i in eachindex(values)
            v = values[i]
            col[i] = (v === missing) ? missing : Float64(v)
        end
        return col
    elseif kind_count == 1 && has_numeric
        col = Vector{Union{Missing,Float64}}(undef, length(values))
        for i in eachindex(values)
            v = values[i]
            col[i] = (v === missing) ? missing : Float64(v)
        end
        return col
    else
        col = Vector{Union{Missing,String}}(undef, length(values))
        for i in eachindex(values)
            v = values[i]
            col[i] = (v === missing) ? missing : string(v)
        end
        return col
    end
end

function summary_rows_to_dataframe(summary_rows::Vector{Dict{String,Any}})
    isempty(summary_rows) && return DataFrame()
    keyset = Set{String}()
    for row in summary_rows
        union!(keyset, keys(row))
    end

    ordered_keys = String[]
    for key in PREFERRED_SUMMARY_COLUMN_ORDER
        if key in keyset
            push!(ordered_keys, key)
            delete!(keyset, key)
        end
    end
    append!(ordered_keys, sort!(collect(keyset)))

    n = length(summary_rows)
    df = DataFrame()

    for key in ordered_keys
        raw = Vector{Any}(undef, n)
        for i in 1:n
            raw[i] = normalize_summary_cell(get(summary_rows[i], key, missing))
        end
        df[!, Symbol(key)] = typed_summary_column(raw)
    end

    return df
end

function summary_table_paths(output_path::AbstractString, campaign_name::AbstractString)
    out_dir = dirname(abspath(output_path))
    return (
        arrow_path=joinpath(out_dir, "Summary_$(campaign_name).arrow"),
        csv_path=joinpath(out_dir, "Summary_$(campaign_name).csv")
    )
end

function aggregate_results(
    campaign_root::AbstractString;
    output_path::AbstractString,
    profile::AbstractDict,
    profile_name::AbstractString,
    pattern::Regex=Regex(raw"^results.*\.(h5|hdf5)$"),
    verbose::Bool=true
)
    campaign_root_abs = abspath(campaign_root)
    isdir(campaign_root_abs) || error("Campaign root is not a directory: $campaign_root_abs")
    campaign_name = basename(normpath(campaign_root_abs))

    all_results_files = discover_results_files(campaign_root_abs, pattern)
    isempty(all_results_files) && error("No results files found under $campaign_root_abs matching regex: $(pattern.pattern)")

    files_by_run_dir = group_results_files_by_run_dir(all_results_files)
    run_dirs = sort!(collect(keys(files_by_run_dir)))

    selected_paths = String[]
    selected_sources = String[]
    selected_convergence = Any[]
    for run_dir in run_dirs
        conv = assess_run(campaign_name, run_dir)
        selection = select_results_file_for_run(run_dir, files_by_run_dir[run_dir], conv)
        push!(selected_paths, selection.path)
        push!(selected_sources, selection.source)
        push!(selected_convergence, conv)
    end

    run_ids = build_unique_run_ids(campaign_root_abs, selected_paths)
    runs = Dict{String,Any}()
    summary_rows = Vector{Dict{String,Any}}()
    schema_counts = Dict{String,Int}()
    n_success = 0
    n_error = 0

    for idx in eachindex(run_ids)
        run_id = run_ids[idx]
        path = selected_paths[idx]
        source = selected_sources[idx]
        verbose && println("[$idx/$(length(selected_paths))] Parsing $(abspath(path))")
        abs_path = abspath(path)
        run_dir = dirname(abs_path)
        conv = selected_convergence[idx]
        try
            parsed = parse_results_file(abs_path; profile=profile)
            schema = get(parsed, "schema", Dict{String,Any}())
            schema_id = get(schema, "results_schema_id", "unknown")
            schema_version = get(schema, "results_schema_version", "unknown")
            schema_key = "$(schema_id)@$(schema_version)"
            schema_counts[schema_key] = get(schema_counts, schema_key, 0) + 1

            run_record = Dict{String,Any}(
                "run_id" => run_id,
                "results_path" => abs_path,
                "run_dir" => run_dir,
                "schema" => schema,
                "meta" => get(parsed, "meta", Dict{String,Any}()),
                "summary" => get(parsed, "summary", Dict{String,Any}()),
                "observables" => get(parsed, "observables", Dict{String,Any}()),
                "issues" => get(parsed, "issues", String[]),
                "results_source" => source,
                "convergence" => conv,
                "status" => "ok"
            )
            runs[run_id] = run_record

            row = Dict{String,Any}()
            merge!(row, run_record["summary"])
            row["campaign_name"] = campaign_name
            row["run_id"] = run_id
            row["status"] = "ok"
            row["issues"] = run_record["issues"]
            row["results_source"] = source
            add_convergence_row_fields!(row, conv)
            push!(summary_rows, row)
            n_success += 1
        catch err
            err_msg = sprint(showerror, err)
            row = Dict{String,Any}(
                "campaign_name" => campaign_name,
                "run_id" => run_id,
                "status" => "error",
                "error" => err_msg,
                "issues" => ["parser_exception"],
                "results_source" => source
            )
            add_convergence_row_fields!(row, conv)
            push!(summary_rows, row)
            runs[run_id] = Dict{String,Any}(
                "run_id" => run_id,
                "results_path" => abs_path,
                "run_dir" => run_dir,
                "results_source" => source,
                "convergence" => conv,
                "status" => "error",
                "error" => err_msg
            )
            n_error += 1
        end
    end

    sort!(summary_rows, by=r -> String(get(r, "run_id", "")))

    table_paths = summary_table_paths(output_path, campaign_name)
    summary_df = summary_rows_to_dataframe(summary_rows)
    Arrow.write(table_paths.arrow_path, summary_df)
    CSV.write(table_paths.csv_path, summary_df)

    manifest = Dict{String,Any}(
        "created_at" => string(Dates.now()),
        "aggregator_name" => "postprocess.aggregate_results",
        "aggregator_version" => "0.1.0",
        "campaign_root" => campaign_root_abs,
        "output_path" => abspath(output_path),
        "campaign_name" => campaign_name,
        "summary_arrow_path" => table_paths.arrow_path,
        "summary_csv_path" => table_paths.csv_path,
        "profile_name" => profile_name,
        "profile" => profile,
        "results_file_pattern" => pattern.pattern,
        "n_discovered_files" => length(all_results_files),
        "n_selected_runs" => length(selected_paths),
        "n_success" => n_success,
        "n_error" => n_error,
        "schema_counts" => schema_counts,
        "results_files" => abspath.(all_results_files),
        "selected_results_files" => abspath.(selected_paths),
        "selection_policy" => "run-level selection: converged->results.h5, else prefer results_from_checkpoint.h5/results_checkpoint.h5"
    )

    mkpath(dirname(abspath(output_path)))
    JLD2.jldopen(output_path, "w") do f
        f["manifest"] = to_serializable(manifest)
        f["summary"] = to_serializable(summary_rows)
        f["runs"] = to_serializable(runs)
    end

    return (
        output_path=abspath(output_path),
        summary_arrow_path=table_paths.arrow_path,
        summary_csv_path=table_paths.csv_path,
        n_discovered=length(selected_paths),
        n_discovered_files=length(all_results_files),
        n_success=n_success,
        n_error=n_error,
        schema_counts=schema_counts
    )
end

function main()
    opts = parse_args(ARGS)
    if opts.show_help
        print_help()
        return nothing
    end
    opts.campaign_root === nothing && error("campaign_root is required. Use --help for usage.")

    pattern = Regex(opts.pattern_str)
    profile, available_profiles = load_profile(opts.profile_file, opts.profile_name)
    campaign_root_abs = abspath(opts.campaign_root)
    campaign_name = basename(normpath(campaign_root_abs))
    output_path = opts.output_path === nothing ?
                  joinpath(campaign_root_abs, "Results_$(campaign_name).jld2") :
                  opts.output_path

    if opts.verbose
        println("Campaign root      : $(campaign_root_abs)")
        println("Campaign name      : $(campaign_name)")
        println("Output JLD         : $(abspath(output_path))")
        println("Profile            : $(opts.profile_name)")
        println("Profile file       : $(abspath(opts.profile_file))")
        println("Available profiles : $(join(sort(String.(available_profiles)), ", "))")
        println("Pattern            : $(opts.pattern_str)")
    end

    result = aggregate_results(
        opts.campaign_root;
        output_path=output_path,
        profile=profile,
        profile_name=opts.profile_name,
        pattern=pattern,
        verbose=opts.verbose
    )

    println("Wrote aggregate: $(result.output_path)")
    println("Wrote summary arrow: $(result.summary_arrow_path)")
    println("Wrote summary csv: $(result.summary_csv_path)")
    println("Runs discovered: $(result.n_discovered), results files discovered: $(result.n_discovered_files), parsed: $(result.n_success), errors: $(result.n_error)")
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
