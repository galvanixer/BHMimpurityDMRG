# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

using CSV
using HDF5
using YAML

if !isdefined(@__MODULE__, :_POSTPROCESS_COMMON_INCLUDED)
    const _POSTPROCESS_COMMON_INCLUDED = true
    include(joinpath(@__DIR__, "common.jl"))
end

const EARLY_STOP_MARKER = "Early stopping DMRG"
const WROTE_RESULTS_MARKER = "Wrote results"
const LOG_CHECKPOINT_SWEEP_REGEX = r"Wrote DMRG checkpoint at sweep\s+([0-9]+)"
const ERROR_PATTERNS = [
    r"(?m)^ERROR:",
    r"(?m)^Stacktrace:",
    r"(?i)segmentation fault",
    r"(?i)\bkilled\b"
]

"""
    as_int_or_nothing(x)

Best-effort coercion of `x` to `Int`.
Returns `nothing` when conversion is not possible.
"""
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

"""
    clean_error(err)

Render an exception/error object as a single-line message suitable for CSV/log fields.
"""
@inline function clean_error(err)
    return replace(sprint(showerror, err), '\n' => ' ')
end

"""
    looks_like_run_dir_name(name)

Return `true` when `name` matches launch-style run directories like `run_0001`.
"""
@inline function looks_like_run_dir_name(name::AbstractString)
    return occursin(r"^run_[0-9]+$", name)
end

"""
    hparam_or_missing(cfg, key)

Read `hamiltonian.<key>` from normalized params config `cfg`.
Returns `missing` when unavailable.
"""
@inline function hparam_or_missing(cfg, key::AbstractString)
    cfg === nothing && return missing
    v = nested_get(cfg, ["hamiltonian", key], nothing)
    return v === nothing ? missing : v
end

"""
    read_params_yaml_from_results(results_path)

Load normalized parameter YAML from `results.h5:/meta/params_yaml`.
Returns `nothing` on missing file/dataset or parse/read failures.
"""
function read_params_yaml_from_results(results_path::AbstractString)
    if !isfile(results_path)
        return nothing
    end
    try
        return HDF5.h5open(results_path, "r") do f
            haskey(f, "meta") || return nothing
            g_meta = f["meta"]
            haskey(g_meta, "params_yaml") || return nothing
            raw = read(g_meta["params_yaml"])
            txt = as_string_or_nothing(raw)
            txt === nothing && return nothing
            parsed = YAML.load(txt)
            return parsed === nothing ? nothing : normalize_yaml(parsed)
        end
    catch
        return nothing
    end
end

"""
    guess_params_file(run_dir)

Heuristically pick a params YAML file from `run_dir`.
Prefers `parameters.yaml`, otherwise chooses among `.yaml/.yml` files.
Returns `nothing` if no YAML file exists.
"""
function guess_params_file(run_dir::AbstractString)
    default_path = joinpath(run_dir, "parameters.yaml")
    if isfile(default_path)
        return default_path
    end

    yamls = String[]
    for fn in sort(readdir(run_dir))
        p = joinpath(run_dir, fn)
        isfile(p) || continue
        low = lowercase(fn)
        (endswith(low, ".yaml") || endswith(low, ".yml")) || continue
        push!(yamls, p)
    end
    isempty(yamls) && return nothing
    if length(yamls) == 1
        return yamls[1]
    end
    for p in yamls
        occursin("param", lowercase(basename(p))) && return p
    end
    return yamls[1]
end

"""
    read_params_yaml_from_file(run_dir)

Load and normalize params YAML from the file selected by `guess_params_file(run_dir)`.
Returns `nothing` when no candidate exists or parsing fails.
"""
function read_params_yaml_from_file(run_dir::AbstractString)
    path = guess_params_file(run_dir)
    path === nothing && return nothing
    try
        parsed = YAML.load_file(path)
        return parsed === nothing ? nothing : normalize_yaml(parsed)
    catch
        return nothing
    end
end

"""
    read_hamiltonian_params(run_dir, results_path)

Read run-level parameter fields used by convergence reporting:
`seed_initial_state`, `t_a`, `t_b`, `U_a`, `U_b`, `U_ab`, `mu_a`, `mu_b`.

Source priority:
1. `results_path:/meta/params_yaml`
2. YAML file in `run_dir`
"""
function read_hamiltonian_params(run_dir::AbstractString, results_path::AbstractString)
    cfg = read_params_yaml_from_results(results_path)
    if cfg === nothing
        cfg = read_params_yaml_from_file(run_dir)
    end

    return (
        seed_initial_state=begin
            v = cfg === nothing ? nothing : nested_get(cfg, ["initial_state", "seed"], nothing)
            v === nothing ? missing : v
        end,
        t_a=hparam_or_missing(cfg, "t_a"),
        t_b=hparam_or_missing(cfg, "t_b"),
        U_a=hparam_or_missing(cfg, "U_a"),
        U_b=hparam_or_missing(cfg, "U_b"),
        U_ab=hparam_or_missing(cfg, "U_ab"),
        mu_a=hparam_or_missing(cfg, "mu_a"),
        mu_b=hparam_or_missing(cfg, "mu_b")
    )
end

"""
    is_campaign_dir(path)

Return `true` if `path` looks like a launch_campaign output directory
(has `run_dirs.txt` / `runs.csv` or `run_XXXX` subdirectories).
"""
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

"""
    discover_campaign_dirs(root)

If `root` is a campaign directory, return `[root]`.
Otherwise return all immediate subdirectories that look like campaign roots.
Errors if none are found.
"""
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

"""
    parse_runs_csv(path)

Read `run_dir` entries from campaign `runs.csv`.
Returns an empty vector on parse failure (caller can fallback to directory scan).
"""
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

"""
    discover_run_dirs(campaign_dir)

Discover absolute run directories for one campaign with fallback order:
1. `run_dirs.txt`
2. `runs.csv` `run_dir` column
3. direct `run_XXXX` directory scan
"""
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

"""
    read_diagnostics_converged(results_path)

Read convergence diagnostics from `results.h5` (`/diagnostics/dmrg`).
Returns a NamedTuple with:
- `results_present`, `diagnostics_present`
- `converged` (`Bool` or `nothing`)
- `last_stored_sweep` (`Int` or `nothing`)
- `read_error` (`String` or `nothing`)
"""
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

"""
    read_checkpoint_sweep(checkpoint_path)

Read `checkpoint_sweep` from checkpoint HDF5 metadata (`/meta/checkpoint_sweep`,
or root-level fallback). Returns presence/value/error fields.
"""
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

"""
    read_log_signals(log_path)

Extract fallback convergence/run signals from `run.log`:
early-stop marker, results-written marker, error signatures, and latest checkpoint sweep.
"""
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

"""
    assess_run(campaign_name, run_dir)

Compute one unified run assessment row used by convergence and aggregation pipelines.

Includes:
- run identity (`campaign_name`, `run_id`)
- key params (`seed_initial_state`, `t_a`, ..., `mu_b`)
- `convergence_status` (`converged|not_converged|unknown`)
- `run_status` (`ok|failed|missing|unknown`)
- sweep/diagnostic flags and `evidence`
"""
function assess_run(
    campaign_name::AbstractString,
    run_dir::AbstractString
)
    run_dir_abs = abspath(run_dir)
    run_id = basename(normpath(run_dir))
    results_path = joinpath(run_dir_abs, "results.h5")
    checkpoint_path = joinpath(run_dir_abs, "dmrg_state_checkpoint.h5")
    log_path = joinpath(run_dir_abs, "run.log")

    hparams = read_hamiltonian_params(run_dir_abs, results_path)
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

    return (
        campaign_name=String(campaign_name),
        run_id=String(run_id),
        seed_initial_state=hparams.seed_initial_state,
        t_a=hparams.t_a,
        t_b=hparams.t_b,
        U_a=hparams.U_a,
        U_b=hparams.U_b,
        U_ab=hparams.U_ab,
        mu_a=hparams.mu_a,
        mu_b=hparams.mu_b,
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
