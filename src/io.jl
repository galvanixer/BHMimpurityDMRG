# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

using Dates
using SHA

const RESULTS_SCHEMA_ID = "bhmimpuritydmrg.results"
const RESULTS_SCHEMA_VERSION = "1.0.0"
const RESULTS_WRITER = "BHMimpurityDMRG"

function ensure_group(parent, name::AbstractString)
    return haskey(parent, name) ? parent[name] : HDF5.create_group(parent, name)
end

function write_or_replace(parent, name::AbstractString, data)
    if haskey(parent, name)
        HDF5.delete_object(parent, name)
    end
    write(parent, name, data)
    return nothing
end

function write_meta!(g_meta; params_path=nothing, params_text=nothing, state_path=nothing, state_params_sha256=nothing)
    if params_text === nothing && params_path !== nothing
        params_text = read(params_path, String)
    end
    if params_text !== nothing
        write_or_replace(g_meta, "params_yaml", params_text)
        write_or_replace(g_meta, "params_sha256", bytes2hex(SHA.sha256(params_text)))
    end
    if params_path !== nothing
        write_or_replace(g_meta, "params_path", abspath(params_path))
    end
    if state_path !== nothing
        write_or_replace(g_meta, "state_path", abspath(state_path))
    end
    if state_params_sha256 !== nothing
        write_or_replace(g_meta, "state_params_sha256", String(state_params_sha256))
    end
    if !haskey(g_meta, "created_at")
        write(g_meta, "created_at", string(Dates.now()))
    end

    # Authors / run_by metadata (best effort)
    if params_text !== nothing
        try
            cfg = YAML.load(params_text)
            if isa(cfg, AbstractDict)
                meta = get(cfg, "meta", Dict{String,Any}())
                authors = get(meta, "authors", nothing)
                run_by = get(meta, "run_by", nothing)
                if authors === nothing
                    authors = AUTHORS
                end
                if authors isa AbstractVector && !isempty(authors)
                    write_or_replace(g_meta, "authors", String.(authors))
                end
                if run_by === nothing || isempty(String(run_by))
                    run_by = get(ENV, "USER", get(ENV, "USERNAME", "unknown"))
                end
                if run_by !== nothing && !isempty(String(run_by))
                    HDF5.attributes(g_meta)["run_by"] = String(run_by)
                end
            end
        catch
            # Best-effort metadata only; ignore parse errors
        end
    end
    return nothing
end

"""
    write_results_schema!(g_meta; schema_id=RESULTS_SCHEMA_ID,
                          schema_version=RESULTS_SCHEMA_VERSION,
                          writer=RESULTS_WRITER)

Write canonical schema metadata for a `results.h5` file under `/meta`.
"""
function write_results_schema!(
    g_meta;
    schema_id::AbstractString=RESULTS_SCHEMA_ID,
    schema_version::AbstractString=RESULTS_SCHEMA_VERSION,
    writer::AbstractString=RESULTS_WRITER
)
    write_or_replace(g_meta, "results_schema_id", String(schema_id))
    write_or_replace(g_meta, "results_schema_version", String(schema_version))
    write_or_replace(g_meta, "results_writer", String(writer))
    return nothing
end

"""
    write_dmrg_diagnostics!(f, diag)

Write DMRG diagnostics to `/diagnostics/dmrg` in an open HDF5 file `f`.
The per-sweep trace is stored as a single compound dataset `sweep_trace`.
"""
function write_dmrg_diagnostics!(f, diag)
    g_diag = ensure_group(f, "diagnostics")
    if haskey(g_diag, "dmrg")
        HDF5.delete_object(g_diag, "dmrg")
    end
    g_dmrg = HDF5.create_group(g_diag, "dmrg")

    sweep_trace = hasproperty(diag, :sweep_trace) ? diag.sweep_trace : DMRGSweepTraceRow[]
    sweep_trace_v = if sweep_trace isa AbstractVector{DMRGSweepTraceRow}
        collect(sweep_trace)
    elseif sweep_trace isa AbstractVector
        DMRGSweepTraceRow[]
    else
        DMRGSweepTraceRow[]
    end
    write_or_replace(g_dmrg, "sweep_trace", sweep_trace_v)

    checkpoint_sweeps = hasproperty(diag, :checkpoint_sweeps) ? Int64.(collect(diag.checkpoint_sweeps)) : Int64[]
    write_or_replace(g_dmrg, "checkpoint_sweeps", checkpoint_sweeps)
    write_or_replace(g_dmrg, "sweeps_completed", Int(hasproperty(diag, :sweeps_completed) ? diag.sweeps_completed : length(sweep_trace_v)))
    write_or_replace(g_dmrg, "converged", Bool(hasproperty(diag, :converged) ? diag.converged : false))
    write_or_replace(g_dmrg, "early_stop_triggered", Bool(hasproperty(diag, :early_stop_triggered) ? diag.early_stop_triggered : false))
    write_or_replace(g_dmrg, "energy_tol", Float64(hasproperty(diag, :energy_tol) ? diag.energy_tol : 0.0))
    write_or_replace(g_dmrg, "trunc_tol", Float64(hasproperty(diag, :trunc_tol) ? diag.trunc_tol : 0.0))
    write_or_replace(g_dmrg, "patience", Int(hasproperty(diag, :patience) ? diag.patience : 1))
    write_or_replace(g_dmrg, "min_sweeps", Int(hasproperty(diag, :min_sweeps) ? diag.min_sweeps : 2))
    write_or_replace(g_dmrg, "resume_mode", String(hasproperty(diag, :resume_mode) ? diag.resume_mode : "unknown"))
    write_or_replace(g_dmrg, "checkpoint_sweep_start", Int(hasproperty(diag, :checkpoint_sweep_start) ? diag.checkpoint_sweep_start : 0))
    return nothing
end

@inline function _diag_row_value(row, name::Symbol, default)
    if hasproperty(row, name)
        return getproperty(row, name)
    elseif row isa AbstractDict
        if haskey(row, name)
            return row[name]
        end
        sname = String(name)
        if haskey(row, sname)
            return row[sname]
        end
    end
    return default
end

function _normalize_sweep_trace(rows)
    rows isa AbstractVector{DMRGSweepTraceRow} && return collect(rows)
    rows isa AbstractVector || return DMRGSweepTraceRow[]

    out = DMRGSweepTraceRow[]
    sizehint!(out, length(rows))
    for row in rows
        push!(out, DMRGSweepTraceRow(
            Float64(_diag_row_value(row, :energy, NaN)),
            Float64(_diag_row_value(row, :delta_energy, NaN)),
            Float64(_diag_row_value(row, :max_truncerr, NaN)),
            Int64(_diag_row_value(row, :maxdim_used, 0)),
            Float64(_diag_row_value(row, :walltime_sec, NaN))
        ))
    end
    return out
end

function load_dmrg_diagnostics(f)
    haskey(f, "diagnostics") || return nothing
    g_diag = f["diagnostics"]
    haskey(g_diag, "dmrg") || return nothing
    g_dmrg = g_diag["dmrg"]

    sweep_trace = haskey(g_dmrg, "sweep_trace") ?
                  _normalize_sweep_trace(read(g_dmrg, "sweep_trace")) :
                  DMRGSweepTraceRow[]
    checkpoint_sweeps = haskey(g_dmrg, "checkpoint_sweeps") ?
                        Int64.(collect(read(g_dmrg, "checkpoint_sweeps"))) :
                        Int64[]

    return (
        sweep_trace=sweep_trace,
        checkpoint_sweeps=checkpoint_sweeps,
        sweeps_completed=haskey(g_dmrg, "sweeps_completed") ?
                         Int(read(g_dmrg, "sweeps_completed")) :
                         length(sweep_trace),
        converged=haskey(g_dmrg, "converged") ?
                  Bool(read(g_dmrg, "converged")) :
                  false,
        early_stop_triggered=haskey(g_dmrg, "early_stop_triggered") ?
                             Bool(read(g_dmrg, "early_stop_triggered")) :
                             false,
        energy_tol=haskey(g_dmrg, "energy_tol") ?
                   Float64(read(g_dmrg, "energy_tol")) :
                   0.0,
        trunc_tol=haskey(g_dmrg, "trunc_tol") ?
                  Float64(read(g_dmrg, "trunc_tol")) :
                  0.0,
        patience=haskey(g_dmrg, "patience") ?
                 Int(read(g_dmrg, "patience")) :
                 1,
        min_sweeps=haskey(g_dmrg, "min_sweeps") ?
                   Int(read(g_dmrg, "min_sweeps")) :
                   2,
        resume_mode=haskey(g_dmrg, "resume_mode") ?
                    String(read(g_dmrg, "resume_mode")) :
                    "unknown",
        checkpoint_sweep_start=haskey(g_dmrg, "checkpoint_sweep_start") ?
                               Int(read(g_dmrg, "checkpoint_sweep_start")) :
                               0
    )
end

"""
    save_state(path::AbstractString, psi::MPS; energy=nothing, sites=siteinds(psi),
               energy_variance=nothing, params_path=nothing, params_text=nothing, na=nothing, nb=nothing,
               init_na=nothing, init_nb=nothing, checkpoint_sweep=nothing, dmrg_diagnostics=nothing)

Save the ground state `psi` (and optionally `energy`, `sites`, YAML parameters, and
site densities `na`, `nb`)
to an HDF5 file.
"""
function save_state(path::AbstractString, psi::MPS; energy=nothing, sites=siteinds(psi),
    energy_variance=nothing, params_path=nothing, params_text=nothing, na=nothing, nb=nothing,
    init_na=nothing, init_nb=nothing, checkpoint_sweep=nothing, dmrg_diagnostics=nothing)
    HDF5.h5open(path, "w") do f
        g_state = HDF5.create_group(f, "state")
        write(g_state, "psi", psi)
        write(g_state, "sites", sites)
        if energy !== nothing
            write(g_state, "energy", energy)
        end
        if energy_variance !== nothing
            write(g_state, "energy_variance", energy_variance)
        end

        g_meta = HDF5.create_group(f, "meta")
        write_meta!(g_meta; params_path=params_path, params_text=params_text)
        if checkpoint_sweep !== nothing
            write_or_replace(g_meta, "checkpoint_sweep", Int(checkpoint_sweep))
        end

        if na !== nothing || nb !== nothing
            g_obs = HDF5.create_group(f, "observables")
            g_den = HDF5.create_group(g_obs, "densities")
            if na !== nothing
                write(g_den, "na", na)
            end
            if nb !== nothing
                write(g_den, "nb", nb)
            end
        end

        if init_na !== nothing || init_nb !== nothing
            g_init = HDF5.create_group(f, "initial_state")
            g_init_obs = HDF5.create_group(g_init, "observables")
            if init_na !== nothing
                write(g_init_obs, "na", init_na)
            end
            if init_nb !== nothing
                write(g_init_obs, "nb", init_nb)
            end
        end

        if dmrg_diagnostics !== nothing
            write_dmrg_diagnostics!(f, dmrg_diagnostics)
        end
    end
    return nothing
end

"""
    load_state(path::AbstractString)

Load a saved MPS ground state from an HDF5 file.

    Returns a NamedTuple `(psi, sites, energy, energy_variance, params_yaml, params_sha256, na, nb, init_na, init_nb, checkpoint_sweep, dmrg_diagnostics)` where
    optional fields may be `nothing` if they were not stored.
"""
function load_state(path::AbstractString)
    HDF5.h5open(path, "r") do f
        if haskey(f, "state")
            g_state = f["state"]
            psi = read(g_state, "psi", MPS)
            sites = read(g_state, "sites", Vector{Index})
            energy = haskey(g_state, "energy") ? read(g_state, "energy") : nothing
            energy_variance = haskey(g_state, "energy_variance") ? read(g_state, "energy_variance") : nothing
        else
            # Backward compatibility: root-level datasets
            psi = read(f, "psi", MPS)
            sites = read(f, "sites", Vector{Index})
            energy = haskey(f, "energy") ? read(f, "energy") : nothing
            energy_variance = haskey(f, "energy_variance") ? read(f, "energy_variance") : nothing
        end

        params_yaml = nothing
        params_sha256 = nothing
        checkpoint_sweep = nothing
        if haskey(f, "meta")
            g_meta = f["meta"]
            if haskey(g_meta, "params_yaml")
                params_yaml = read(g_meta, "params_yaml")
            end
            if haskey(g_meta, "params_sha256")
                params_sha256 = read(g_meta, "params_sha256")
            end
            if haskey(g_meta, "checkpoint_sweep")
                checkpoint_sweep = Int(read(g_meta, "checkpoint_sweep"))
            end
        end
        if params_yaml === nothing && haskey(f, "params_yaml")
            params_yaml = read(f, "params_yaml")
        end

        na = nothing
        nb = nothing
        if haskey(f, "observables") && haskey(f["observables"], "densities")
            g_den = f["observables"]["densities"]
            na = haskey(g_den, "na") ? read(g_den, "na") : nothing
            nb = haskey(g_den, "nb") ? read(g_den, "nb") : nothing
        end

        init_na = nothing
        init_nb = nothing
        if haskey(f, "initial_state") && haskey(f["initial_state"], "observables")
            g_init_obs = f["initial_state"]["observables"]
            init_na = haskey(g_init_obs, "na") ? read(g_init_obs, "na") : nothing
            init_nb = haskey(g_init_obs, "nb") ? read(g_init_obs, "nb") : nothing
        end

        dmrg_diagnostics = load_dmrg_diagnostics(f)

        return (; psi, sites, energy, energy_variance, params_yaml, params_sha256, na, nb, init_na, init_nb, checkpoint_sweep, dmrg_diagnostics)
    end
end
