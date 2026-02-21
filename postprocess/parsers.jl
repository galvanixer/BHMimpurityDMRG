# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

using HDF5
using YAML
using Statistics

const DEFAULT_RESULTS_SCHEMA_ID = "bhmimpuritydmrg.results"
const DEFAULT_RESULTS_SCHEMA_VERSION = "1.0.0"

function normalize_yaml(x)
    if x isa AbstractDict
        return Dict{String,Any}(String(k) => normalize_yaml(v) for (k, v) in x)
    elseif x isa AbstractVector
        return [normalize_yaml(v) for v in x]
    else
        return x
    end
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

@inline function as_string_or_nothing(x)
    x === nothing && return nothing
    try
        return String(x)
    catch
        return string(x)
    end
end

@inline function read_if_exists(parent, name::AbstractString; default=nothing)
    haskey(parent, name) || return default
    return read(parent[name])
end

function read_group_recursive(g)
    out = Dict{String,Any}()
    for key in keys(g)
        obj = g[key]
        skey = String(key)
        if obj isa HDF5.Group
            out[skey] = read_group_recursive(obj)
        else
            try
                out[skey] = read(obj)
            catch
                # Keep parser resilient to unusual dataset types.
            end
        end
    end
    return out
end

@inline function profile_bool(profile::AbstractDict, key::AbstractString, default::Bool)
    return parse_bool(get(profile, key, default), default)
end

@inline function nested_get(d::AbstractDict, path::Vector{String}, default=nothing)
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

@inline function safe_mean(x)
    if x isa AbstractArray && !isempty(x)
        return Float64(mean(x))
    end
    return nothing
end

function maxdim_max_value(x)
    if x === nothing
        return nothing
    elseif x isa Number
        return x
    elseif x isa AbstractVector
        vals = Number[]
        for v in x
            mv = maxdim_max_value(v)
            if mv isa Number
                push!(vals, mv)
            end
        end
        if isempty(vals)
            return nothing
        end
        has_float = any(v -> v isa AbstractFloat, vals)
        if has_float
            return maximum(Float64.(vals))
        end
        return maximum(Int64.(vals))
    elseif x isa AbstractDict
        if haskey(x, "max")
            return maxdim_max_value(x["max"])
        elseif haskey(x, :max)
            return maxdim_max_value(x[:max])
        elseif haskey(x, "values")
            return maxdim_max_value(x["values"])
        elseif haskey(x, :values)
            return maxdim_max_value(x[:values])
        else
            return nothing
        end
    else
        return nothing
    end
end

function extract_key_params(meta::AbstractDict)
    params_yaml = get(meta, "params_yaml", nothing)
    params_yaml isa AbstractString || return Dict{String,Any}()
    cfg = try
        normalize_yaml(YAML.load(params_yaml))
    catch
        Dict{String,Any}()
    end

    out = Dict{String,Any}()
    # Frequently-used campaign knobs for quick summary filtering.
    out["cfg_lattice_L"] = nested_get(cfg, ["lattice", "L"], nothing)
    out["cfg_lattice_periodic"] = nested_get(cfg, ["lattice", "periodic"], nothing)
    out["cfg_initial_Na_total"] = nested_get(cfg, ["initial_state", "Na_total"], nothing)
    out["cfg_initial_Nb_total"] = nested_get(cfg, ["initial_state", "Nb_total"], nothing)
    impurity_distribution = nested_get(cfg, ["initial_state", "impurity_distribution"], nothing)
    out["cfg_initial_impurity_distribution"] = impurity_distribution
    out["cfg_initial_seed"] = begin
        if impurity_distribution === nothing
            nothing
        elseif lowercase(string(impurity_distribution)) == "random"
            nested_get(cfg, ["initial_state", "seed"], nothing)
        else
            nothing
        end
    end
    out["cfg_local_nmax_a"] = nested_get(cfg, ["local_hilbert", "nmax_a"], nothing)
    out["cfg_local_nmax_b"] = nested_get(cfg, ["local_hilbert", "nmax_b"], nothing)
    out["cfg_t_a"] = nested_get(cfg, ["hamiltonian", "t_a"], nothing)
    out["cfg_t_b"] = nested_get(cfg, ["hamiltonian", "t_b"], nothing)
    out["cfg_U_a"] = nested_get(cfg, ["hamiltonian", "U_a"], nothing)
    out["cfg_U_b"] = nested_get(cfg, ["hamiltonian", "U_b"], nothing)
    out["cfg_U_ab"] = nested_get(cfg, ["hamiltonian", "U_ab"], nothing)
    out["cfg_mu_a"] = nested_get(cfg, ["hamiltonian", "mu_a"], nothing)
    out["cfg_mu_b"] = nested_get(cfg, ["hamiltonian", "mu_b"], nothing)
    out["cfg_dmrg_nsweeps"] = nested_get(cfg, ["dmrg", "nsweeps"], nothing)
    out["cfg_dmrg_maxdim_max"] = maxdim_max_value(nested_get(cfg, ["dmrg", "maxdim"], nothing))
    out["cfg_meta_run_name"] = nested_get(cfg, ["meta", "run_name"], nothing)
    return out
end

function summarize_observables(g_obs)
    summary = Dict{String,Any}()
    summary["has_observables"] = true

    if haskey(g_obs, "energy")
        g_energy = g_obs["energy"]
        summary["E0"] = read_if_exists(g_energy, "E0"; default=nothing)
    else
        summary["E0"] = nothing
    end

    if haskey(g_obs, "totals")
        g_tot = g_obs["totals"]
        summary["Na"] = read_if_exists(g_tot, "Na"; default=nothing)
        summary["Nb"] = read_if_exists(g_tot, "Nb"; default=nothing)
    else
        summary["Na"] = nothing
        summary["Nb"] = nothing
    end

    if haskey(g_obs, "densities")
        g_den = g_obs["densities"]
        na = read_if_exists(g_den, "na"; default=nothing)
        nb = read_if_exists(g_den, "nb"; default=nothing)
        summary["L_from_density"] = na isa AbstractArray ? length(na) : nothing
        summary["na_mean"] = safe_mean(na)
        summary["nb_mean"] = safe_mean(nb)
    else
        summary["L_from_density"] = nothing
        summary["na_mean"] = nothing
        summary["nb_mean"] = nothing
    end

    summary["has_density_density"] = haskey(g_obs, "density_density")
    summary["has_structure_factor"] = haskey(g_obs, "structure_factor")
    summary["has_triple_corr"] = haskey(g_obs, "triple_corr")
    summary["has_sampled_configs"] = haskey(g_obs, "sampled_configs")

    return summary
end

function collect_observables_payload(g_obs, profile::AbstractDict)
    payload = Dict{String,Any}()

    # Always keep cheap scalar groups.
    if haskey(g_obs, "energy")
        payload["energy"] = read_group_recursive(g_obs["energy"])
    end
    if haskey(g_obs, "totals")
        payload["totals"] = read_group_recursive(g_obs["totals"])
    end

    include_arrays = profile_bool(profile, "include_arrays", false)
    include_density_vectors = include_arrays && profile_bool(profile, "include_density_vectors", false)
    include_density_density = include_arrays && profile_bool(profile, "include_density_density", false)
    include_structure_factor = include_arrays && profile_bool(profile, "include_structure_factor", false)
    include_triple_corr = include_arrays && profile_bool(profile, "include_triple_corr", false)
    include_sampled_configs = include_arrays && profile_bool(profile, "include_sampled_configs", false)

    if include_density_vectors && haskey(g_obs, "densities")
        payload["densities"] = read_group_recursive(g_obs["densities"])
    end
    if include_density_density && haskey(g_obs, "density_density")
        payload["density_density"] = read_group_recursive(g_obs["density_density"])
    end
    if include_structure_factor && haskey(g_obs, "structure_factor")
        payload["structure_factor"] = read_group_recursive(g_obs["structure_factor"])
    end
    if include_triple_corr && haskey(g_obs, "triple_corr")
        payload["triple_corr"] = read_group_recursive(g_obs["triple_corr"])
    end
    if include_sampled_configs && haskey(g_obs, "sampled_configs")
        payload["sampled_configs"] = read_group_recursive(g_obs["sampled_configs"])
    end

    return payload
end

function parse_common_results(
    results_path::AbstractString,
    profile::AbstractDict;
    schema_id::AbstractString,
    schema_version::AbstractString,
    parser_name::AbstractString,
    writer=nothing
)
    issues = String[]
    meta = Dict{String,Any}()
    observables = Dict{String,Any}()
    summary = Dict{String,Any}()

    HDF5.h5open(results_path, "r") do f
        if haskey(f, "meta")
            g_meta = f["meta"]
            for key in keys(g_meta)
                skey = String(key)
                obj = g_meta[key]
                if obj isa HDF5.Group
                    continue
                end
                try
                    meta[skey] = read(obj)
                catch
                    push!(issues, "failed_to_read_meta_dataset:$skey")
                end
            end
            attrs = HDF5.attributes(g_meta)
            if haskey(attrs, "run_by")
                try
                    meta["run_by"] = read(attrs["run_by"])
                catch
                    push!(issues, "failed_to_read_meta_attr:run_by")
                end
            end
        else
            push!(issues, "missing_group:/meta")
        end

        if haskey(f, "observables")
            g_obs = f["observables"]
            summary = summarize_observables(g_obs)
            observables = collect_observables_payload(g_obs, profile)
        else
            summary = Dict{String,Any}(
                "has_observables" => false,
                "E0" => nothing,
                "Na" => nothing,
                "Nb" => nothing,
                "na_mean" => nothing,
                "nb_mean" => nothing
            )
            push!(issues, "missing_group:/observables")
        end
    end

    merge!(summary, extract_key_params(meta))
    summary["schema_id"] = schema_id
    summary["schema_version"] = schema_version
    summary["status"] = "ok"

    return Dict{String,Any}(
        "schema" => Dict{String,Any}(
            "results_schema_id" => schema_id,
            "results_schema_version" => schema_version,
            "results_writer" => writer,
            "parser_name" => parser_name
        ),
        "meta" => meta,
        "summary" => summary,
        "observables" => observables,
        "issues" => issues
    )
end

function parse_v1_results(results_path::AbstractString, profile::AbstractDict; schema_id::AbstractString, schema_version::AbstractString, writer=nothing)
    return parse_common_results(
        results_path,
        profile;
        schema_id=schema_id,
        schema_version=schema_version,
        parser_name="v1",
        writer=writer
    )
end

function parse_legacy_results(results_path::AbstractString, profile::AbstractDict)
    return parse_common_results(
        results_path,
        profile;
        schema_id="legacy",
        schema_version="0",
        parser_name="legacy_fallback",
        writer=nothing
    )
end

function detect_results_schema(results_path::AbstractString)
    return HDF5.h5open(results_path, "r") do f
        haskey(f, "meta") || return (nothing, nothing, nothing)
        g_meta = f["meta"]
        schema_id = as_string_or_nothing(read_if_exists(g_meta, "results_schema_id"; default=nothing))
        schema_version = as_string_or_nothing(read_if_exists(g_meta, "results_schema_version"; default=nothing))
        writer = as_string_or_nothing(read_if_exists(g_meta, "results_writer"; default=nothing))
        return (schema_id, schema_version, writer)
    end
end

function parse_results_file(results_path::AbstractString; profile::AbstractDict=Dict{String,Any}())
    schema_id, schema_version, writer = detect_results_schema(results_path)

    if schema_id === nothing || schema_version === nothing
        return parse_legacy_results(results_path, profile)
    end

    if schema_id == DEFAULT_RESULTS_SCHEMA_ID && schema_version == DEFAULT_RESULTS_SCHEMA_VERSION
        return parse_v1_results(
            results_path,
            profile;
            schema_id=schema_id,
            schema_version=schema_version,
            writer=writer
        )
    end

    parsed = parse_common_results(
        results_path,
        profile;
        schema_id=schema_id,
        schema_version=schema_version,
        parser_name="unknown_schema_fallback",
        writer=writer
    )
    push!(parsed["issues"], "unknown_schema:$(schema_id)@$(schema_version)")
    return parsed
end
