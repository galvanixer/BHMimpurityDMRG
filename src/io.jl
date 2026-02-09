# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

using Dates
using SHA

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

function write_meta!(g_meta; params_path=nothing, params_text=nothing)
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
    save_state(path::AbstractString, psi::MPS; energy=nothing, sites=siteinds(psi),
               params_path=nothing, params_text=nothing, na=nothing, nb=nothing)

Save the ground state `psi` (and optionally `energy`, `sites`, YAML parameters, and
site densities `na`, `nb`)
to an HDF5 file.
"""
function save_state(path::AbstractString, psi::MPS; energy=nothing, sites=siteinds(psi),
    params_path=nothing, params_text=nothing, na=nothing, nb=nothing)
    HDF5.h5open(path, "w") do f
        g_state = HDF5.create_group(f, "state")
        write(g_state, "psi", psi)
        write(g_state, "sites", sites)
        if energy !== nothing
            write(g_state, "energy", energy)
        end

        g_meta = HDF5.create_group(f, "meta")
        write_meta!(g_meta; params_path=params_path, params_text=params_text)

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
    end
    return nothing
end

"""
    load_state(path::AbstractString)

Load a saved MPS ground state from an HDF5 file.

Returns a NamedTuple `(psi, sites, energy, params_yaml, params_sha256, na, nb)` where
`energy`, `params_yaml`, `params_sha256`, `na`, and `nb` may be `nothing` if they were not stored.
"""
function load_state(path::AbstractString)
    HDF5.h5open(path, "r") do f
        if haskey(f, "state")
            g_state = f["state"]
            psi = read(g_state, "psi", MPS)
            sites = read(g_state, "sites", Vector{Index})
            energy = haskey(g_state, "energy") ? read(g_state, "energy") : nothing
        else
            # Backward compatibility: root-level datasets
            psi = read(f, "psi", MPS)
            sites = read(f, "sites", Vector{Index})
            energy = haskey(f, "energy") ? read(f, "energy") : nothing
        end

        params_yaml = nothing
        params_sha256 = nothing
        if haskey(f, "meta")
            g_meta = f["meta"]
            if haskey(g_meta, "params_yaml")
                params_yaml = read(g_meta, "params_yaml")
            end
            if haskey(g_meta, "params_sha256")
                params_sha256 = read(g_meta, "params_sha256")
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
        return (; psi, sites, energy, params_yaml, params_sha256, na, nb)
    end
end
