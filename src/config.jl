# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

"""
    load_params(path::AbstractString)

Load parameters from a YAML file using YAML.jl.
"""
function load_params(path::AbstractString)
    return YAML.load_file(path)
end

"""
    dict_to_namedtuple(d::AbstractDict)

Convert a Dict with string/symbol keys to a NamedTuple suitable for keyword splatting.
"""
function dict_to_namedtuple(d::AbstractDict)
    return (; (Symbol(k) => v for (k, v) in d)...)
end

"""
    get_section(cfg::AbstractDict, name::AbstractString)

Fetch a section from a config dict with a helpful error if missing.
"""
function get_section(cfg::AbstractDict, name::AbstractString)
    haskey(cfg, name) || throw(ArgumentError("missing section '$name' in config"))
    return cfg[name]
end

"""
    merge_sections(cfg::AbstractDict, names::Vector{String})

Merge multiple config sections into a single NamedTuple. Later sections win on key conflicts.
Missing sections are skipped.
"""
function merge_sections(cfg::AbstractDict, names::Vector{String})
    nt = (;)
    for name in names
        if haskey(cfg, name)
            nt = merge(nt, dict_to_namedtuple(cfg[name]))
        end
    end
    return nt
end

@inline function _cfg_get(cfg::AbstractDict, key::String, default=nothing)
    if haskey(cfg, key)
        return cfg[key]
    end
    sk = Symbol(key)
    if haskey(cfg, sk)
        return cfg[sk]
    end
    return default
end

"""
    default_observables_path()

Return the default observables config path (`configs/observables.yaml`).
"""
function default_observables_path()
    return normpath(joinpath(@__DIR__, "..", "configs", "observables.yaml"))
end

"""
    with_observables_config(base_cfg; observables_path=nothing)

Return a NamedTuple with:
- `cfg`: copy of `base_cfg` with `"observables"` overridden from observables file when available.
- `observables_path`: resolved observables file path.
- `observables_loaded`: whether observables file was found and loaded.

If the observables file has a top-level `observables` section, that section is used.
Otherwise, the entire observables file is treated as the observables section.
If the observables file is missing, `base_cfg` is returned unchanged.
"""
function with_observables_config(base_cfg::AbstractDict; observables_path::Union{Nothing,AbstractString}=nothing)
    cfg = Dict{String,Any}(String(k) => v for (k, v) in base_cfg)
    obs_path = observables_path === nothing ?
               get(ENV, "OBSERVABLES", default_observables_path()) :
               String(observables_path)
    loaded = false

    if isfile(obs_path)
        raw_obs_cfg = load_params(obs_path)
        if raw_obs_cfg isa AbstractDict
            section = _cfg_get(raw_obs_cfg, "observables", nothing)
            if section isa AbstractDict
                cfg["observables"] = Dict{String,Any}(String(k) => v for (k, v) in section)
            else
                cfg["observables"] = Dict{String,Any}(String(k) => v for (k, v) in raw_obs_cfg)
            end
            loaded = true
        end
    end

    return (; cfg, observables_path=obs_path, observables_loaded=loaded)
end
