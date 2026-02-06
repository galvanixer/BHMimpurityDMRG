# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

"""
    load_params(path::AbstractString)

Load parameters from a YAML file using YAML.jl.

Requires `YAML.jl` to be installed:
    import Pkg; Pkg.add("YAML")
"""
function load_params(path::AbstractString)
    @eval using YAML
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
