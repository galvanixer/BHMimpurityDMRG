# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

function normalize_yaml(x)
    if x isa AbstractDict
        return Dict{String,Any}(String(k) => normalize_yaml(v) for (k, v) in x)
    elseif x isa AbstractVector
        return [normalize_yaml(v) for v in x]
    else
        return x
    end
end

@inline function as_string_or_nothing(x)
    x === nothing && return nothing
    try
        return String(x)
    catch
        return string(x)
    end
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
