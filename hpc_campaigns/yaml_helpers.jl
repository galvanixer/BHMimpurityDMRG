# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

module YAMLHelpers

using YAML
using OrderedCollections: OrderedDict

export as_dict
export load_ordered_yaml
export order_like
export set_nested!
export sweep_assignments
export write_yaml_ordered

function normalize_yaml(x)
    if x isa AbstractDict
        d = OrderedDict{String,Any}()
        for (k, v) in x
            d[String(k)] = normalize_yaml(v)
        end
        return d
    elseif x isa AbstractVector
        return [normalize_yaml(v) for v in x]
    else
        return x
    end
end

function as_dict(x)
    x isa AbstractDict || throw(ArgumentError("expected dictionary, got $(typeof(x))"))
    return normalize_yaml(x)
end

function load_ordered_yaml(path::AbstractString)
    return as_dict(YAML.load_file(path; dicttype=OrderedDict{Any,Any}))
end

function order_like(reference::AbstractDict, data::AbstractDict)
    out = OrderedDict{String,Any}()

    # Keep keys that appear in reference first, recursively preserving sub-order.
    for (k_ref, v_ref) in reference
        key = String(k_ref)
        if haskey(data, key)
            v = data[key]
            if v isa AbstractDict && v_ref isa AbstractDict
                out[key] = order_like(v_ref, v)
            else
                out[key] = v
            end
        end
    end

    # Append keys that are not in reference (for newly added fields).
    for (k, v) in data
        key = String(k)
        if !haskey(out, key)
            out[key] = v
        end
    end

    return out
end

function set_nested!(d::AbstractDict, dotted_key::AbstractString, value)
    parts = split(String(dotted_key), ".")
    isempty(parts) && throw(ArgumentError("empty key"))
    cur = d
    for p in parts[1:(end - 1)]
        if !haskey(cur, p) || !(cur[p] isa AbstractDict)
            cur[p] = OrderedDict{String,Any}()
        end
        cur = cur[p]
    end
    cur[parts[end]] = value
    return d
end

as_vector(v) = v isa AbstractVector ? collect(v) : [v]

function sweep_assignments(sweep::AbstractDict)
    sweepd = OrderedDict{String,Any}(String(k) => v for (k, v) in sweep)
    keys_sorted = sort(collect(keys(sweepd)))
    if isempty(keys_sorted)
        return [OrderedDict{String,Any}()]
    end
    values = [as_vector(sweepd[k]) for k in keys_sorted]
    assigns = OrderedDict{String,Any}[]
    for tuple_vals in Iterators.product(values...)
        a = OrderedDict{String,Any}()
        for (k, v) in zip(keys_sorted, tuple_vals)
            a[k] = v
        end
        push!(assigns, a)
    end
    return assigns
end

function _inline_block_array(text::AbstractString, key::AbstractString)
    s = String(text)
    key_quoted = "\\Q" * String(key) * "\\E"
    key_line = Regex("^([ \\t]*)$(key_quoted):[ \\t]*\$")
    item_line = Regex("^([ \\t]*)  -[ \\t]*(.*)\$")

    lines = split(s, '\n'; keepempty=true)
    out = String[]
    i = 1
    while i <= length(lines)
        m_key = match(key_line, lines[i])
        if m_key === nothing
            push!(out, lines[i])
            i += 1
            continue
        end

        indent = m_key.captures[1]
        vals = String[]
        j = i + 1
        while j <= length(lines)
            m_item = match(item_line, lines[j])
            if m_item === nothing || m_item.captures[1] != indent
                break
            end
            push!(vals, strip(m_item.captures[2]))
            j += 1
        end

        if isempty(vals)
            push!(out, lines[i])
            i += 1
        else
            push!(out, indent * key * ": [" * join(vals, ", ") * "]")
            i = j
        end
    end

    return join(out, "\n")
end

function _compact_block_list_of_lists(text::AbstractString, key::AbstractString)
    s = String(text)
    key_quoted = "\\Q" * String(key) * "\\E"
    key_line = Regex("^([ \\t]*)$(key_quoted):[ \\t]*\$")
    row_start = Regex("^([ \\t]*)  -[ \\t]*\$")
    row_item = Regex("^([ \\t]*)    -[ \\t]*(.*)\$")

    lines = split(s, '\n'; keepempty=true)
    out = String[]
    i = 1
    while i <= length(lines)
        m_key = match(key_line, lines[i])
        if m_key === nothing
            push!(out, lines[i])
            i += 1
            continue
        end

        indent = m_key.captures[1]
        rows = Vector{Vector{String}}()
        j = i + 1
        while j <= length(lines)
            m_row_start = match(row_start, lines[j])
            if m_row_start === nothing || m_row_start.captures[1] != indent
                break
            end
            j += 1
            vals = String[]
            while j <= length(lines)
                m_row_item = match(row_item, lines[j])
                if m_row_item === nothing || m_row_item.captures[1] != indent
                    break
                end
                push!(vals, strip(m_row_item.captures[2]))
                j += 1
            end
            isempty(vals) && break
            push!(rows, vals)
        end

        if isempty(rows)
            push!(out, lines[i])
            i += 1
        else
            push!(out, lines[i])
            for row in rows
                push!(out, indent * "  - [" * join(row, ", ") * "]")
            end
            i = j
        end
    end

    return join(out, "\n")
end

function rewrite_yaml_layout!(
    params_path::AbstractString;
    inline_keys::Vector{String}=String["maxdim"],
    compact_row_keys::Vector{String}=String["pairs"]
)
    txt = read(params_path, String)
    new_txt = txt
    for key in inline_keys
        new_txt = _inline_block_array(new_txt, key)
    end
    for key in compact_row_keys
        new_txt = _compact_block_list_of_lists(new_txt, key)
    end
    if new_txt != txt
        write(params_path, new_txt)
    end
    return nothing
end

function write_yaml_ordered(
    path::AbstractString,
    data::AbstractDict;
    reference::Union{Nothing,AbstractDict}=nothing,
    inline_keys::Vector{String}=String["maxdim"],
    compact_row_keys::Vector{String}=String["pairs"]
)
    to_write = isnothing(reference) ? data : order_like(reference, data)
    YAML.write_file(path, to_write)
    rewrite_yaml_layout!(path; inline_keys=inline_keys, compact_row_keys=compact_row_keys)
    return nothing
end

end
