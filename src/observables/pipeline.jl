# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

# ----------------------------
# Centralized observables compute/write pipeline
# ----------------------------

const OBSERVABLES_SCHEMA_VERSION = "1.0.0"

"""
    _cfg_get(d, key, default=nothing)

Read a config value from `d` by trying both string and symbol keys.
Returns `default` when `key` is not present.
"""
@inline function _cfg_get(d::AbstractDict, key::String, default=nothing)
    if haskey(d, key)
        return d[key]
    end
    sk = Symbol(key)
    if haskey(d, sk)
        return d[sk]
    end
    return default
end

"""
    _cfg_section(cfg, key)

Return the subsection `cfg[key]` as an `AbstractDict`, or an empty
`Dict{String,Any}` if the section is missing or not a dictionary.
"""
@inline function _cfg_section(cfg::AbstractDict, key::String)
    sec = _cfg_get(cfg, key, Dict{String,Any}())
    return sec isa AbstractDict ? sec : Dict{String,Any}()
end

"""
    _parse_bool(x, default=false)

Parse common boolean-like values (`"true"`, `"false"`, `"1"`, `"0"`, etc.).
If parsing fails, return `default`.
"""
@inline function _parse_bool(x, default::Bool=false)
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

"""
    _parse_species(species, section_name)

Map species selector `"a"`, `"b"`, or `"both"` to operator names
`["Na"]`, `["Nb"]`, or `["Na","Nb"]`. Throws on invalid input.
"""
function _parse_species(species, section_name::String)
    s = lowercase(String(species))
    if s == "a"
        return ["Na"]
    elseif s == "b"
        return ["Nb"]
    elseif s == "both"
        return ["Na", "Nb"]
    end
    error("$section_name.species must be one of: \"a\", \"b\", \"both\"")
end

"""
    _resolve_triple_pairs(tc_cfg, L)

Resolve `(r,s)` pairs for triple-correlation evaluation from config.
Supports explicit `pairs`, `all_pairs`, and optional `rmax`/`smax`.
"""
function _resolve_triple_pairs(tc_cfg::AbstractDict, L::Int)
    pairs_spec = _cfg_get(tc_cfg, "pairs", nothing)
    all_pairs = _parse_bool(_cfg_get(tc_cfg, "all_pairs", false), false) ||
                (pairs_spec isa AbstractString && lowercase(pairs_spec) == "all")
    rmax = _cfg_get(tc_cfg, "rmax", nothing)
    smax = _cfg_get(tc_cfg, "smax", nothing)

    pairs = Vector{Tuple{Int,Int}}()
    if all_pairs || rmax !== nothing || smax !== nothing
        rmax_eff = rmax === nothing ? (L - 1) : min(Int(rmax), L - 1)
        smax_eff = smax === nothing ? (L - 1) : min(Int(smax), L - 1)
        for r in 0:rmax_eff, s in 0:smax_eff
            push!(pairs, (r, s))
        end
    else
        raw_pairs = pairs_spec === nothing ?
                    [[0, 0], [0, 1], [0, 2], [1, 1], [1, 2], [1, 3], [2, 3]] :
                    pairs_spec
        pairs = [(Int(p[1]), Int(p[2])) for p in raw_pairs]
    end
    return pairs
end

"""
    _pairs_matrix(pairs)

Convert a vector of `(r,s)` tuples to an `N x 2` integer matrix for HDF5 output.
Returns an empty `0 x 2` matrix when `pairs` is empty.
"""
@inline function _pairs_matrix(pairs::Vector{Tuple{Int,Int}})
    isempty(pairs) && return zeros(Int, 0, 2)
    return hcat([p[1] for p in pairs], [p[2] for p in pairs])
end

"""
    compute_observables(psi, sites; ...)

Centralized observable evaluation from a state `psi`.

Returns a NamedTuple with canonical sections:
`energy`, `densities`, `totals`, `density_density`, `structure_factor`, `triple_corr`.
"""
function compute_observables(
    psi::MPS,
    sites;
    energy=nothing,
    na=nothing,
    nb=nothing,
    cfg::AbstractDict=Dict{String,Any}(),
    periodic::Union{Bool,Nothing}=nothing,
    compute_density_density::Bool=true,
    compute_structure_factor::Bool=true,
    compute_triple_corr::Bool=false,
    progress::Bool=false
)
    obs_cfg = _cfg_section(cfg, "observables")
    lattice_cfg = _cfg_section(cfg, "lattice")
    periodic_eff = periodic === nothing ?
                   _parse_bool(_cfg_get(lattice_cfg, "periodic", true), true) :
                   periodic

    if progress
        println("Evaluating observable: energy")
    end

    if progress
        println("Evaluating observable: densities / totals")
    end
    na_v, nb_v = if na !== nothing && nb !== nothing
        na, nb
    else
        measure_densities(psi, sites)
    end
    Na, Nb = total_numbers(na_v, nb_v)

    dd_cfg = _cfg_section(obs_cfg, "density_density")
    sf_cfg = _cfg_section(obs_cfg, "structure_factor")
    tc_cfg = _cfg_section(obs_cfg, "triple_corr")

    dd_species = compute_density_density ?
                 _parse_species(_cfg_get(dd_cfg, "species", "both"), "observables.density_density") :
                 String[]
    sf_species = compute_structure_factor ?
                 _parse_species(_cfg_get(sf_cfg, "species", "both"), "observables.structure_factor") :
                 String[]

    same_site_convention = "factorial"
    if compute_density_density || compute_structure_factor
        same_site_convention = lowercase(String(_cfg_get(dd_cfg, "same_site_convention", "factorial")))
        same_site_convention in ("factorial", "plain") ||
            error("observables.density_density.same_site_convention must be \"factorial\" or \"plain\"")
    end
    sf_factorial_diagonal = (same_site_convention == "factorial")
    max_r = _cfg_get(dd_cfg, "max_r", nothing)
    max_r = max_r === nothing ? nothing : Int(max_r)
    fold_min_image = _parse_bool(_cfg_get(dd_cfg, "fold_min_image", false), false)

    needed_ops = union(dd_species, sf_species)
    nvec_a = nothing
    nn_a = nothing
    nvec_b = nothing
    nn_b = nothing
    if !isempty(needed_ops)
        if "Na" in needed_ops
            nvec_a, nn_a = density_density_matrix(
                psi,
                sites,
                "Na";
                same_site_convention=same_site_convention
            )
        end
        if "Nb" in needed_ops
            nvec_b, nn_b = density_density_matrix(
                psi,
                sites,
                "Nb";
                same_site_convention=same_site_convention
            )
        end
    end

    cnn_a = nothing
    r_a = nothing
    g_a = nothing
    c_a = nothing
    anchors_a = nothing
    cnn_b = nothing
    r_b = nothing
    g_b = nothing
    c_b = nothing
    anchors_b = nothing
    if compute_density_density
        if "Na" in dd_species
            progress && println("Evaluating observable: density_density (species=a)")
            cnn_a = connected_density_density_matrix(nvec_a, nn_a)
            r_a, g_a, c_a, anchors_a = transl_avg_density_density(
                nvec_a,
                nn_a;
                periodic=periodic_eff,
                max_r=max_r,
                fold_min_image=fold_min_image
            )
        end
        if "Nb" in dd_species
            progress && println("Evaluating observable: density_density (species=b)")
            cnn_b = connected_density_density_matrix(nvec_b, nn_b)
            r_b, g_b, c_b, anchors_b = transl_avg_density_density(
                nvec_b,
                nn_b;
                periodic=periodic_eff,
                max_r=max_r,
                fold_min_image=fold_min_image
            )
        end
    elseif progress
        println("Skipping observable: density_density (not requested)")
    end

    k_a = nothing
    sf_a = nothing
    sfc_a = nothing
    k_b = nothing
    sf_b = nothing
    sfc_b = nothing
    if compute_structure_factor
        if "Na" in sf_species
            progress && println("Evaluating observable: structure_factor (species=a)")
            k_a, sf_a = structure_factor_from_nn(
                nvec_a,
                nn_a;
                connected=false,
                factorial_diagonal=sf_factorial_diagonal
            )
            _, sfc_a = structure_factor_from_nn(
                nvec_a,
                nn_a;
                connected=true,
                factorial_diagonal=sf_factorial_diagonal
            )
        end
        if "Nb" in sf_species
            progress && println("Evaluating observable: structure_factor (species=b)")
            k_b, sf_b = structure_factor_from_nn(
                nvec_b,
                nn_b;
                connected=false,
                factorial_diagonal=sf_factorial_diagonal
            )
            _, sfc_b = structure_factor_from_nn(
                nvec_b,
                nn_b;
                connected=true,
                factorial_diagonal=sf_factorial_diagonal
            )
        end
    elseif progress
        println("Skipping observable: structure_factor (not requested)")
    end

    tc_pairs = Tuple{Int,Int}[]
    tc_anchors = nothing
    tc_a = nothing
    tc_b = nothing
    if compute_triple_corr
        tc_species = _parse_species(_cfg_get(tc_cfg, "species", "both"), "observables.triple_corr")
        tc_pairs = _resolve_triple_pairs(tc_cfg, length(sites))
        tc_precompute = _parse_bool(_cfg_get(tc_cfg, "precompute", true), true)

        tc_anchors = zeros(Int, length(tc_pairs))
        tc_a = ("Na" in tc_species) ? zeros(Float64, length(tc_pairs)) : nothing
        tc_b = ("Nb" in tc_species) ? zeros(Float64, length(tc_pairs)) : nothing

        npre_a = nothing
        nnpre_a = nothing
        npre_b = nothing
        nnpre_b = nothing
        if tc_precompute
            if tc_a !== nothing
                npre_a = precompute_n(psi, sites, "Na")
                nnpre_a = precompute_nn(psi, sites, "Na")
            end
            if tc_b !== nothing
                npre_b = precompute_n(psi, sites, "Nb")
                nnpre_b = precompute_nn(psi, sites, "Nb")
            end
        end

        for (idx, (r, s)) in enumerate(tc_pairs)
            if tc_a !== nothing
                if tc_precompute && npre_a !== nothing && nnpre_a !== nothing
                    C, N = transl_avg_connected_nnn_no_cached(
                        psi,
                        sites,
                        "Na",
                        r,
                        s;
                        periodic=periodic_eff,
                        nvec=npre_a,
                        nnmat=nnpre_a
                    )
                else
                    C, N = transl_avg_connected_nnn_no(psi, sites, "Na", r, s; periodic=periodic_eff)
                end
                tc_a[idx] = C
                tc_anchors[idx] = N
            end
            if tc_b !== nothing
                if tc_precompute && npre_b !== nothing && nnpre_b !== nothing
                    C, N = transl_avg_connected_nnn_no_cached(
                        psi,
                        sites,
                        "Nb",
                        r,
                        s;
                        periodic=periodic_eff,
                        nvec=npre_b,
                        nnmat=nnpre_b
                    )
                else
                    C, N = transl_avg_connected_nnn_no(psi, sites, "Nb", r, s; periodic=periodic_eff)
                end
                tc_b[idx] = C
                tc_anchors[idx] = N
            end
        end
    end

    return (
        energy=energy,
        densities=(na=na_v, nb=nb_v),
        totals=(Na=Na, Nb=Nb),
        density_density=(
            requested=compute_density_density,
            same_site_convention=same_site_convention,
            nn_a=nn_a,
            connected_nn_a=cnn_a,
            r_a=r_a,
            transl_avg_nn_a=g_a,
            transl_avg_connected_nn_a=c_a,
            anchors_a=anchors_a,
            nn_b=nn_b,
            connected_nn_b=cnn_b,
            r_b=r_b,
            transl_avg_nn_b=g_b,
            transl_avg_connected_nn_b=c_b,
            anchors_b=anchors_b
        ),
        structure_factor=(
            requested=compute_structure_factor,
            k_a=k_a,
            S_a=sf_a,
            S_connected_a=sfc_a,
            k_b=k_b,
            S_b=sf_b,
            S_connected_b=sfc_b
        ),
        triple_corr=(
            requested=compute_triple_corr,
            pairs=tc_pairs,
            pairs_matrix=_pairs_matrix(tc_pairs),
            anchors=tc_anchors,
            C_no_a=tc_a,
            C_no_b=tc_b
        )
    )
end

"""
    _replace_group(parent, name)

Delete group `name` under `parent` if it exists, then recreate it.
Used to prevent stale observables datasets across rewrites.
"""
@inline function _replace_group(parent, name::AbstractString)
    if haskey(parent, name)
        HDF5.delete_object(parent, name)
    end
    return HDF5.create_group(parent, name)
end

"""
    write_observables_hdf5!(f, obs; schema_version=OBSERVABLES_SCHEMA_VERSION)

Write canonical observables content to an open HDF5 file handle `f`.
The `/observables` group is replaced atomically to avoid stale datasets.
"""
function write_observables_hdf5!(f, obs; schema_version::AbstractString=OBSERVABLES_SCHEMA_VERSION)
    g_meta = ensure_group(f, "meta")
    write_or_replace(g_meta, "observables_schema_version", String(schema_version))

    g_obs = _replace_group(f, "observables")

    g_energy = HDF5.create_group(g_obs, "energy")
    if obs.energy !== nothing
        write_or_replace(g_energy, "E0", obs.energy)
    end

    g_den = HDF5.create_group(g_obs, "densities")
    write_or_replace(g_den, "na", obs.densities.na)
    write_or_replace(g_den, "nb", obs.densities.nb)

    g_tot = HDF5.create_group(g_obs, "totals")
    write_or_replace(g_tot, "Na", obs.totals.Na)
    write_or_replace(g_tot, "Nb", obs.totals.Nb)

    dd = obs.density_density
    if dd.requested
        g_dd = HDF5.create_group(g_obs, "density_density")
        write_or_replace(g_dd, "same_site_convention", dd.same_site_convention)
        if dd.nn_a !== nothing
            write_or_replace(g_dd, "nn_a", dd.nn_a)
            write_or_replace(g_dd, "connected_nn_a", dd.connected_nn_a)
            write_or_replace(g_dd, "r_a", dd.r_a)
            write_or_replace(g_dd, "transl_avg_nn_a", dd.transl_avg_nn_a)
            write_or_replace(g_dd, "transl_avg_connected_nn_a", dd.transl_avg_connected_nn_a)
            write_or_replace(g_dd, "anchors_a", dd.anchors_a)
        end
        if dd.nn_b !== nothing
            write_or_replace(g_dd, "nn_b", dd.nn_b)
            write_or_replace(g_dd, "connected_nn_b", dd.connected_nn_b)
            write_or_replace(g_dd, "r_b", dd.r_b)
            write_or_replace(g_dd, "transl_avg_nn_b", dd.transl_avg_nn_b)
            write_or_replace(g_dd, "transl_avg_connected_nn_b", dd.transl_avg_connected_nn_b)
            write_or_replace(g_dd, "anchors_b", dd.anchors_b)
        end
    end

    sf = obs.structure_factor
    if sf.requested
        g_sf = HDF5.create_group(g_obs, "structure_factor")
        if sf.S_a !== nothing
            write_or_replace(g_sf, "k_a", sf.k_a)
            write_or_replace(g_sf, "S_a", sf.S_a)
            write_or_replace(g_sf, "S_connected_a", sf.S_connected_a)
        end
        if sf.S_b !== nothing
            write_or_replace(g_sf, "k_b", sf.k_b)
            write_or_replace(g_sf, "S_b", sf.S_b)
            write_or_replace(g_sf, "S_connected_b", sf.S_connected_b)
        end
    end

    tc = obs.triple_corr
    if tc.requested
        g_tc = HDF5.create_group(g_obs, "triple_corr")
        write_or_replace(g_tc, "pairs", tc.pairs_matrix)
        write_or_replace(g_tc, "anchors", tc.anchors)
        if tc.C_no_a !== nothing
            write_or_replace(g_tc, "C_no_a", tc.C_no_a)
        end
        if tc.C_no_b !== nothing
            write_or_replace(g_tc, "C_no_b", tc.C_no_b)
        end
    end

    return nothing
end
