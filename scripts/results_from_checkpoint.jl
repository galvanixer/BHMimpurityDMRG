# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

include(joinpath(@__DIR__, "..", "src", "BHMimpurityDMRG.jl"))
using .BHMimpurityDMRG
using HDF5
using YAML
using ITensors
using ITensorMPS

function normalize_yaml(x)
    if x isa AbstractDict
        return Dict{String,Any}(String(k) => normalize_yaml(v) for (k, v) in x)
    elseif x isa AbstractVector
        return [normalize_yaml(v) for v in x]
    else
        return x
    end
end

function parse_bool(x, default::Bool=false)
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

function parse_species(species, section_name::String)
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

function compute_energy_from_state(st, cfg::AbstractDict)
    if st.energy !== nothing
        return Float64(real(st.energy))
    end
    ham_cfg = get(cfg, "hamiltonian", Dict{String,Any}())
    lattice_cfg = get(cfg, "lattice", Dict{String,Any}())
    H = build_hamiltonian(
        st.sites;
        t_a=Float64(get(ham_cfg, "t_a", 1.0)),
        t_b=Float64(get(ham_cfg, "t_b", 1.0)),
        U_a=Float64(get(ham_cfg, "U_a", 10.0)),
        U_b=Float64(get(ham_cfg, "U_b", 0.0)),
        U_ab=Float64(get(ham_cfg, "U_ab", 5.0)),
        mu_a=Float64(get(ham_cfg, "mu_a", 0.0)),
        mu_b=Float64(get(ham_cfg, "mu_b", 0.0)),
        periodic=parse_bool(get(lattice_cfg, "periodic", true), true)
    )
    return Float64(real(inner(st.psi, Apply(H, st.psi))))
end

function main()
    checkpoint_path = length(ARGS) >= 1 ? ARGS[1] : "dmrg_state_checkpoint.h5"
    results_path = length(ARGS) >= 2 ? ARGS[2] : "results_checkpoint.h5"

    isfile(checkpoint_path) || error("Checkpoint not found: $checkpoint_path")
    st = load_state(checkpoint_path)

    cfg = st.params_yaml === nothing ? Dict{String,Any}() : normalize_yaml(YAML.load(st.params_yaml))
    obs_cfg = get(cfg, "observables", Dict{String,Any}())
    lattice_cfg = get(cfg, "lattice", Dict{String,Any}())
    periodic = parse_bool(get(lattice_cfg, "periodic", true), true)

    psi = st.psi
    sites = st.sites
    energy = compute_energy_from_state(st, cfg)

    println("Evaluating observable: energy")

    println("Evaluating observable: densities / totals")
    na, nb = if st.na !== nothing && st.nb !== nothing
        st.na, st.nb
    else
        measure_densities(psi, sites)
    end
    Na, Nb = total_numbers(na, nb)

    dd_requested = haskey(obs_cfg, "density_density")
    sf_requested = haskey(obs_cfg, "structure_factor")
    if haskey(obs_cfg, "triple_corr")
        println("Note: results_from_checkpoint.jl ignores observables.triple_corr by design.")
    end

    dd_cfg = get(obs_cfg, "density_density", Dict{String,Any}())
    sf_cfg = get(obs_cfg, "structure_factor", Dict{String,Any}())

    dd_species = dd_requested ?
        parse_species(get(dd_cfg, "species", "both"), "observables.density_density") :
        String[]
    sf_species = sf_requested ?
        parse_species(get(sf_cfg, "species", "both"), "observables.structure_factor") :
        String[]

    needed_dd_ops = union(dd_species, sf_species)
    same_site_convention = lowercase(String(get(dd_cfg, "same_site_convention", "factorial")))
    same_site_convention in ("factorial", "plain") ||
        error("observables.density_density.same_site_convention must be \"factorial\" or \"plain\"")
    sf_factorial_diagonal = (same_site_convention == "factorial")
    max_r = get(dd_cfg, "max_r", nothing)
    max_r = max_r === nothing ? nothing : Int(max_r)
    fold_min_image = parse_bool(get(dd_cfg, "fold_min_image", false), false)

    nvec_a = nothing
    nn_a = nothing
    nvec_b = nothing
    nn_b = nothing
    if dd_requested || sf_requested
        if "Na" in needed_dd_ops
            nvec_a, nn_a = density_density_matrix(
                psi,
                sites,
                "Na";
                same_site_convention=same_site_convention
            )
        end
        if "Nb" in needed_dd_ops
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

    if dd_requested && "Na" in dd_species
        println("Evaluating observable: density_density (species=a)")
        cnn_a = connected_density_density_matrix(nvec_a, nn_a)
        r_a, g_a, c_a, anchors_a = transl_avg_density_density(
            nvec_a,
            nn_a;
            periodic=periodic,
            max_r=max_r,
            fold_min_image=fold_min_image
        )
    end
    if dd_requested && "Nb" in dd_species
        println("Evaluating observable: density_density (species=b)")
        cnn_b = connected_density_density_matrix(nvec_b, nn_b)
        r_b, g_b, c_b, anchors_b = transl_avg_density_density(
            nvec_b,
            nn_b;
            periodic=periodic,
            max_r=max_r,
            fold_min_image=fold_min_image
        )
    end
    if !dd_requested
        println("Skipping observable: density_density (not requested)")
    end

    k_a = nothing
    sf_a = nothing
    sfc_a = nothing
    k_b = nothing
    sf_b = nothing
    sfc_b = nothing
    if sf_requested && "Na" in sf_species
        println("Evaluating observable: structure_factor (species=a)")
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
    if sf_requested && "Nb" in sf_species
        println("Evaluating observable: structure_factor (species=b)")
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
    if !sf_requested
        println("Skipping observable: structure_factor (not requested)")
    end

    HDF5.h5open(results_path, "w") do f
        g_meta = ensure_group(f, "meta")
        write_meta!(
            g_meta;
            params_text=st.params_yaml,
            state_path=checkpoint_path,
            state_params_sha256=st.params_sha256
        )
        if st.checkpoint_sweep !== nothing
            write_or_replace(g_meta, "checkpoint_sweep", Int(st.checkpoint_sweep))
        end

        g_obs = ensure_group(f, "observables")
        g_energy = ensure_group(g_obs, "energy")
        g_den = ensure_group(g_obs, "densities")
        g_tot = ensure_group(g_obs, "totals")

        write_or_replace(g_energy, "E0", energy)
        write_or_replace(g_den, "na", na)
        write_or_replace(g_den, "nb", nb)
        write_or_replace(g_tot, "Na", Na)
        write_or_replace(g_tot, "Nb", Nb)

        if dd_requested
            g_dd = ensure_group(g_obs, "density_density")
            write_or_replace(g_dd, "same_site_convention", same_site_convention)
            if nn_a !== nothing
                write_or_replace(g_dd, "nn_a", nn_a)
                write_or_replace(g_dd, "connected_nn_a", cnn_a)
                write_or_replace(g_dd, "r_a", r_a)
                write_or_replace(g_dd, "transl_avg_nn_a", g_a)
                write_or_replace(g_dd, "transl_avg_connected_nn_a", c_a)
                write_or_replace(g_dd, "anchors_a", anchors_a)
            end
            if nn_b !== nothing
                write_or_replace(g_dd, "nn_b", nn_b)
                write_or_replace(g_dd, "connected_nn_b", cnn_b)
                write_or_replace(g_dd, "r_b", r_b)
                write_or_replace(g_dd, "transl_avg_nn_b", g_b)
                write_or_replace(g_dd, "transl_avg_connected_nn_b", c_b)
                write_or_replace(g_dd, "anchors_b", anchors_b)
            end
        end

        if sf_requested
            g_sf = ensure_group(g_obs, "structure_factor")
            if sf_a !== nothing
                write_or_replace(g_sf, "k_a", k_a)
                write_or_replace(g_sf, "S_a", sf_a)
                write_or_replace(g_sf, "S_connected_a", sfc_a)
            end
            if sf_b !== nothing
                write_or_replace(g_sf, "k_b", k_b)
                write_or_replace(g_sf, "S_b", sf_b)
                write_or_replace(g_sf, "S_connected_b", sfc_b)
            end
        end

    end

    println("Wrote checkpoint-derived observables to: $(abspath(results_path))")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
