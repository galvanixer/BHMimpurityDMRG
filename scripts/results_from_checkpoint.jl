# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

include(joinpath(@__DIR__, "..", "src", "BHMimpurityDMRG.jl"))
using .BHMimpurityDMRG
using HDF5
using YAML
using SHA
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
    observables_path_arg = length(ARGS) >= 3 ? ARGS[3] : nothing

    isfile(checkpoint_path) || error("Checkpoint not found: $checkpoint_path")
    st = load_state(checkpoint_path)

    cfg_base = st.params_yaml === nothing ? Dict{String,Any}() : normalize_yaml(YAML.load(st.params_yaml))
    merged_cfg = with_observables_config(cfg_base; observables_path=observables_path_arg)
    cfg = merged_cfg.cfg
    observables_path = merged_cfg.observables_path
    observables_loaded = merged_cfg.observables_loaded
    params_has_observables = haskey(cfg_base, "observables") || haskey(cfg_base, :observables)

    obs_cfg = get(cfg, "observables", Dict{String,Any}())
    lattice_cfg = get(cfg, "lattice", Dict{String,Any}())
    periodic = parse_bool(get(lattice_cfg, "periodic", true), true)

    psi = st.psi
    sites = st.sites
    energy = compute_energy_from_state(st, cfg)

    if observables_loaded
        println("Using observables config from: $(abspath(observables_path))")
    elseif params_has_observables
        println("Observables file not found ($observables_path); using observables from checkpoint parameters.")
    else
        println("Observables file not found ($observables_path); using default observables settings.")
    end

    na, nb = if st.na !== nothing && st.nb !== nothing
        st.na, st.nb
    else
        measure_densities(psi, sites)
    end

    dd_requested = haskey(obs_cfg, "density_density")
    sf_requested = haskey(obs_cfg, "structure_factor")
    if haskey(obs_cfg, "triple_corr")
        println("Note: results_from_checkpoint.jl ignores observables.triple_corr by design.")
    end

    obs = compute_observables(
        psi,
        sites;
        energy=energy,
        na=na,
        nb=nb,
        cfg=cfg,
        periodic=periodic,
        compute_density_density=dd_requested,
        compute_structure_factor=sf_requested,
        compute_triple_corr=false,
        progress=true
    )

    HDF5.h5open(results_path, "w") do f
        g_meta = ensure_group(f, "meta")
        write_meta!(
            g_meta;
            params_text=st.params_yaml,
            state_path=checkpoint_path,
            state_params_sha256=st.params_sha256
        )
        write_results_schema!(g_meta)
        if observables_loaded
            write_or_replace(g_meta, "observables_path", abspath(observables_path))
            write_or_replace(g_meta, "observables_sha256", bytes2hex(SHA.sha256(read(observables_path, String))))
        end
        if st.checkpoint_sweep !== nothing
            write_or_replace(g_meta, "checkpoint_sweep", Int(st.checkpoint_sweep))
        end

        write_observables_hdf5!(f, obs)
    end

    println("Wrote checkpoint-derived observables to: $(abspath(results_path))")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
