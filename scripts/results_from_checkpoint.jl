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

function compute_energy_stats_from_state(st, cfg::AbstractDict)
    need_hamiltonian = st.energy === nothing || st.energy_variance === nothing
    H = need_hamiltonian ? build_hamiltonian_from_config(st.sites, cfg) : nothing
    energy = st.energy !== nothing ? Float64(real(st.energy)) : expect_operator(st.psi, H)
    energy_variance = st.energy_variance !== nothing ?
                      Float64(real(st.energy_variance)) :
                      operator_variance(st.psi, H; expectation=energy)
    return energy, energy_variance
end

function compute_initial_densities_from_state(st, cfg::AbstractDict)
    if st.init_na !== nothing && st.init_nb !== nothing
        return st.init_na, st.init_nb
    end
    init_cfg = merge_sections(cfg, ["lattice", "local_hilbert", "initial_state"])
    _, init_na, init_nb = dmrg_initial_configuration(; init_cfg...)
    return init_na, init_nb
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
    energy, energy_variance = compute_energy_stats_from_state(st, cfg)
    init_na, init_nb = compute_initial_densities_from_state(st, cfg)

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
    pdd_requested = haskey(obs_cfg, "pair_distance_distribution")
    spdm_requested = haskey(obs_cfg, "single_particle_density_matrix")
    if haskey(obs_cfg, "triple_corr")
        println("Note: results_from_checkpoint.jl ignores observables.triple_corr by design.")
    end

    obs = compute_observables(
        psi,
        sites;
        energy=energy,
        energy_variance=energy_variance,
        init_na=init_na,
        init_nb=init_nb,
        na=na,
        nb=nb,
        cfg=cfg,
        periodic=periodic,
        compute_density_density=dd_requested,
        compute_structure_factor=sf_requested,
        compute_pair_distance_distribution=pdd_requested,
        compute_single_particle_density_matrix=spdm_requested,
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
        if st.dmrg_diagnostics !== nothing
            write_dmrg_diagnostics!(f, st.dmrg_diagnostics)
        end
    end

    println("Wrote checkpoint-derived observables to: $(abspath(results_path))")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
