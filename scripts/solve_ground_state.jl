# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

include(joinpath(@__DIR__, "..", "src", "BHMimpurityDMRG.jl"))
using .BHMimpurityDMRG
using HDF5
using Logging
using SHA

function main()
    params_path = length(ARGS) >= 1 ? ARGS[1] :
                  get(ENV, "PARAMS", joinpath(@__DIR__, "..", "configs", "parameters.yaml"))

    # Load parameters with fallback to empty dict if file not found or loading fails
    cfg = isfile(params_path) ? load_params(params_path) : Dict{String,Any}()
    dmrg_cfg = merge_sections(cfg, ["lattice", "local_hilbert", "initial_state", "hamiltonian", "dmrg"])
    io_cfg = get(cfg, "io", Dict{String,Any}())
    results_path = get(io_cfg, "results_path", "results.h5")
    state_path = get(io_cfg, "state_save_path", "dmrg_state.h5")
    save_state_flag = get(io_cfg, "save_state", true)
    logcfg = setup_logger(io_cfg; default_log_path="run.log")
    logger = logcfg.logger

    with_logger(logger) do
        @info "Solving for ground state using DMRG" params_path=params_path results_path=results_path
        params_text = isfile(params_path) ? read(params_path, String) : nothing
        current_hash = params_text === nothing ? nothing : bytes2hex(SHA.sha256(params_text))
        _, init_na, init_nb = dmrg_initial_configuration(; dmrg_cfg...)

        na = nothing
        nb = nothing
        if isfile(state_path)
            st = load_state(state_path)
            if current_hash !== nothing && st.params_sha256 !== nothing &&
               st.params_sha256 == current_hash && st.energy !== nothing
                psi = st.psi
                sites = st.sites
                energy = st.energy
                na = st.na
                nb = st.nb
                @info "Using cached state (hash matched)" state_path=state_path
            else
                @info "Cached state missing/mismatched hash; running DMRG" state_path=state_path
                energy, psi, sites, _ = run_dmrg(; checkpoint_params_path=params_path, dmrg_cfg...)
                if save_state_flag
                    na_tmp, nb_tmp = measure_densities(psi, sites)
                    save_state(
                        state_path,
                        psi;
                        energy=energy,
                        params_path=params_path,
                        na=na_tmp,
                        nb=nb_tmp,
                        init_na=init_na,
                        init_nb=init_nb
                    )
                    na = na_tmp
                    nb = nb_tmp
                    @info "Saved state" state_path=state_path
                end
            end
        else
            @info "No cached state found; running DMRG" state_path=state_path
            energy, psi, sites, _ = run_dmrg(; checkpoint_params_path=params_path, dmrg_cfg...)
            if save_state_flag
                na_tmp, nb_tmp = measure_densities(psi, sites)
                save_state(
                    state_path,
                    psi;
                    energy=energy,
                    params_path=params_path,
                    na=na_tmp,
                    nb=nb_tmp,
                    init_na=init_na,
                    init_nb=init_nb
                )
                na = na_tmp
                nb = nb_tmp
                @info "Saved state" state_path=state_path
            end
        end

        if na === nothing || nb === nothing
            na, nb = measure_densities(psi, sites)
        end
        Na, Nb = total_numbers(na, nb)
        obs_cfg = get(cfg, "observables", Dict{String,Any}())
        dd_cfg = get(obs_cfg, "density_density", Dict{String,Any}())
        sf_cfg = get(obs_cfg, "structure_factor", Dict{String,Any}())
        periodic = get(get(cfg, "lattice", Dict{String,Any}()), "periodic", true)
        species = lowercase(String(get(dd_cfg, "species", "both")))
        sf_species = lowercase(String(get(sf_cfg, "species", "both")))
        max_r = get(dd_cfg, "max_r", nothing)
        max_r = max_r === nothing ? nothing : Int(max_r)
        fold_min_image = get(dd_cfg, "fold_min_image", false)
        same_site_convention = lowercase(String(get(dd_cfg, "same_site_convention", "factorial")))
        same_site_convention in ("factorial", "plain") ||
            error("observables.density_density.same_site_convention must be \"factorial\" or \"plain\"")
        sf_factorial_diagonal = (same_site_convention == "factorial")
        dd_opnames = species == "a" ? ["Na"] :
                     species == "b" ? ["Nb"] :
                     species == "both" ? ["Na", "Nb"] :
                     error("observables.density_density.species must be one of: \"a\", \"b\", \"both\"")
        sf_opnames = sf_species == "a" ? ["Na"] :
                     sf_species == "b" ? ["Nb"] :
                     sf_species == "both" ? ["Na", "Nb"] :
                     error("observables.structure_factor.species must be one of: \"a\", \"b\", \"both\"")
        needed_opnames = union(dd_opnames, sf_opnames)

        nvec_a = nothing    # Density vector for species a
        nn_a = nothing      # Density-density matrix for species a
        cnn_a = nothing     # Connected density-density matrix for species a
        r_a = nothing       # Distances for translational averaging for species a
        g_a = nothing       # Translationally averaged density-density for species a
        c_a = nothing       # Translationally averaged connected density-density for species a
        anchors_a = nothing # Anchor site indices for translational averaging for species a

        nvec_b = nothing    # Density vector for species b
        nn_b = nothing      # Density-density matrix for species b
        cnn_b = nothing     # Connected density-density matrix for species b
        r_b = nothing       # Distances for translational averaging for species b
        g_b = nothing       # Translationally averaged density-density for species b
        c_b = nothing       # Translationally averaged connected density-density for species b
        anchors_b = nothing # Anchor site indices for translational averaging for species b

        if "Na" in needed_opnames
            nvec_a, nn_a = density_density_matrix(
                psi,
                sites,
                "Na";
                same_site_convention=same_site_convention
            )
        end
        if "Nb" in needed_opnames
            nvec_b, nn_b = density_density_matrix(
                psi,
                sites,
                "Nb";
                same_site_convention=same_site_convention
            )
        end

        if "Na" in dd_opnames
            cnn_a = connected_density_density_matrix(nvec_a, nn_a)
            r_a, g_a, c_a, anchors_a = transl_avg_density_density(
                nvec_a,
                nn_a;
                periodic=periodic,
                max_r=max_r,
                fold_min_image=fold_min_image
            )
        end
        if "Nb" in dd_opnames
            cnn_b = connected_density_density_matrix(nvec_b, nn_b)
            r_b, g_b, c_b, anchors_b = transl_avg_density_density(
                nvec_b,
                nn_b;
                periodic=periodic,
                max_r=max_r,
                fold_min_image=fold_min_image
            )
        end

        k_a = nothing
        sf_a = nothing
        sfc_a = nothing
        if "Na" in sf_opnames
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
        k_b = nothing
        sf_b = nothing
        sfc_b = nothing
        if "Nb" in sf_opnames
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

        mode = isfile(results_path) ? "r+" : "w"
        HDF5.h5open(results_path, mode) do f
            g_meta = ensure_group(f, "meta")
            write_meta!(g_meta;
                params_path=params_path,
                state_path=state_path,
                state_params_sha256=current_hash
            )
            g_obs = ensure_group(f, "observables")
            g_energy = ensure_group(g_obs, "energy")
            g_den = ensure_group(g_obs, "densities")
            g_tot = ensure_group(g_obs, "totals")
            write_or_replace(g_energy, "E0", energy)
            write_or_replace(g_den, "na", na)
            write_or_replace(g_den, "nb", nb)
            write_or_replace(g_tot, "Na", Na)
            write_or_replace(g_tot, "Nb", Nb)

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
        @info "Wrote results" results_path=results_path
    end

    close(logcfg.logio)

    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
