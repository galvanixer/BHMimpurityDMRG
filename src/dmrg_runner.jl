# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

# ----------------------------
# 4) Run DMRG
# ----------------------------

mutable struct EarlyStopDMRGObserver <: ITensorMPS.AbstractObserver
    energies::Vector{Float64}
    truncerrs::Vector{Float64}
    energy_tol::Float64
    trunc_tol::Float64
    min_sweeps::Int
    patience::Int
    streak::Int
    sites::Any
    checkpoint_every::Int
    checkpoint_path::Union{Nothing,String}
    checkpoint_params_path::Union{Nothing,String}
    checkpoint_save_densities::Bool
    checkpoint_density_every::Int
end

function EarlyStopDMRGObserver(;
    energy_tol::Real=0.0,
    trunc_tol::Real=0.0,
    min_sweeps::Int=2,
    patience::Int=1,
    sites=nothing,
    checkpoint_every::Int=0,
    checkpoint_path::Union{Nothing,AbstractString}=nothing,
    checkpoint_params_path::Union{Nothing,AbstractString}=nothing,
    checkpoint_save_densities::Bool=false,
    checkpoint_density_every::Int=1
)
    return EarlyStopDMRGObserver(
        Float64[],
        Float64[],
        Float64(energy_tol),
        Float64(trunc_tol),
        min_sweeps,
        max(1, patience),
        0,
        sites,
        checkpoint_every,
        checkpoint_path === nothing ? nothing : String(checkpoint_path),
        checkpoint_params_path === nothing ? nothing : String(checkpoint_params_path),
        checkpoint_save_densities,
        max(1, checkpoint_density_every)
    )
end

function checkpoint_densities(psi::MPS, sites)
    L = length(sites)
    na = zeros(Float64, L)
    nb = zeros(Float64, L)
    for i in 1:L
        na[i] = expect_n(psi, sites, "Na", i)
        nb[i] = expect_n(psi, sites, "Nb", i)
    end
    return na, nb
end

function maybe_checkpoint!(
    obs::EarlyStopDMRGObserver;
    psi,
    energy,
    sweep::Int,
    outputlevel::Integer=0
)
    if obs.checkpoint_every <= 0 || obs.checkpoint_path === nothing
        return nothing
    end
    if sweep % obs.checkpoint_every != 0
        return nothing
    end

    path = obs.checkpoint_path
    tmp_path = path * ".tmp"
    checkpoint_index = sweep ÷ obs.checkpoint_every
    save_densities = obs.checkpoint_save_densities &&
        (checkpoint_index % obs.checkpoint_density_every == 0)
    try
        na = nothing
        nb = nothing
        if save_densities
            sites = obs.sites === nothing ? siteinds(psi) : obs.sites
            na, nb = checkpoint_densities(psi, sites)
        end
        save_state(
            tmp_path,
            psi;
            energy=energy,
            sites=obs.sites === nothing ? siteinds(psi) : obs.sites,
            params_path=obs.checkpoint_params_path,
            na=na,
            nb=nb
        )
        mv(tmp_path, path; force=true)
        if outputlevel > 0
            println("Wrote DMRG checkpoint at sweep $sweep to $path (densities_saved=$save_densities)")
        end
    catch err
        if outputlevel > 0
            println("Warning: failed to write DMRG checkpoint at sweep $sweep: $err")
        end
        if isfile(tmp_path)
            rm(tmp_path; force=true)
        end
    end
    return nothing
end

function ITensorMPS.measure!(obs::EarlyStopDMRGObserver; kwargs...)
    half_sweep = kwargs[:half_sweep]
    b = kwargs[:bond]
    psi = kwargs[:psi]
    truncerr = kwargs[:spec].truncerr

    if half_sweep == 2
        N = length(psi)
        if b == (N - 1)
            push!(obs.truncerrs, 0.0)
        end
        if isempty(obs.truncerrs)
            push!(obs.truncerrs, truncerr)
        else
            obs.truncerrs[end] = max(obs.truncerrs[end], truncerr)
        end
    end
    return nothing
end

function ITensorMPS.checkdone!(obs::EarlyStopDMRGObserver; outputlevel=0, energy=nothing, sweep=nothing, kwargs...)
    energy === nothing && return false
    push!(obs.energies, Float64(real(energy)))

    sw = sweep === nothing ? length(obs.energies) : Int(sweep)
    psi = haskey(kwargs, :psi) ? kwargs[:psi] : nothing
    if psi !== nothing
        maybe_checkpoint!(obs; psi=psi, energy=energy, sweep=sw, outputlevel=outputlevel)
    end

    energy_active = obs.energy_tol > 0.0
    trunc_active = obs.trunc_tol > 0.0
    if !(energy_active || trunc_active)
        return false
    end
    if sw < obs.min_sweeps || length(obs.energies) < 2
        obs.streak = 0
        return false
    end

    dE = abs(obs.energies[end] - obs.energies[end - 1])
    energy_ok = !energy_active || (dE < obs.energy_tol)

    maxerr = isempty(obs.truncerrs) ? Inf : obs.truncerrs[end]
    trunc_ok = !trunc_active || (maxerr < obs.trunc_tol)

    if energy_ok && trunc_ok
        obs.streak += 1
    else
        obs.streak = 0
    end

    if obs.streak >= obs.patience
        if outputlevel > 0
            println(
                "Early stopping DMRG at sweep $sw: " *
                "dE=$(dE), maxerr=$(maxerr), streak=$(obs.streak)"
            )
        end
        return true
    end
    return false
end

function dmrg_initial_configuration(; L=12,
    nmax_a=3, nmax_b=1,
    Na_total=12, Nb_total=1,
    impurity_distribution::Union{Symbol,AbstractString}=:center,
    seed::Union{Int,Nothing}=nothing,
    kwargs...)
    impdist = impurity_distribution isa Symbol ? impurity_distribution : Symbol(impurity_distribution)
    conf = initial_configuration(L;
        Na_total=Na_total, Nb_total=Nb_total,
        impurity_distribution=impdist,
        nmax_a=nmax_a,
        nmax_b=nmax_b,
        seed=seed
    )
    na0 = Float64[first(x) for x in conf]
    nb0 = Float64[last(x) for x in conf]
    return conf, na0, nb0
end

function run_dmrg(; L=12,
    nmax_a=3, nmax_b=1,
    conserve_qns=true,
    Na_total=12, Nb_total=1,
    impurity_distribution::Union{Symbol,AbstractString}=:center,
    seed::Union{Int,Nothing}=nothing,
    t_a=1.0, t_b=1.0,
    U_a=10.0, U_b=0.0, U_ab=5.0,
    mu_a=0.0, mu_b=0.0,
    nsweeps=12, periodic=true,
    cutoff=1e-10,
    maxdim=[50, 100, 200, 400, 600, 800, 800, 800, 800, 800, 800, 800],
    energy_tol=0.0,
    trunc_tol=0.0,
    min_sweeps=2,
    patience=1,
    checkpoint_every=0,
    checkpoint_path=nothing,
    checkpoint_params_path=nothing,
    checkpoint_save_densities=false,
    checkpoint_density_every=1,
    outputlevel=1,
    saveresults=false, savepath="results.h5", kwargs...)

    sites = two_boson_siteinds(L; nmax_a=nmax_a, nmax_b=nmax_b, conserve_qns=conserve_qns)

    # Initial state in the correct (Na, Nb) sector
    conf, _, _ = dmrg_initial_configuration(;
        L=L,
        nmax_a=nmax_a, nmax_b=nmax_b,
        Na_total=Na_total, Nb_total=Nb_total,
        impurity_distribution=impurity_distribution,
        seed=seed
    )
    psi0 = product_state_mps(sites, conf; nmax_b=nmax_b)

    H = build_hamiltonian(
        sites;
        t_a=t_a, t_b=t_b,
        U_a=U_a, U_b=U_b, U_ab=U_ab,
        mu_a=mu_a, mu_b=mu_b,
        periodic=periodic
    )
    sweeps = Sweeps(nsweeps)
    # A robust default schedule (adjust as needed)
    if maxdim isa AbstractVector
        maxdim!(sweeps, (Int(x) for x in maxdim)...)
    else
        maxdim!(sweeps, Int(maxdim))
    end
    cutoff!(sweeps, cutoff)
    # Uncomment noise if you see convergence to excited states/local minima:
    # noise!(sweeps, 1e-6, 1e-7, 1e-8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    observer = EarlyStopDMRGObserver(;
        energy_tol=energy_tol,
        trunc_tol=trunc_tol,
        min_sweeps=min_sweeps,
        patience=patience,
        sites=sites,
        checkpoint_every=Int(checkpoint_every),
        checkpoint_path=checkpoint_path,
        checkpoint_params_path=checkpoint_params_path,
        checkpoint_save_densities=checkpoint_save_densities,
        checkpoint_density_every=Int(checkpoint_density_every)
    )

    energy, psi = dmrg(H, psi0, sweeps; outputlevel=outputlevel, observer=observer)
    return energy, psi, sites, H
end
