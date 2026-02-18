# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

# ----------------------------
# 4) Run DMRG
# ----------------------------

using LinearAlgebra

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
    sweep_offset::Int
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
    checkpoint_density_every::Int=1,
    sweep_offset::Int=0
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
        max(1, checkpoint_density_every),
        max(0, sweep_offset)
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
            nb=nb,
            checkpoint_sweep=sweep
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

    sw_local = sweep === nothing ? length(obs.energies) : Int(sweep)
    sw = sw_local + obs.sweep_offset
    psi = haskey(kwargs, :psi) ? kwargs[:psi] : nothing
    if psi !== nothing
        maybe_checkpoint!(obs; psi=psi, energy=energy, sweep=sw, outputlevel=outputlevel)
    end

    energy_active = obs.energy_tol > 0.0
    trunc_active = obs.trunc_tol > 0.0
    if !(energy_active || trunc_active)
        return false
    end
    if sw_local < max(2, obs.min_sweeps)
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
                "Early stopping DMRG at local sweep $sw_local (global sweep $sw): " *
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

function expand_maxdim_schedule(maxdim, nsweeps::Int)
    nsweeps >= 0 || error("nsweeps must be non-negative, got $nsweeps")
    if nsweeps == 0
        return Int[]
    end

    function _dict_get(d::AbstractDict, keys::Vector{String}, default=nothing)
        for k in keys
            if haskey(d, k)
                return d[k]
            end
            ks = Symbol(k)
            if haskey(d, ks)
                return d[ks]
            end
        end
        return default
    end

    function _auto_warmup_schedule(spec::AbstractDict)
        maxv_raw = _dict_get(spec, ["max", "maximum"], nothing)
        maxv_raw === nothing && error("maxdim warmup spec requires key \"max\" (or \"maximum\")")
        maxv = Int(maxv_raw)
        maxv >= 1 || error("maxdim max must be >= 1, got $maxv")

        minv = Int(_dict_get(spec, ["min", "start"], min(50, maxv)))
        minv >= 1 || error("maxdim min/start must be >= 1, got $minv")
        minv <= maxv || error("maxdim min/start must be <= max, got min=$minv max=$maxv")

        warmup_raw = _dict_get(spec, ["warmup_sweeps", "steps", "length"], nothing)
        if warmup_raw !== nothing
            nwarm = Int(warmup_raw)
            nwarm >= 1 || error("maxdim warmup_sweeps/steps/length must be >= 1, got $nwarm")
            nwarm = min(nwarm, nsweeps)

            if nwarm == 1
                return [maxv]
            end

            # Geometric ramp from min -> max, then `expand_maxdim_schedule` pads with max.
            ratio = maxv == minv ? 1.0 : (maxv / minv)^(1 / (nwarm - 1))
            vals = Vector{Int}(undef, nwarm)
            vals[1] = minv
            for i in 2:nwarm
                vals[i] = max(vals[i - 1], Int(round(minv * ratio^(i - 1))))
            end
            vals[end] = maxv
            return vals
        end

        # Default warmup: doubling ramp to max, then pad with max.
        vals = Int[minv]
        while vals[end] < maxv
            push!(vals, min(maxv, 2 * vals[end]))
        end
        return vals
    end

    sched = if maxdim isa AbstractVector
        vals = Int.(collect(maxdim))
        isempty(vals) && error("maxdim vector cannot be empty")
        vals
    elseif maxdim isa AbstractDict
        spec = maxdim
        mode_raw = _dict_get(spec, ["mode"], "warmup")
        mode = lowercase(String(mode_raw))
        mode in ("warmup", "auto", "automatic") ||
            error("maxdim.mode must be one of: warmup, auto, automatic (got: $mode_raw)")
        _auto_warmup_schedule(spec)
    else
        [Int(maxdim)]
    end

    if length(sched) >= nsweeps
        return sched[1:nsweeps]
    end
    return vcat(sched, fill(sched[end], nsweeps - length(sched)))
end

function current_params_sha256(params_path)::Union{Nothing,String}
    if params_path === nothing || !isfile(params_path)
        return nothing
    end
    return bytes2hex(SHA.sha256(read(params_path, String)))
end

function parse_bool(x)
    x isa Bool && return x
    s = lowercase(strip(String(x)))
    if s in ("1", "true", "yes", "y", "on")
        return true
    elseif s in ("0", "false", "no", "n", "off")
        return false
    end
    error("Cannot parse boolean value from: $x")
end

function sites_compatible_for_resume(sites_requested, sites_checkpoint)
    if length(sites_requested) != length(sites_checkpoint)
        return false
    end
    for i in eachindex(sites_requested)
        si = sites_requested[i]
        sj = sites_checkpoint[i]
        if dim(si) != dim(sj)
            return false
        end
        if hasqns(si) != hasqns(sj)
            return false
        end
    end
    return true
end

function load_checkpoint_for_resume(;
    checkpoint_path,
    checkpoint_params_path,
    checkpoint_require_hash::Bool,
    requested_sites,
    outputlevel::Integer=0
)
    if checkpoint_path === nothing || !isfile(checkpoint_path)
        return nothing
    end
    try
        st = load_state(checkpoint_path)
        if checkpoint_require_hash
            current_hash = current_params_sha256(checkpoint_params_path)
            if current_hash === nothing
                outputlevel > 0 && println(
                    "Found checkpoint at $checkpoint_path but current params hash is unavailable; not resuming."
                )
                return nothing
            end
            if st.params_sha256 === nothing || st.params_sha256 != current_hash
                outputlevel > 0 && println(
                    "Checkpoint hash mismatch at $checkpoint_path; not resuming."
                )
                return nothing
            end
        end
        if !sites_compatible_for_resume(requested_sites, st.sites)
            outputlevel > 0 && println(
                "Checkpoint site structure incompatible with requested setup; not resuming."
            )
            return nothing
        end
        return st
    catch err
        outputlevel > 0 && println("Failed to load checkpoint at $checkpoint_path: $err")
        return nothing
    end
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
    resume_from_checkpoint=true,
    resume_mode::Union{Symbol,AbstractString}=:remaining,
    checkpoint_require_hash=true,
    checkpoint_save_densities=false,
    checkpoint_density_every=1,
    outputlevel=1,
    saveresults=false, savepath="results.h5", kwargs...)

    if outputlevel > 0
        @info "Thread config" julia_threads=Threads.nthreads() gc_threads=Threads.ngcthreads() blas_threads=BLAS.get_num_threads() julia_num_threads_env=get(ENV, "JULIA_NUM_THREADS", "unset") openblas_num_threads_env=get(ENV, "OPENBLAS_NUM_THREADS", "unset")
    end

    resume_mode_sym = resume_mode isa Symbol ? resume_mode : Symbol(lowercase(String(resume_mode)))
    resume_mode_sym in (:remaining, :warm_start) ||
        error("resume_mode must be :remaining or :warm_start (got: $resume_mode)")

    requested_sites = two_boson_siteinds(L; nmax_a=nmax_a, nmax_b=nmax_b, conserve_qns=conserve_qns)
    resume_from_checkpoint_flag = parse_bool(resume_from_checkpoint)
    checkpoint_require_hash_flag = parse_bool(checkpoint_require_hash)

    st_checkpoint = resume_from_checkpoint_flag ? load_checkpoint_for_resume(
        checkpoint_path=checkpoint_path,
        checkpoint_params_path=checkpoint_params_path,
        checkpoint_require_hash=checkpoint_require_hash_flag,
        requested_sites=requested_sites,
        outputlevel=outputlevel
    ) : nothing

    checkpoint_sweep = 0
    if st_checkpoint !== nothing
        sites = st_checkpoint.sites
        psi0 = st_checkpoint.psi
        checkpoint_sweep = st_checkpoint.checkpoint_sweep === nothing ? 0 : Int(st_checkpoint.checkpoint_sweep)
        if outputlevel > 0
            println("Loaded checkpoint from $checkpoint_path (checkpoint_sweep=$checkpoint_sweep)")
        end
    else
        sites = requested_sites
        # Initial state in the correct (Na, Nb) sector
        conf, _, _ = dmrg_initial_configuration(;
            L=L,
            nmax_a=nmax_a, nmax_b=nmax_b,
            Na_total=Na_total, Nb_total=Nb_total,
            impurity_distribution=impurity_distribution,
            seed=seed
        )
        psi0 = product_state_mps(sites, conf; nmax_b=nmax_b)
    end

    H = build_hamiltonian(
        sites;
        t_a=t_a, t_b=t_b,
        U_a=U_a, U_b=U_b, U_ab=U_ab,
        mu_a=mu_a, mu_b=mu_b,
        periodic=periodic
    )
    maxdim_schedule = expand_maxdim_schedule(maxdim, Int(nsweeps))
    sweep_offset = 0
    if st_checkpoint !== nothing && resume_mode_sym == :remaining && checkpoint_sweep > 0
        if checkpoint_sweep >= nsweeps
            if outputlevel > 0
                println(
                    "Checkpoint already reached requested nsweeps (checkpoint_sweep=$checkpoint_sweep, nsweeps=$nsweeps). Skipping DMRG."
                )
            end
            energy = st_checkpoint.energy
            if energy === nothing
                energy = real(inner(psi0, Apply(H, psi0)))
            end
            return energy, psi0, sites, H
        end
        maxdim_schedule = maxdim_schedule[(checkpoint_sweep + 1):end]
        sweep_offset = checkpoint_sweep
        outputlevel > 0 && println("Resuming remaining sweeps: $(length(maxdim_schedule)) (offset=$sweep_offset)")
    elseif st_checkpoint !== nothing && resume_mode_sym == :warm_start
        outputlevel > 0 && println("Warm-starting from checkpoint with full sweep schedule.")
    end
    isempty(maxdim_schedule) && error("No sweeps to run. Set nsweeps > 0.")

    sweeps = Sweeps(length(maxdim_schedule))
    maxdim!(sweeps, maxdim_schedule...)
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
        checkpoint_density_every=Int(checkpoint_density_every),
        sweep_offset=sweep_offset
    )

    energy, psi = dmrg(H, psi0, sweeps; outputlevel=outputlevel, observer=observer)
    return energy, psi, sites, H
end
