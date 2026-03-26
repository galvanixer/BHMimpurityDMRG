# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

# ----------------------------
# 4) Run DMRG
# ----------------------------

using LinearAlgebra

struct DMRGSweepTraceRow
    energy::Float64
    delta_energy::Float64
    max_truncerr::Float64
    maxdim_used::Int64
    walltime_sec::Float64
end

mutable struct EarlyStopDMRGObserver <: ITensorMPS.AbstractObserver
    energies::Vector{Float64}
    truncerrs::Vector{Float64}
    maxdim_used::Vector{Int64}
    walltime_sec::Vector{Float64}
    checkpoint_sweeps::Vector{Int64}
    energy_tol::Float64
    trunc_tol::Float64
    min_sweeps::Int
    patience::Int
    streak::Int
    sites::Any
    hamiltonian::Any
    checkpoint_every::Int
    checkpoint_path::Union{Nothing,String}
    checkpoint_params_path::Union{Nothing,String}
    checkpoint_save_densities::Bool
    checkpoint_density_every::Int
    sweep_offset::Int
    resume_mode::Symbol
    checkpoint_sweep_start::Int
    init_na::Union{Nothing,Vector{Float64}}
    init_nb::Union{Nothing,Vector{Float64}}
    early_stop_triggered::Bool
    run_start_time_ns::Int64
    last_sweep_time_ns::Int64
end

function EarlyStopDMRGObserver(;
    energy_tol::Real=0.0,
    trunc_tol::Real=0.0,
    min_sweeps::Int=2,
    patience::Int=1,
    sites=nothing,
    hamiltonian=nothing,
    checkpoint_every::Int=0,
    checkpoint_path::Union{Nothing,AbstractString}=nothing,
    checkpoint_params_path::Union{Nothing,AbstractString}=nothing,
    checkpoint_save_densities::Bool=false,
    checkpoint_density_every::Int=1,
    sweep_offset::Int=0,
    resume_mode::Symbol=:unknown,
    checkpoint_sweep_start::Int=0,
    init_na=nothing,
    init_nb=nothing
)
    init_na_v = init_na === nothing ? nothing : Float64.(collect(init_na))
    init_nb_v = init_nb === nothing ? nothing : Float64.(collect(init_nb))
    t0 = Int64(time_ns())
    return EarlyStopDMRGObserver(
        Float64[],
        Float64[],
        Int64[],
        Float64[],
        Int64[],
        Float64(energy_tol),
        Float64(trunc_tol),
        min_sweeps,
        max(1, patience),
        0,
        sites,
        hamiltonian,
        checkpoint_every,
        checkpoint_path === nothing ? nothing : String(checkpoint_path),
        checkpoint_params_path === nothing ? nothing : String(checkpoint_params_path),
        checkpoint_save_densities,
        max(1, checkpoint_density_every),
        max(0, sweep_offset),
        resume_mode,
        max(0, checkpoint_sweep_start),
        init_na_v,
        init_nb_v,
        false,
        t0,
        t0
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
        return false
    end
    if sweep % obs.checkpoint_every != 0
        return false
    end

    path = obs.checkpoint_path
    tmp_path = path * ".tmp"
    checkpoint_index = sweep ÷ obs.checkpoint_every
    save_densities = obs.checkpoint_save_densities &&
        (checkpoint_index % obs.checkpoint_density_every == 0)
    try
        na = nothing
        nb = nothing
        energy_variance = obs.hamiltonian === nothing ? nothing :
                          operator_variance(psi, obs.hamiltonian; expectation=energy)
        if save_densities
            sites = obs.sites === nothing ? siteinds(psi) : obs.sites
            na, nb = checkpoint_densities(psi, sites)
        end
        save_state(
            tmp_path,
            psi;
            energy=energy,
            energy_variance=energy_variance,
            sites=obs.sites === nothing ? siteinds(psi) : obs.sites,
            params_path=obs.checkpoint_params_path,
            na=na,
            nb=nb,
            init_na=obs.init_na,
            init_nb=obs.init_nb,
            checkpoint_sweep=sweep
        )
        mv(tmp_path, path; force=true)
        if outputlevel > 0
            println("Wrote DMRG checkpoint at sweep $sweep to $path (densities_saved=$save_densities)")
        end
        return true
    catch err
        if outputlevel > 0
            println("Warning: failed to write DMRG checkpoint at sweep $sweep: $err")
        end
        if isfile(tmp_path)
            rm(tmp_path; force=true)
        end
        return false
    end
end

function persist_checkpoint_diagnostics!(
    obs::EarlyStopDMRGObserver;
    outputlevel::Integer=0
)
    path = obs.checkpoint_path
    path === nothing && return false
    try
        diag = build_dmrg_diagnostics(
            obs;
            resume_mode=obs.resume_mode,
            checkpoint_sweep_start=obs.checkpoint_sweep_start
        )
        HDF5.h5open(path, "r+") do f
            write_dmrg_diagnostics!(f, diag)
        end
        return true
    catch err
        if outputlevel > 0
            println("Warning: failed to append DMRG diagnostics to checkpoint $path: $err")
        end
        return false
    end
end

@inline function max_bond_dim(psi::MPS)
    N = length(psi)
    N <= 1 && return 1
    md = 1
    for b in 1:(N - 1)
        li = linkind(psi, b)
        if li !== nothing
            md = max(md, dim(li))
        end
    end
    return md
end

function build_dmrg_diagnostics(
    obs::EarlyStopDMRGObserver;
    resume_mode::Symbol,
    checkpoint_sweep_start::Int
)
    n = length(obs.energies)
    delta_energy = fill(NaN, n)
    for i in 2:n
        delta_energy[i] = abs(obs.energies[i] - obs.energies[i - 1])
    end

    trunc = length(obs.truncerrs) == n ? obs.truncerrs : begin
        out = fill(NaN, n)
        m = min(length(obs.truncerrs), n)
        for i in 1:m
            out[i] = obs.truncerrs[i]
        end
        out
    end
    maxdim_used = length(obs.maxdim_used) == n ? obs.maxdim_used : begin
        out = fill(Int64(0), n)
        m = min(length(obs.maxdim_used), n)
        for i in 1:m
            out[i] = obs.maxdim_used[i]
        end
        out
    end
    walltime = length(obs.walltime_sec) == n ? obs.walltime_sec : begin
        out = fill(NaN, n)
        m = min(length(obs.walltime_sec), n)
        for i in 1:m
            out[i] = obs.walltime_sec[i]
        end
        out
    end

    rows = Vector{DMRGSweepTraceRow}(undef, n)
    for i in 1:n
        rows[i] = DMRGSweepTraceRow(
            obs.energies[i],
            delta_energy[i],
            trunc[i],
            Int64(maxdim_used[i]),
            walltime[i]
        )
    end

    converged = false
    if n >= max(2, obs.min_sweeps)
        energy_ok = obs.energy_tol <= 0.0 || delta_energy[end] < obs.energy_tol
        trunc_ok = obs.trunc_tol <= 0.0 || (!isempty(trunc) && isfinite(trunc[end]) && trunc[end] < obs.trunc_tol)
        converged = energy_ok && trunc_ok
    end

    return (
        sweep_trace=rows,
        checkpoint_sweeps=Int64.(obs.checkpoint_sweeps),
        sweeps_completed=n,
        converged=converged,
        early_stop_triggered=obs.early_stop_triggered,
        energy_tol=obs.energy_tol,
        trunc_tol=obs.trunc_tol,
        patience=obs.patience,
        min_sweeps=obs.min_sweeps,
        resume_mode=String(resume_mode),
        checkpoint_sweep_start=Int(checkpoint_sweep_start)
    )
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
    checkpoint_written = false
    if psi !== nothing
        checkpoint_written = maybe_checkpoint!(obs; psi=psi, energy=energy, sweep=sw, outputlevel=outputlevel)
        push!(obs.maxdim_used, Int64(max_bond_dim(psi)))
    else
        push!(obs.maxdim_used, Int64(0))
    end
    now_ns = Int64(time_ns())
    push!(obs.walltime_sec, (now_ns - obs.last_sweep_time_ns) / 1e9)
    obs.last_sweep_time_ns = now_ns
    if checkpoint_written
        push!(obs.checkpoint_sweeps, Int64(sw))
        persist_checkpoint_diagnostics!(obs; outputlevel=outputlevel)
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
        obs.early_stop_triggered = true
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
    impurity_distribution::Union{Symbol,AbstractString}=:centered_pileup,
    seed::Union{Int,Nothing}=nothing,
    impurity_sites=nothing,
    kwargs...)
    impdist = impurity_distribution isa Symbol ? impurity_distribution : Symbol(impurity_distribution)
    conf = initial_configuration(L;
        Na_total=Na_total, Nb_total=Nb_total,
        impurity_distribution=impdist,
        nmax_a=nmax_a,
        nmax_b=nmax_b,
        seed=seed,
        impurity_sites=impurity_sites
    )
    na0 = Float64[first(x) for x in conf]
    nb0 = Float64[last(x) for x in conf]
    return conf, na0, nb0
end

function _spec_get(d::AbstractDict, key::String, default=nothing)
    if haskey(d, key)
        return d[key]
    end
    ks = Symbol(key)
    if haskey(d, ks)
        return d[ks]
    end
    return default
end

function _spec_keys(d::AbstractDict)
    keys_out = String[]
    for k in keys(d)
        push!(keys_out, k isa Symbol ? String(k) : string(k))
    end
    return keys_out
end

function _require_allowed_spec_keys(d::AbstractDict, allowed_keys::Vector{String}, what::AbstractString)
    allowed = Set(allowed_keys)
    unexpected = sort!(filter(k -> !(k in allowed), _spec_keys(d)))
    allowed_list = join(sort!(copy(allowed_keys)), ", ")
    unexpected_list = join(unexpected, ", ")
    isempty(unexpected) || error(
        "$what only accepts keys $allowed_list; unexpected key(s): $unexpected_list"
    )
    return nothing
end

function expand_maxdim_schedule(maxdim, nsweeps::Int)
    nsweeps >= 0 || error("nsweeps must be non-negative, got $nsweeps")
    if nsweeps == 0
        return Int[]
    end

    function _warmup_schedule(spec::AbstractDict)
        _require_allowed_spec_keys(spec, ["mode", "min", "max", "sweeps"], "maxdim dict")

        mode_raw = _spec_get(spec, "mode", nothing)
        mode_raw === nothing && error("maxdim dict requires mode=\"warmup\"")
        mode = lowercase(String(mode_raw))
        mode == "warmup" || error("maxdim.mode must be \"warmup\" (got: $mode_raw)")

        maxv_raw = _spec_get(spec, "max", nothing)
        maxv_raw === nothing && error("maxdim dict requires key \"max\"")
        maxv = Int(maxv_raw)
        maxv >= 1 || error("maxdim max must be >= 1, got $maxv")

        minv = Int(_spec_get(spec, "min", min(50, maxv)))
        minv >= 1 || error("maxdim min must be >= 1, got $minv")
        minv <= maxv || error("maxdim min must be <= max, got min=$minv max=$maxv")

        sweeps_raw = _spec_get(spec, "sweeps", nothing)
        if sweeps_raw !== nothing
            nsched = Int(sweeps_raw)
            nsched >= 1 || error("maxdim sweeps must be >= 1, got $nsched")
            nsched = min(nsched, nsweeps)

            if nsched == 1
                return [maxv]
            end

            # Geometric ramp from min -> max, then `expand_maxdim_schedule` pads with max.
            ratio = maxv == minv ? 1.0 : (maxv / minv)^(1 / (nsched - 1))
            vals = Vector{Int}(undef, nsched)
            vals[1] = minv
            for i in 2:nsched
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
        _warmup_schedule(maxdim)
    else
        [Int(maxdim)]
    end

    all(v -> v >= 1, sched) || error("all maxdim values must be >= 1")
    if length(sched) >= nsweeps
        return sched[1:nsweeps]
    end
    return vcat(sched, fill(sched[end], nsweeps - length(sched)))
end

function expand_cutoff_schedule(cutoff, nsweeps::Int)
    nsweeps >= 0 || error("nsweeps must be non-negative, got $nsweeps")
    if nsweeps == 0
        return Float64[]
    end

    function _geometric_schedule(spec::AbstractDict)
        _require_allowed_spec_keys(spec, ["mode", "start", "stop", "sweeps"], "cutoff dict")

        mode_raw = _spec_get(spec, "mode", nothing)
        mode_raw === nothing && error("cutoff dict requires mode=\"geometric\"")
        mode = lowercase(String(mode_raw))
        mode == "geometric" || error("cutoff.mode must be \"geometric\" (got: $mode_raw)")

        start_raw = _spec_get(spec, "start", nothing)
        stop_raw = _spec_get(spec, "stop", nothing)
        start_raw === nothing && error("cutoff dict requires key \"start\"")
        stop_raw === nothing && error("cutoff dict requires key \"stop\"")

        start_v = Float64(start_raw)
        stop_v = Float64(stop_raw)
        start_v > 0 || error("cutoff start must be > 0, got $start_v")
        stop_v > 0 || error("cutoff stop must be > 0, got $stop_v")

        sweeps_raw = _spec_get(spec, "sweeps", nothing)
        nsched = sweeps_raw === nothing ? nsweeps : Int(sweeps_raw)
        nsched >= 1 || error("cutoff sweeps must be >= 1, got $nsched")
        nsched = min(nsched, nsweeps)

        if nsched == 1
            return [stop_v]
        end

        # Geometric ramp for cutoff so multi-decade schedules are natural.
        ratio = (stop_v / start_v)^(1 / (nsched - 1))
        vals = Vector{Float64}(undef, nsched)
        vals[1] = start_v
        for i in 2:nsched
            vals[i] = vals[i - 1] * ratio
        end
        vals[end] = stop_v
        return vals
    end

    sched = if cutoff isa AbstractVector
        vals = Float64.(collect(cutoff))
        isempty(vals) && error("cutoff vector cannot be empty")
        vals
    elseif cutoff isa AbstractDict
        _geometric_schedule(cutoff)
    else
        [Float64(cutoff)]
    end

    all(v -> v > 0, sched) || error("all cutoff values must be > 0")
    if length(sched) >= nsweeps
        return sched[1:nsweeps]
    end
    return vcat(sched, fill(sched[end], nsweeps - length(sched)))
end

function expand_noise_schedule(noise, nsweeps::Int)
    nsweeps >= 0 || error("nsweeps must be non-negative, got $nsweeps")
    if nsweeps == 0
        return Float64[]
    end
    if noise === nothing
        return fill(0.0, nsweeps)
    end

    function _linear_schedule(start_v::Float64, stop_v::Float64, nwarm::Int)
        if nwarm == 1
            return [stop_v]
        end
        vals = Vector{Float64}(undef, nwarm)
        for i in 1:nwarm
            t = (i - 1) / (nwarm - 1)
            vals[i] = (1 - t) * start_v + t * stop_v
        end
        vals[end] = stop_v
        return vals
    end

    function _geometric_schedule(start_v::Float64, stop_v::Float64, nwarm::Int)
        start_v > 0 || error("noise geometric start must be > 0, got $start_v")
        stop_v > 0 || error("noise geometric stop must be > 0, got $stop_v")
        if nwarm == 1
            return [stop_v]
        end
        ratio = (stop_v / start_v)^(1 / (nwarm - 1))
        vals = Vector{Float64}(undef, nwarm)
        vals[1] = start_v
        for i in 2:nwarm
            vals[i] = vals[i - 1] * ratio
        end
        vals[end] = stop_v
        return vals
    end

    function _bursts_schedule(spec::AbstractDict)
        _require_allowed_spec_keys(spec, ["mode", "bursts"], "noise dict")

        bursts_raw = _spec_get(spec, "bursts", nothing)
        bursts_raw isa AbstractVector || error("noise bursts mode requires key \"bursts\" as a non-empty array")
        isempty(bursts_raw) && error("noise bursts array cannot be empty")

        sched = fill(0.0, nsweeps)
        occupied = falses(nsweeps)

        for (i, burst_raw) in enumerate(bursts_raw)
            burst_raw isa AbstractDict || error("noise burst #$i must be a dictionary")
            burst = burst_raw
            _require_allowed_spec_keys(
                burst,
                ["start_sweep", "start", "stop", "sweeps"],
                "noise burst #$i"
            )

            start_sweep_raw = _spec_get(burst, "start_sweep", nothing)
            start_sweep_raw === nothing && error("noise burst #$i requires key \"start_sweep\"")
            start_sweep = Int(start_sweep_raw)
            start_sweep >= 1 || error("noise burst #$i start_sweep must be >= 1, got $start_sweep")
            start_sweep <= nsweeps || error(
                "noise burst #$i start_sweep must be <= nsweeps=$nsweeps, got $start_sweep"
            )

            start_raw = _spec_get(burst, "start", nothing)
            stop_raw = _spec_get(burst, "stop", nothing)
            start_raw === nothing && error("noise burst #$i requires key \"start\"")
            stop_raw === nothing && error("noise burst #$i requires key \"stop\"")

            start_v = Float64(start_raw)
            stop_v = Float64(stop_raw)
            start_v >= 0 || error("noise burst #$i start must be >= 0, got $start_v")
            stop_v >= 0 || error("noise burst #$i stop must be >= 0, got $stop_v")

            sweeps_raw = _spec_get(burst, "sweeps", nothing)
            sweeps_raw === nothing && error("noise burst #$i requires key \"sweeps\"")
            burst_sweeps = Int(sweeps_raw)
            burst_sweeps >= 1 || error("noise burst #$i sweeps must be >= 1, got $burst_sweeps")

            stop_sweep = min(nsweeps, start_sweep + burst_sweeps - 1)
            any(@view occupied[start_sweep:stop_sweep]) && error(
                "noise bursts must not overlap; burst #$i overlaps a previous burst"
            )

            local_sched = _linear_schedule(start_v, stop_v, stop_sweep - start_sweep + 1)
            sched[start_sweep:stop_sweep] .= local_sched
            occupied[start_sweep:stop_sweep] .= true
        end

        return sched
    end

    function _generated_noise_schedule(spec::AbstractDict)
        mode_raw = _spec_get(spec, "mode", nothing)
        mode_raw === nothing && error("noise dict requires mode=\"linear\", mode=\"geometric\", or mode=\"bursts\"")
        mode = lowercase(String(mode_raw))
        mode in ("linear", "geometric", "bursts") ||
            error("noise.mode must be \"linear\", \"geometric\", or \"bursts\" (got: $mode_raw)")

        if mode == "bursts"
            return _bursts_schedule(spec)
        end

        _require_allowed_spec_keys(spec, ["mode", "start", "stop", "sweeps"], "noise dict")

        start_raw = _spec_get(spec, "start", nothing)
        stop_raw = _spec_get(spec, "stop", nothing)
        start_raw === nothing && error("noise dict requires key \"start\"")
        stop_raw === nothing && error("noise dict requires key \"stop\"")

        start_v = Float64(start_raw)
        stop_v = Float64(stop_raw)
        start_v >= 0 || error("noise start must be >= 0, got $start_v")
        stop_v >= 0 || error("noise stop must be >= 0, got $stop_v")

        sweeps_raw = _spec_get(spec, "sweeps", nothing)
        nsched = sweeps_raw === nothing ? nsweeps : Int(sweeps_raw)
        nsched >= 1 || error("noise sweeps must be >= 1, got $nsched")
        nsched = min(nsched, nsweeps)

        if mode == "linear"
            return _linear_schedule(start_v, stop_v, nsched)
        end
        return _geometric_schedule(start_v, stop_v, nsched)
    end

    sched = if noise isa AbstractVector
        vals = Float64.(collect(noise))
        isempty(vals) && error("noise vector cannot be empty")
        vals
    elseif noise isa AbstractDict
        _generated_noise_schedule(noise)
    else
        [Float64(noise)]
    end

    all(v -> v >= 0, sched) || error("all noise values must be >= 0")
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

"""
    sites_compatible_for_resume(sites_requested, sites_checkpoint)

Check if two sets of sites are compatible for resuming a DMRG calculation.

This function verifies that the requested sites and checkpoint sites have matching
properties, ensuring that a previous DMRG calculation can be safely resumed with
the new site configuration.

# Arguments
- `sites_requested`: The sites configuration for the current/new DMRG run
- `sites_checkpoint`: The sites configuration from a previous checkpoint

# Returns
- `Bool`: `true` if the sites are compatible for resume, `false` otherwise

# Compatibility criteria
Two site configurations are considered compatible if:
1. They have the same number of sites
2. Each corresponding site has the same physical dimension
3. Each corresponding site has the same quantum number structure (both with or without QNs)

"""
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
    impurity_distribution::Union{Symbol,AbstractString}=:centered_pileup,
    seed::Union{Int,Nothing}=nothing,
    impurity_sites=nothing,
    t_a=1.0, t_b=1.0,
    U_a=10.0, U_b=0.0, U_ab=5.0,
    mu_a=0.0, mu_b=0.0,
    nsweeps=12, periodic=true,
    cutoff=1e-10,
    noise=0.0,
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
    return_diagnostics=false,
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
    init_na = nothing
    init_nb = nothing
    if st_checkpoint !== nothing
        sites = st_checkpoint.sites
        psi0 = st_checkpoint.psi
        checkpoint_sweep = st_checkpoint.checkpoint_sweep === nothing ? 0 : Int(st_checkpoint.checkpoint_sweep)
        init_na = st_checkpoint.init_na
        init_nb = st_checkpoint.init_nb
        if init_na === nothing || init_nb === nothing
            _, init_na, init_nb = dmrg_initial_configuration(;
                L=L,
                nmax_a=nmax_a, nmax_b=nmax_b,
                Na_total=Na_total, Nb_total=Nb_total,
                impurity_distribution=impurity_distribution,
                seed=seed,
                impurity_sites=impurity_sites
            )
        end
        if outputlevel > 0
            println("Loaded checkpoint from $checkpoint_path (checkpoint_sweep=$checkpoint_sweep)")
        end
    else
        sites = requested_sites
        # Initial state in the correct (Na, Nb) sector
        conf, init_na, init_nb = dmrg_initial_configuration(;
            L=L,
            nmax_a=nmax_a, nmax_b=nmax_b,
            Na_total=Na_total, Nb_total=Nb_total,
            impurity_distribution=impurity_distribution,
            seed=seed,
            impurity_sites=impurity_sites
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
    cutoff_schedule = expand_cutoff_schedule(cutoff, Int(nsweeps))
    noise_schedule = expand_noise_schedule(noise, Int(nsweeps))
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
            if parse_bool(return_diagnostics)
                diag = (
                    sweep_trace=DMRGSweepTraceRow[],
                    checkpoint_sweeps=Int64[],
                    sweeps_completed=0,
                    converged=false,
                    early_stop_triggered=false,
                    energy_tol=Float64(energy_tol),
                    trunc_tol=Float64(trunc_tol),
                    patience=Int(patience),
                    min_sweeps=Int(min_sweeps),
                    resume_mode=String(resume_mode_sym),
                    checkpoint_sweep_start=Int(checkpoint_sweep)
                )
                return energy, psi0, sites, H, diag
            end
            return energy, psi0, sites, H
        end
        maxdim_schedule = maxdim_schedule[(checkpoint_sweep + 1):end]
        cutoff_schedule = cutoff_schedule[(checkpoint_sweep + 1):end]
        noise_schedule = noise_schedule[(checkpoint_sweep + 1):end]
        sweep_offset = checkpoint_sweep
        outputlevel > 0 && println("Resuming remaining sweeps: $(length(maxdim_schedule)) (offset=$sweep_offset)")
    elseif st_checkpoint !== nothing && resume_mode_sym == :warm_start
        outputlevel > 0 && println("Warm-starting from checkpoint with full sweep schedule.")
    end
    isempty(maxdim_schedule) && error("No sweeps to run. Set nsweeps > 0.")
    length(maxdim_schedule) == length(cutoff_schedule) ||
        error(
            "internal schedule mismatch: maxdim has $(length(maxdim_schedule)) sweeps, cutoff has $(length(cutoff_schedule)) sweeps"
        )
    length(maxdim_schedule) == length(noise_schedule) ||
        error(
            "internal schedule mismatch: maxdim has $(length(maxdim_schedule)) sweeps, noise has $(length(noise_schedule)) sweeps"
        )

    sweeps = Sweeps(length(maxdim_schedule))
    maxdim!(sweeps, maxdim_schedule...)
    cutoff!(sweeps, cutoff_schedule...)
    noise!(sweeps, noise_schedule...)

    observer = EarlyStopDMRGObserver(;
        energy_tol=energy_tol,
        trunc_tol=trunc_tol,
        min_sweeps=min_sweeps,
        patience=patience,
        sites=sites,
        hamiltonian=H,
        checkpoint_every=Int(checkpoint_every),
        checkpoint_path=checkpoint_path,
        checkpoint_params_path=checkpoint_params_path,
        checkpoint_save_densities=checkpoint_save_densities,
        checkpoint_density_every=Int(checkpoint_density_every),
        sweep_offset=sweep_offset,
        resume_mode=resume_mode_sym,
        checkpoint_sweep_start=checkpoint_sweep,
        init_na=init_na,
        init_nb=init_nb
    )

    energy, psi = dmrg(H, psi0, sweeps; outputlevel=outputlevel, observer=observer)
    if parse_bool(return_diagnostics)
        diag = build_dmrg_diagnostics(
            observer;
            resume_mode=resume_mode_sym,
            checkpoint_sweep_start=checkpoint_sweep
        )
        return energy, psi, sites, H, diag
    end
    return energy, psi, sites, H
end
