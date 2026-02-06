# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

include(joinpath(@__DIR__, "..", "src", "BHMimpurityDMRG.jl"))
using .BHMimpurityDMRG

function main()
    params_path = get(ENV, "PARAMS", joinpath(@__DIR__, "..", "parameters.yaml"))
    cfg = isfile(params_path) ? load_params(params_path) : Dict{String,Any}()

    dmrg_cfg = haskey(cfg, "dmrg") ? dict_to_namedtuple(cfg["dmrg"]) : (;)
    energy, psi, sites, H = run_dmrg(; dmrg_cfg...)

    println("Ground-state energy = ", energy)

    # Measure translationally averaged connected 3-point correlator:
    # C^(3)(r, s) = (1/L) Σ_i ⟨n_i n_{i+r} n_{i+s}⟩_c
    #
    # Since `run_dmrg` currently builds a periodic Hamiltonian, use periodic=true here.
    periodic = get(get(cfg, "triple_corr", Dict{String,Any}()), "periodic", true)

    # Example separations (r, s). Add more as needed.
    raw_pairs = get(get(cfg, "triple_corr", Dict{String,Any}()), "pairs",
        [[0, 0], [0, 1], [0, 2], [1, 1], [1, 2], [1, 3], [2, 3]])
    pairs = [(Int(p[1]), Int(p[2])) for p in raw_pairs]

    for (r, s) in pairs
        # Normal-ordered definition:
        # if indices coincide (e.g. r=0 or s=0), uses n(n-1) and n(n-1)(n-2).
        Caa, N = transl_avg_connected_nnn_no(psi, sites, "Na", r, s; periodic=periodic)
        Cbb, _ = transl_avg_connected_nnn_no(psi, sites, "Nb", r, s; periodic=periodic)
        println("C_no^(3)_a(r=$r, s=$s) = $Caa   (anchors=$N)")
        println("C_no^(3)_b(r=$r, s=$s) = $Cbb   (anchors=$N)")
    end

    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
