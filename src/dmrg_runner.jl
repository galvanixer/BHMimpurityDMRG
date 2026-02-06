# Author: Tanul Gupta <tanulgupta123@gmail.com>
# License: MIT
# Copyright (c) 2026 Tanul Gupta

# ----------------------------
# 4) Run DMRG
# ----------------------------

function run_dmrg(; L=12,
    nmax_a=3, nmax_b=1,
    conserve_qns=true,
    Na_total=12, Nb_total=1,
    t_a=1.0, t_b=1.0,
    U_a=10.0, U_b=0.0, U_ab=5.0,
    mu_a=0.0, mu_b=0.0,
    nsweeps=12)

    sites = two_boson_siteinds(L; nmax_a=nmax_a, nmax_b=nmax_b, conserve_qns=conserve_qns)

    # Initial state in the correct (Na, Nb) sector
    conf = initial_configuration(L; Na_total=Na_total, Nb_total=Nb_total)
    psi0 = product_state_mps(sites, conf; nmax_b=nmax_b)

    H = build_hamiltonian(
        sites;
        t_a=t_a, t_b=t_b,
        U_a=U_a, U_b=U_b, U_ab=U_ab,
        mu_a=mu_a, mu_b=mu_b,
        periodic=true
    )
    sweeps = Sweeps(nsweeps)
    # A robust default schedule (adjust as needed)
    maxdim!(sweeps, 50, 100, 200, 400, 600, 800, 800, 800, 800, 800, 800, 800)
    cutoff!(sweeps, 1e-10)
    # Uncomment noise if you see convergence to excited states/local minima:
    # noise!(sweeps, 1e-6, 1e-7, 1e-8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    energy, psi = dmrg(H, psi0, sweeps; outputlevel=1)
    return energy, psi, sites, H
end
