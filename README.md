# BHMimpurityDMRG

DMRG simulations for a two-species Bose-Hubbard model with impurities, custom boson site types, observables, and binding-energy analysis. Built on ITensors/ITensorMPS.

## Features
- Custom two-boson site type with U(1) × U(1) QN conservation.
- Bose-Hubbard Hamiltonian with A/B species, interspecies coupling, and chemical potentials.
- DMRG driver with reasonable sweep defaults.
- Observables including densities and connected triple correlators (normal-ordered).
- Binding-energy estimators for impurity sectors.
- YAML-based parameter loading.
- Optional HDF5 save/load of ground states.

## Repository Layout
- `src/` — library code (module `BHMimpurityDMRG`)
- `src/observables/` — observables and correlators
- `scripts/` — runnable examples / apps
- `configs/parameters.yaml` — example configuration

## Install Dependencies
This is a Julia project intended to be used in a Julia environment with:
- `ITensors`
- `ITensorMPS`
- `YAML` (for config files)
- `HDF5` (optional, for saving/loading states)

Install as needed:
```julia
import Pkg
Pkg.add(["ITensors", "ITensorMPS", "YAML", "HDF5"])
```

## Usage
Run the example apps:
```bash
julia scripts/impuritydmrg.jl
julia scripts/binding_energy_app.jl
julia scripts/triple_corr_app.jl
```

### YAML Parameters
By default, scripts load `configs/parameters.yaml`. You can point to a different file:
```bash
PARAMS=/path/to/params.yaml julia scripts/triple_corr_app.jl
```

Example `parameters.yaml`:
```yaml
meta:
  description: "Default parameters for 1D impurity Bose-Hubbard DMRG"
  date: "2026-02-06"
  run_by: ""

lattice:
  L: 12
  periodic: true

local_hilbert:
  nmax_a: 2
  nmax_b: 3
  conserve_qns: true

initial_state:
  Na_total: 12
  Nb_total: 3
  impurity_distribution: "center"  # "center" or "random"
  seed: 123

hamiltonian:
  t_a: 1.0
  t_b: 1.0
  U_a: 16.0
  U_b: 0.0
  U_ab: -2.0
  mu_a: 6.627
  mu_b: 4.2

dmrg:
  nsweeps: 15
  cutoff: 1e-10
  maxdim: [50, 100, 200, 400, 600, 800, 800, 800, 800, 800, 800, 800]

observables:
  triple_corr:
    species: "both"   # "a", "b", or "both"
    all_pairs: false  # If true, compute on a grid. If false, use the specified pairs below.
    rmax: 6           # Optional: restrict r in 0..rmax (used if all_pairs=true or rmax/smax set)
    smax: 6           # Optional: restrict s in 0..smax
    precompute: true  # Precompute ⟨n_i⟩ and ⟨n_i n_j⟩ for speed
    pairs:
      - [0, 0]
      - [0, 1]
      - [1, 2]

io:
  save_state: true
  state_save_path: "dmrg_state.h5"
  results_path: "results.h5"
  log_path: "run.log"

```

## Saving and Loading Ground States
```julia
save_state("gs_state.h5", psi; energy=energy, params_path="configs/parameters.yaml", na=na, nb=nb)
st = load_state("gs_state.h5")
psi = st.psi
sites = st.sites
params_hash = st.params_sha256
```

## Notes
- The code is organized for 1D chains now, but lattice utilities are in `src/lattice.jl` to enable future 2D support.
- For large or 2D systems, expect to increase DMRG bond dimensions significantly.

## License
MIT. See `LICENSE`.
