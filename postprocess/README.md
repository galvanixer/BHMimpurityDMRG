# Postprocess

Campaign-level aggregation utilities live here so the core `src/` DMRG code stays solver-focused.

## Purpose
- Read many per-run `results.h5` files (for example `runs/<campaign>/run_XXXX/results.h5`).
- Aggregate key metadata and observables into a single `all_results.jld2`.
- Keep parser dispatch schema-aware via `/meta/results_schema_id` and `/meta/results_schema_version`.

## Usage
From repository root:

```bash
julia --startup-file=no --project=postprocess postprocess/aggregate_results.jl runs/<campaign_name>
```

Or with the launcher:

```bash
./bin/aggregate_results runs/<campaign_name>
```

Default outputs:

```text
runs/<campaign_name>/Results_<campaign_name>.jld2
runs/<campaign_name>/Summary_<campaign_name>.arrow
runs/<campaign_name>/Summary_<campaign_name>.csv
```

### Choose extraction profile

```bash
./bin/aggregate_results --profile summary_only runs/<campaign_name>
./bin/aggregate_results --profile full runs/<campaign_name> runs/<campaign_name>/Results_<campaign_name>_full.jld2
```

Profiles are defined in `postprocess/extract_profile.yaml`.

## Convergence report
Use `check_convergence.jl` to summarize which runs converged.

It uses this source priority:
1. `results.h5:/diagnostics/dmrg/converged` (authoritative when present)
2. `run.log` fallback markers (for legacy runs missing diagnostics)

Examples:

```bash
./bin/check_convergence runs/<campaign_name>
./bin/check_convergence runs
```

Default outputs:
- Single campaign input: `runs/<campaign_name>/Convergence_<campaign_name>.csv`
- Multi-campaign input: `runs/Convergence_all_campaigns.csv`

The CSV contains tri-state `convergence_status`:
- `converged`
- `not_converged`
- `unknown`

It also includes `run_status` and `evidence` columns to explain classification.
The CSV also includes per-run Hamiltonian parameters:
`t_a`, `t_b`, `U_a`, `U_b`, `U_ab`, `mu_a`, `mu_b`.
It also includes `seed_initial_state` from `initial_state.seed`.
Parameter source priority:
1. `results.h5:/meta/params_yaml`
2. run-directory YAML (for example `parameters.yaml`)
`last_stored_sweep` source priority:
1. `results.h5:/diagnostics/dmrg` sweep metadata
2. `dmrg_state_checkpoint.h5:/meta/checkpoint_sweep`
3. `run.log` checkpoint lines (`Wrote DMRG checkpoint at sweep ...`)

## Output layout
The JLD2 file contains:
- `manifest`: aggregation config, discovered files, schema counts.
- `summary`: one dictionary row per run (good for quick filtering/plot prep).
- `runs`: per-run detailed records (`meta`, selected `observables`, `issues`).

Tabular exports:
- `*_summary.arrow`: typed columnar table for fast filtering/analysis.
- `*_summary.csv`: plain-text table export.

## Schema contract
Core writers stamp:
- `meta/results_schema_id`
- `meta/results_schema_version`
- `meta/results_writer`

This directory only consumes that contract and can keep evolving independently.
