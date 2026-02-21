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

Default output path:

```text
runs/<campaign_name>/all_results.jld2
```

### Choose extraction profile

```bash
./bin/aggregate_results --profile summary_only runs/<campaign_name>
./bin/aggregate_results --profile full runs/<campaign_name> runs/<campaign_name>/all_results_full.jld2
```

Profiles are defined in `postprocess/extract_profile.yaml`.

## Output layout
The JLD2 file contains:
- `manifest`: aggregation config, discovered files, schema counts.
- `summary`: one dictionary row per run (good for quick filtering/plot prep).
- `runs`: per-run detailed records (`meta`, selected `observables`, `issues`).

## Schema contract
Core writers stamp:
- `meta/results_schema_id`
- `meta/results_schema_version`
- `meta/results_writer`

This directory only consumes that contract and can keep evolving independently.
