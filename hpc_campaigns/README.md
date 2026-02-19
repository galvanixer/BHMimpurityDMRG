# HPC Campaigns

This folder contains HPC tooling to generate and submit many DMRG runs with a dynamic Slurm launcher.

## Files
- `hpc_campaigns/launch_campaign.jl`: generates run directories and per-run `parameters.yaml`.
- `hpc_campaigns/yaml_helpers.jl`: ordered YAML load/write helpers used by campaign tooling.
- `hpc_campaigns/templates/base_campaign.yaml`: example campaign definition.
- `hpc_campaigns/slurm/job.slurm`: Slurm array worker that calls the dynamic launcher.
- `hpc_campaigns/slurm/dynamic_multilauncher.sh`: dynamic queue launcher for `jobfile` commands.
- `hpc_campaigns/slurm/submit_multilauncher.sh`: submits the Slurm array.

## Workflow
1. Edit `hpc_campaigns/templates/base_campaign.yaml`.
2. Generate runs:
   ```bash
   julia hpc_campaigns/launch_campaign.jl hpc_campaigns/templates/base_campaign.yaml
   ```
   To force a custom absolute run root:
   ```bash
   julia hpc_campaigns/launch_campaign.jl hpc_campaigns/templates/base_campaign.yaml /absolute/path/to/runs
   ```
3. Submit with the dynamic launcher:
   ```bash
   bash hpc_campaigns/slurm/submit_multilauncher.sh /absolute/path/to/runs/<campaign_name>
   ```
   Optional arguments:
   ```bash
   bash hpc_campaigns/slurm/submit_multilauncher.sh /absolute/path/to/runs/<campaign_name> [num_array_tasks] [threads_per_run] [cpus_per_task] [partition] [account_name] [job_script]
   ```

## Generated Structure
For each run, the launcher creates:
- `runs/<campaign_name>/run_XXXX/parameters.yaml`
- `runs/<campaign_name>/run_XXXX/run.log`
- `runs/<campaign_name>/run_XXXX/results.h5`
- `runs/<campaign_name>/run_XXXX/dmrg_state.h5`
- `runs/<campaign_name>/run_XXXX/dmrg_state_checkpoint.h5`

It also creates:
- `runs/<campaign_name>/runs.csv`
- `runs/<campaign_name>/jobfile`
- `runs/<campaign_name>/run_dirs.txt`

## Seed Generation
- `sweep.initial_state.seed` supports:
  - explicit list (existing behavior), e.g. `[14061990, 14061991]`
  - `AUTO` (new behavior)
- When set to `AUTO`, define `seed_generation` in campaign YAML:
  - `mode: random|sequential`
  - `count: <integer>|prompt`
  - random mode keys: `min`, `max`, optional `generator_seed`
  - sequential mode keys: `start`, `step`
- If `count: prompt`, `launch_campaign.jl` asks interactively for number of seeds.
- For non-interactive runs, set `count` to an integer.

## Linked Sweep (Non-Cartesian)
- Use `linked_sweep` for variables that must move together by index.
- All arrays in `linked_sweep` must have the same length.
- `sweep` is still Cartesian.
- Final runs are: `Cartesian(sweep) x linked_sweep_rows`.

Example:
```yaml
sweep:
  hamiltonian.mu_b: [3.8, 4.0]

linked_sweep:
  lattice.L: [32, 64]
  initial_state.Na_total: [32, 64]
```

This gives 4 runs total:
- (`mu_b=3.8`, `L=32`, `Na_total=32`)
- (`mu_b=3.8`, `L=64`, `Na_total=64`)
- (`mu_b=4.0`, `L=32`, `Na_total=32`)
- (`mu_b=4.0`, `L=64`, `Na_total=64`)

## Notes
- You can set run-root in 3 ways (priority order):
  1. CLI arg `output_root_override`
  2. `campaign.output_root_abs` in campaign YAML
  3. `campaign.output_root` in campaign YAML
- `base_config` path resolution:
  1. absolute path as-is
  2. relative to campaign YAML directory (if file exists there)
  3. otherwise relative to repository root
- Generated `parameters.yaml` preserves the section/key order of `base_config` when possible.
- `dmrg.maxdim` is emitted in inline YAML list form for readability.
- You can also auto-generate a warmup schedule by setting `dmrg.maxdim` as a dict:
  - `dmrg.maxdim.max` (required), optional `dmrg.maxdim.min`, `dmrg.maxdim.warmup_sweeps`.
- `observables.triple_corr.pairs` is emitted in compact row form (`- [r, s]`).
- `jobfile` contains one command per run:
  - `julia --project=<repo_root> <app_script> <run_dir>/parameters.yaml`
- Paths for state/results/log/checkpoint are written as absolute paths per run.
- `meta.run_name` in each generated config is set to `run_XXXX`.
- If base config has `meta.date: AUTO`, `launch_campaign.jl` replaces it with the
  current launch timestamp (`YYYY-mm-dd HH:MM:SS`) in generated run configs.
- This tooling is intentionally outside `src/` to keep the library scheduler-agnostic.
