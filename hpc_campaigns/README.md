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
   bash hpc_campaigns/slurm/submit_multilauncher.sh /absolute/path/to/runs/<campaign_name> <n_subjobs> <threads_per_run> [partition] [cpus_per_task]
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
- `observables.triple_corr.pairs` is emitted in compact row form (`- [r, s]`).
- `jobfile` contains one command per run:
  - `julia --project=<repo_root> <app_script> <run_dir>/parameters.yaml`
- Paths for state/results/log/checkpoint are written as absolute paths per run.
- `meta.run_name` in each generated config is set to `run_XXXX`.
- This tooling is intentionally outside `src/` to keep the library scheduler-agnostic.
