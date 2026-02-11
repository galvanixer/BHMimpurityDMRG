# HPC Campaigns

This folder contains Slurm-first tooling to generate and submit many DMRG runs.

## Files
- `hpc_campaigns/launch_campaign.jl`: generates run directories and per-run `parameters.yaml`.
- `hpc_campaigns/templates/base_campaign.yaml`: example campaign definition.
- `hpc_campaigns/slurm/job.slurm`: one array-task job runner.
- `hpc_campaigns/slurm/submit_array.sh`: submits all generated runs as a Slurm array.

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
3. Submit array:
   ```bash
   bash hpc_campaigns/slurm/submit_array.sh /absolute/path/to/runs/<campaign_name>
   ```

## Generated Structure
For each run, the launcher creates:
- `runs/<campaign_name>/run_XXXX/parameters.yaml`
- `runs/<campaign_name>/run_XXXX/run.log`
- `runs/<campaign_name>/run_XXXX/results.h5`
- `runs/<campaign_name>/run_XXXX/dmrg_state.h5`
- `runs/<campaign_name>/run_XXXX/dmrg_state_checkpoint.h5`

It also creates:
- `runs/<campaign_name>/index.csv`
- `runs/<campaign_name>/run_dirs.txt`

## Notes
- You can set run-root in 3 ways (priority order):
  1. CLI arg `output_root_override`
  2. `campaign.output_root_abs` in campaign YAML
  3. `campaign.output_root` in campaign YAML
- Paths for state/results/log/checkpoint are written as absolute paths per run.
- `meta.run_name` in each generated config is set to `run_XXXX`.
- This tooling is intentionally outside `src/` to keep the library scheduler-agnostic.
