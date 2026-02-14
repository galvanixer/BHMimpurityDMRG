# Slurm Submission (hpc_multilauncher)

This folder contains Slurm scripts for running campaign `jobfile` commands through `hpc_multilauncher`.

## Files
- `job.slurm`: Slurm array worker. Each array task calls:
  - `hpc_multilauncher <run_root>/jobfile <threads_per_run>`
- `submit_multilauncher.sh`: helper that computes array size and submits `job.slurm`.

## Prerequisites
- `hpc_multilauncher` available in `PATH` on compute nodes.
- `GNU parallel` available (required by `hpc_multilauncher`).
- Campaign already generated with:
  - `hpc_campaigns/launch_campaign.jl`
  - and containing `<run_root>/jobfile`.

## Usage
Submit with defaults:
```bash
bash hpc_campaigns/slurm/submit_multilauncher.sh /absolute/path/to/runs/<campaign_name>
```

Submit with explicit array size and threads per run:
```bash
bash hpc_campaigns/slurm/submit_multilauncher.sh /absolute/path/to/runs/<campaign_name> 16 2
```

Arguments:
- `run_root` (required): campaign directory containing `jobfile`.
- `n_subjobs` (optional): number of Slurm array tasks. `0` means auto.
- `threads_per_run` (optional): passed to `hpc_multilauncher`.
- `job_script` (optional): alternate Slurm script path.

## Cluster Settings
Edit `job.slurm` header for your cluster:
- partition/account (`#SBATCH -p`, `#SBATCH -A`)
- walltime (`#SBATCH --time`)
- CPU/memory (`#SBATCH --cpus-per-task`, `#SBATCH --mem`)
- output pattern (`#SBATCH --output`)
