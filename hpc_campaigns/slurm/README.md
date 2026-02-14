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

Submit to `grant` with account credential:
```bash
bash hpc_campaigns/slurm/submit_multilauncher.sh /absolute/path/to/runs/<campaign_name> 16 2 hpc_campaigns/slurm/job.slurm grant
```

Submit to `grant` choosing a named account profile:
```bash
bash hpc_campaigns/slurm/submit_multilauncher.sh /absolute/path/to/runs/<campaign_name> 16 2 hpc_campaigns/slurm/job.slurm grant francesco
```

Arguments:
- `run_root` (required): campaign directory containing `jobfile`.
- `n_subjobs` (optional): number of Slurm array tasks. `0` means auto.
- `threads_per_run` (optional): passed to `hpc_multilauncher`.
- `job_script` (optional): alternate Slurm script path.
- `partition` (optional): submit-time partition override (`public`, `grant`, `publicgpu`, `grantgpu`, ...).
- `account_name` (optional): selector for account mapping when using `grant`/`grantgpu`.
  - `francesco` -> looks for `CAIUS_GRANT_ACCOUNT_FRANCESCO` (or legacy equivalent)
  - `acct:g2025a457b` -> uses direct account id without lookup

## Credentials For `grant` / `grantgpu`
`submit_multilauncher.sh` can auto-attach account credentials when you pass:
- `partition=grant` -> uses `CAIUS_GRANT_ACCOUNT` (or `GRANT_ACCOUNT`)
- `partition=grantgpu` -> uses `CAIUS_GRANTGPU_ACCOUNT` (or `GRANTGPU_ACCOUNT`)

For multi-user setups, user-specific variables are supported:
- `CAIUS_GRANT_ACCOUNT_<USER_KEY>`
- `CAIUS_GRANTGPU_ACCOUNT_<USER_KEY>`

`<USER_KEY>` is your username uppercased, with non `[A-Z0-9]` converted to `_`.
Example: `tanul` -> `TANUL`.

Resolution priority is:
1. If `account_name` is passed:
   - `CAIUS_<PARTITION>_ACCOUNT_<ACCOUNT_NAME_KEY>`
   - `<PARTITION>_ACCOUNT_<ACCOUNT_NAME_KEY>`
2. Current-user keys:
   - `CAIUS_<PARTITION>_ACCOUNT_<USER_KEY>`
   - `<PARTITION>_ACCOUNT_<USER_KEY>`
3. Global keys:
   - `CAIUS_<PARTITION>_ACCOUNT`
   - `<PARTITION>_ACCOUNT`

`<ACCOUNT_NAME_KEY>` is `account_name` uppercased with non `[A-Z0-9]` converted to `_`.

You can set these either:
1. In your shell environment:
```bash
export CAIUS_GRANT_ACCOUNT="your_grant_account"
export CAIUS_GRANTGPU_ACCOUNT="your_grantgpu_account"
# Optional per-user overrides
export CAIUS_GRANT_ACCOUNT_TANUL="tanul_grant_account"
export CAIUS_GRANTGPU_ACCOUNT_TANUL="tanul_grantgpu_account"
# Optional named selectors
export CAIUS_GRANT_ACCOUNT_FRANCESCO="g2025a457b"
```
2. In local file `hpc_campaigns/slurm/credentials.local.sh` (not tracked by git), e.g.:
```bash
GRANT_ACCOUNT="your_grant_account"
GRANTGPU_ACCOUNT="your_grantgpu_account"
# Optional per-user overrides
CAIUS_GRANT_ACCOUNT_TANUL="tanul_grant_account"
CAIUS_GRANTGPU_ACCOUNT_TANUL="tanul_grantgpu_account"
# Optional named selectors
CAIUS_GRANT_ACCOUNT_FRANCESCO="g2025a457b"
```

Template file: `hpc_campaigns/slurm/credentials.local.sh.example`

## Cluster Settings
Edit `job.slurm` header for your cluster:
- partition/account (`#SBATCH -p`, `#SBATCH -A`)
- walltime (`#SBATCH --time`)
- CPU/memory (`#SBATCH --cpus-per-task`, `#SBATCH --mem`)
- output pattern (`#SBATCH --output`)
