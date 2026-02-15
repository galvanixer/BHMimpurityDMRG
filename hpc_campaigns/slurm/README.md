# Slurm Submission (dynamic launcher)

This folder contains Slurm scripts for running campaign `jobfile` commands through an in-repo dynamic launcher.

## Files
- `job.slurm`: Slurm array worker. Each array task calls:
  - `hpc_campaigns/slurm/dynamic_multilauncher.sh <run_root>/jobfile <threads_per_run>`
- `dynamic_multilauncher.sh`: shared-queue launcher with dynamic work stealing across array tasks.
- `submit_multilauncher.sh`: helper that computes array size and submits `job.slurm`.

## Prerequisites
- Bash available on compute nodes.
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
bash hpc_campaigns/slurm/submit_multilauncher.sh /absolute/path/to/runs/<campaign_name> 16 2 grant
```

Submit to `grant` choosing a named account profile:
```bash
bash hpc_campaigns/slurm/submit_multilauncher.sh /absolute/path/to/runs/<campaign_name> 16 2 grant francesco
```

Arguments:
- `run_root` (required): campaign directory containing `jobfile`.
- `n_subjobs` (optional): number of Slurm array tasks. `0` means auto full-node packing (`1` task).
- `threads_per_run` (optional): expected CPU threads per command, used to set per-node parallel slots.
- `partition` (optional): submit-time partition override (`public`, `grant`, `publicgpu`, `grantgpu`, ...).
- `account_name` (optional): selector for account mapping when using `grant`/`grantgpu`.
  - `francesco` -> looks for `CAIUS_GRANT_ACCOUNT_FRANCESCO` (or legacy equivalent)
  - `acct:g2025a457b` -> uses direct account id without lookup
- `job_script` (optional): alternate Slurm script path.

Preferred argument order:
`<run_root> [n_subjobs] [threads_per_run] [partition] [account_name] [job_script]`

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

## Full-Node Scheduling Note
On CAIUS partitions that allocate full nodes per array task, prefer:
- `n_subjobs=0` (auto => 1 task, pack work inside node)
- tune `threads_per_run` to match your per-command threading needs
- increase `n_subjobs` only when you intentionally want multiple full nodes
- dynamic launcher scheduling avoids static chunk imbalance between heterogenous nodes

## Cluster Settings
Edit `job.slurm` header for your cluster:
- partition/account (`#SBATCH -p`, `#SBATCH -A`)
- walltime (`#SBATCH --time`)
- CPU/memory (`#SBATCH --cpus-per-task`, `#SBATCH --mem`)
- output pattern (`#SBATCH --output`)
