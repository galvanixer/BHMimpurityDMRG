# Slurm Submission (dynamic launcher)

This folder contains Slurm scripts for running campaign `jobfile` commands through an in-repo dynamic launcher.

## Files
- `job.slurm`: Slurm array worker. Each array task calls:
  - `hpc_campaigns/slurm/dynamic_multilauncher.sh <run_root>/jobfile <threads_per_run>`
- `dynamic_multilauncher.sh`: shared-queue launcher with dynamic work stealing across array tasks.
- `submit_multilauncher.sh`: helper that computes array size and submits `job.slurm`.
- `precompile_project.sh`: one-shot helper to instantiate/precompile with a persistent depot.

## Prerequisites
- Bash available on compute nodes.
- Campaign already generated with:
  - `hpc_campaigns/launch_campaign.jl`
  - and containing `<run_root>/jobfile`.

## Persistent Julia Cache (recommended)
To avoid repeated precompilation on random nodes, use a persistent shared depot:

```bash
export BHM_JULIA_DEPOT=/shared/path/julia_depot_bhm
export JULIA_CPU_TARGET=generic
export JULIA_PKG_PRECOMPILE_AUTO=0
```

Warm it once:

```bash
bash hpc_campaigns/slurm/precompile_project.sh /path/to/BHMimpurityDMRG /shared/path/julia_depot_bhm
```

`job.slurm` uses:
- `JULIA_PROJECT=$REPO_ROOT` (detected on HPC filesystem)
- `JULIA_DEPOT_PATH=${JULIA_DEPOT_PATH:-${BHM_JULIA_DEPOT:-$HOME/.julia_bhmimpuritydmrg}}`

### How `precompile_project.sh` works
`precompile_project.sh` is a one-shot warmup helper for package/cache setup.

Usage:
```bash
bash hpc_campaigns/slurm/precompile_project.sh [repo_root] [julia_depot]
```

Arguments:
- `repo_root` (optional):
  - default: `${BHM_REPO_ROOT}` if set, else repository inferred from script location.
- `julia_depot` (optional):
  - default: `${BHM_JULIA_DEPOT}` if set, else `${JULIA_DEPOT_PATH}`, else `$HOME/.julia_bhmimpuritydmrg`.

What it does:
1. Validates that `Project.toml` exists under `repo_root`.
2. Exports:
   - `JULIA_PROJECT=repo_root`
   - `JULIA_DEPOT_PATH=julia_depot`
   - `JULIA_CPU_TARGET=${JULIA_CPU_TARGET:-generic}`
   - `JULIA_PKG_PRECOMPILE_AUTO=0`
3. Creates the primary depot directory if needed.
4. Runs:
   - `Pkg.instantiate()`
   - `Pkg.precompile()`

Result:
- Dependencies are installed and precompiled into the persistent depot, so production jobs on random nodes avoid repeated cold precompile.

## Usage
Submit with defaults:
```bash
bash hpc_campaigns/slurm/submit_multilauncher.sh /absolute/path/to/runs/<campaign_name>
```

Submit with explicit array size and threads per run:
```bash
bash hpc_campaigns/slurm/submit_multilauncher.sh /absolute/path/to/runs/<campaign_name> 16 2
```

Submit with explicit threading env overrides:
```bash
OPENBLAS_NUM_THREADS=4 JULIA_NUM_THREADS=1 JULIA_NUM_GC_THREADS=1 OPENBLAS_DYNAMIC=0 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
bash hpc_campaigns/slurm/submit_multilauncher.sh /absolute/path/to/runs/<campaign_name> 0 5 24 public
```

Submit with explicit CPU request per array task:
```bash
bash hpc_campaigns/slurm/submit_multilauncher.sh /absolute/path/to/runs/<campaign_name> 1 1 25 public
```

Submit to `grant` with account credential:
```bash
bash hpc_campaigns/slurm/submit_multilauncher.sh /absolute/path/to/runs/<campaign_name> 16 2 16 grant
```

Submit to `grant` choosing a named account profile:
```bash
bash hpc_campaigns/slurm/submit_multilauncher.sh /absolute/path/to/runs/<campaign_name> 16 2 16 grant francesco
```

Arguments:
- `run_root` (required): campaign directory containing `jobfile`.
- `num_array_tasks` (optional): number of Slurm array tasks.
  - default (or `0`): `ceil(jobfile_lines / 24)`
  - on CAIUS full-node partitions this is effectively the number of nodes.
- `threads_per_run` (optional): expected CPU threads per command, used to compute launcher slots.
- Thread env defaults (can be overridden via environment):
  - `JULIA_NUM_THREADS=1`
  - `OPENBLAS_NUM_THREADS=max(1, threads_per_run-1)`
  - `JULIA_NUM_GC_THREADS=1`
  - `OPENBLAS_DYNAMIC=0`
  - `OMP_NUM_THREADS=1`
  - `MKL_NUM_THREADS=1`
- `cpus_per_task` (optional): Slurm `--cpus-per-task` for each array task (default: `24`).
- `partition` (optional): submit-time partition override (`public`, `grant`, `publicgpu`, `grantgpu`, ...).
  - walltime defaults by partition: `public/publicgpu -> 1-00:00:00`, `grant/grantgpu -> 4-00:00:00`
- `account_name` (optional): selector for account mapping when using `grant`/`grantgpu`.
  - `francesco` -> looks for `CAIUS_GRANT_ACCOUNT_FRANCESCO` (or legacy equivalent)
  - `acct:g2025a457b` -> uses direct account id without lookup
- `job_script` (optional): alternate Slurm script path.

Preferred argument order:
`<run_root> [num_array_tasks] [threads_per_run] [cpus_per_task] [partition] [account_name] [job_script]`

## Default Request Profile
If you run:
```bash
bash hpc_campaigns/slurm/submit_multilauncher.sh <run_root>
```
the effective request is:
- array: `--array=1-K`, where `K=ceil(n_commands/24)` (`num_array_tasks` omitted or `0`)
- partition: `public` (from `job.slurm`)
- time limit: `1-00:00:00` (submit default for `public`; `grant`/`grantgpu` use `4-00:00:00`)
- cpus per task: `24` (submit default)
- memory: `--mem-per-cpu=4G` from `job.slurm` (total per task = `cpus_per_task * 4G`)
- launcher slots per array task: `floor(cpus_per_task / threads_per_run)`
- thread env:
  - `JULIA_NUM_THREADS=1`
  - `OPENBLAS_NUM_THREADS=1` (derived from `threads_per_run=1`)
  - `JULIA_NUM_GC_THREADS=1`
  - `OPENBLAS_DYNAMIC=0`
  - `OMP_NUM_THREADS=1`
  - `MKL_NUM_THREADS=1`

To request a different CPU count explicitly, pass `cpus_per_task`:
```bash
bash hpc_campaigns/slurm/submit_multilauncher.sh <run_root> 0 1 25 public
```
This adds `--cpus-per-task=25` at submit time.

## Threading Environment
`submit_multilauncher.sh` now exports threading vars explicitly to `sbatch`:
- `JULIA_NUM_THREADS`
- `JULIA_NUM_GC_THREADS`
- `OPENBLAS_NUM_THREADS`
- `OPENBLAS_DYNAMIC`
- `OMP_NUM_THREADS`
- `MKL_NUM_THREADS`

Defaults are derived from `threads_per_run` as:
- `JULIA_NUM_THREADS=1`
- `OPENBLAS_NUM_THREADS=max(1, threads_per_run-1)`

The script warns if:
- `JULIA_NUM_THREADS>1` and `OPENBLAS_NUM_THREADS>1` (nested threading risk)
- `threads_per_run` does not match `JULIA_NUM_THREADS + OPENBLAS_NUM_THREADS - 1`

## Partition Time Defaults
`submit_multilauncher.sh` applies submit-time walltime defaults:
- `public`: `1-00:00:00`
- `publicgpu`: `1-00:00:00`
- `grant`: `4-00:00:00`
- `grantgpu`: `4-00:00:00`

For other partitions, no time override is added and `job.slurm` header is used.

## Example: 25 Jobfile Lines
If `<run_root>/jobfile` has 25 lines and you want up to 25 concurrent single-thread runs on one array task:
```bash
bash hpc_campaigns/slurm/submit_multilauncher.sh <run_root> 1 1 25 public
```

What this requests:
- `--array=1-1`
- `--cpus-per-task=25`
- `--partition=public`
- `--mem-per-cpu=4G` (from `job.slurm`, so with `--cpus-per-task=25` total is about `100G`)

What launcher does:
- `threads_per_run=1`
- `cpu_budget=25`
- `slots=floor(25/1)=25`
- up to 25 commands can run concurrently; extra commands queue.

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
- omitted `num_array_tasks` (or `0`) for auto `ceil(jobfile_lines/24)` node count
- set `cpus_per_task` to the number of concurrent CPU threads you want the launcher to use
- tune `threads_per_run` to match your per-command threading needs
- increase `num_array_tasks` only when you intentionally want multiple full nodes
- dynamic launcher scheduling avoids static chunk imbalance between heterogenous nodes

## Cluster Settings
Edit `job.slurm` header for your cluster:
- partition/account (`#SBATCH -p`, `#SBATCH -A`)
- walltime (`#SBATCH --time`)
- CPU/memory (`#SBATCH --cpus-per-task`, `#SBATCH --mem`)
- output pattern (`#SBATCH --output`)
