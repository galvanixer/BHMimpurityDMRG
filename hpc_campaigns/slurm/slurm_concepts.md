# Slurm Concepts For This Project

This document explains the Slurm terms used by the scripts in this folder.

## 0) Node, CPU, Core, Thread, Task (and relationship)

- **Node**: one physical compute machine in the cluster.
- **CPU (socket/chip)**: a processor package inside a node. A node can have 1 or more CPUs.
- **Core**: an independent compute unit inside a CPU.
- **Thread (hardware thread)**: an execution context on a core (with SMT/Hyper-Threading, one core can expose multiple threads).
- **Task (Slurm task, in this project)**: one launched Slurm process context. For array jobs here, one array index = one task running `job.slurm`.

Typical hierarchy:

`node -> cpu(s) -> core(s) -> hardware thread(s)`

Slurm allocates CPU resources to a task (for example via `--cpus-per-task`).  
Inside each task, `dynamic_multilauncher.sh` runs worker slots and pulls commands from a shared queue.

## 1) What is a Slurm "task" here?

In this project, a **task** is one Slurm array element (`SLURM_ARRAY_TASK_ID`).

If you submit:
```bash
sbatch --array=1-N hpc_campaigns/slurm/job.slurm <run_root> <threads_per_run>
```
then Slurm launches `N` tasks.  
Each task runs `job.slurm` once.

Inside each task, the dynamic launcher takes the next pending command from `<run_root>/jobfile` until the queue is empty.

## 2) What does `--cpus-per-task` mean?

Current submit defaults request `--cpus-per-task=24` (24-core standard node assumption).

The dynamic launcher uses available cores (via Slurm env vars) to decide how many commands to run in parallel:
- first cap: `SLURM_CPUS_PER_TASK` when set
- fallback: `SLURM_CPUS_ON_NODE`
- parallel jobs per task ≈ `cpu_budget / threads_per_run`

So increasing `--cpus-per-task` can increase per-task parallel throughput.

## 3) What does memory request mean?

Memory requests are RAM requests, not disk.

Common forms:
- `#SBATCH --mem=16G` : total RAM per task.
- `#SBATCH --mem-per-cpu=4000` : RAM per allocated CPU core (typically MB).

If memory is too low, Slurm may kill the task for OOM.  
If memory is too high, scheduling can be slower.

For the dynamic launcher, memory should cover all commands running concurrently inside one task.

## 4) What does `-N 1 --exclusive` mean?

- `-N 1`: request exactly 1 node.
- `--exclusive`: no other jobs share that node.

This gives your job full-node ownership (all node cores/memory), useful for dense internal packing with the dynamic launcher.

## 5) Current project script vs classic full-node style

Current `job.slurm` is portable and generic:
- uses `--cpus-per-task` and `--mem`
- takes `<run_root>` and optional `<threads_per_run>`
- validates prerequisites (`jobfile`, dynamic launcher script)

A cluster-specific full-node style often looks like:
```bash
#SBATCH -N 1 --exclusive
#SBATCH --mem-per-cpu=4000
#SBATCH -p <partition> -A <account>
#SBATCH -t 4-00:00:00
```

Use whichever matches your cluster policy and workload shape.

## 6) Practical tuning knobs

- `num_array_tasks` (submit script arg): number of array tasks (effectively node count on CAIUS full-node partitions); default is `ceil(jobfile_lines/24)`.
- `threads_per_run` (submit script arg): cores needed by one command.
- `cpus_per_task` (submit script arg): cores available inside each array task; default is `24`.

These three control throughput and memory pressure together.

## 7) Rule of thumb

Choose settings so that:
- commands do not OOM
- node CPU stays busy
- queue wait remains acceptable

Then adjust `num_array_tasks` and `threads_per_run` based on observed runtime and memory usage.
