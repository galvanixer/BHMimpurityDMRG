# Current DMRG Settings

This note explains the DMRG control parameters currently used in this repository and the simplified schedule syntax that the code now accepts.

As of March 15, 2026, the main runtime configs are:
- `configs/parameters.yaml` for ground-state and triple-correlation runs
- `configs/binding_energy.yaml` for binding-energy runs

These settings are passed into `run_dmrg` by:
- `scripts/solve_ground_state.jl`
- `scripts/triple_corr_app.jl`
- `scripts/binding_energy_app.jl`

## Active Defaults

Both runtime configs currently use:

```yaml
dmrg:
  nsweeps: 15
  maxdim: [50, 100, 200, 400, 600, 800, 800, 800, 800, 800, 800, 800]
  cutoff: 1e-10
  noise: 0.0

  min_sweeps: 8
  energy_tol: 1e-7
  trunc_tol: 1e-7
  patience: 3

  checkpoint_every: 0
  checkpoint_path: "dmrg_state_checkpoint.h5"
  resume_from_checkpoint: true
  resume_mode: "remaining"
  checkpoint_require_hash: true
  checkpoint_save_densities: false
  checkpoint_density_every: 1

  outputlevel: 1
```

## What Actually Runs

`run_dmrg` expands scalars and short schedules to length `nsweeps`.

With the current defaults:
- `nsweeps = 15`
- `maxdim` becomes

```yaml
[50, 100, 200, 400, 600, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800]
```

- `cutoff` is `1e-10` on every sweep
- `noise` is `0.0` on every sweep

So the current strategy is:
- ramp bond dimension during early sweeps
- keep a fixed tight cutoff
- run with no artificial noise
- stop early only if the convergence checks are satisfied

## Early Stopping

The run can stop before sweep 15 because these checks are active:

- `min_sweeps = 8`
- `energy_tol = 1e-7`
- `trunc_tol = 1e-7`
- `patience = 3`

This means:
- no early-stop checks are applied before sweep 8
- starting at sweep 8, the code checks both:
  - `|E_s - E_{s-1}| < 1e-7`
  - maximum truncation error on that sweep `< 1e-7`
- both conditions must hold for 3 consecutive sweeps

The earliest possible early stop is sweep 10.

## Checkpoint and Resume

- `checkpoint_every = 0` disables writing new checkpoints
- `resume_from_checkpoint = true` still allows loading an existing checkpoint file
- `checkpoint_require_hash = true` requires the checkpoint parameter hash to match the current parameter file
- site compatibility is also checked before resuming
- `resume_mode = "remaining"` means completed sweeps are skipped and only the remaining suffix of each schedule is used

Important consequence:

With the current defaults, no new checkpoints are written. Resume only matters if a matching checkpoint file already exists.

## Simplified Schedule Syntax

The schedule API has been simplified on purpose.

The parser now accepts only:
- exact canonical mode names
- exact canonical dictionary keys

Old aliases are no longer supported.

For all three parameters, there are only three supported top-level forms:

### 1. Scalar

Use one value for all sweeps.

Examples:

```yaml
maxdim: 800
cutoff: 1e-10
noise: 0.0
```

### 2. Explicit Vector

Provide one value per sweep.

Examples:

```yaml
maxdim: [50, 100, 200, 400, 800]
cutoff: [1e-8, 1e-8, 1e-9, 1e-10, 1e-10]
noise: [1e-6, 1e-7, 1e-8, 0.0, 0.0]
```

Length rule:
- if the vector is shorter than `nsweeps`, the last value is repeated
- if the vector is longer than `nsweeps`, extra values are ignored

### 3. Generated Dictionary

Use a compact specification that the code expands into a per-sweep schedule.

## `maxdim`

### Canonical Generated Form

```yaml
dmrg:
  maxdim:
    mode: "warmup"
    min: 50
    max: 1200
    sweeps: 6
```

Accepted dictionary keys:
- `mode`
- `min`
- `max`
- `sweeps`

Rules:
- `mode` must be exactly `"warmup"`
- `max` is required
- `min` is optional; if omitted, the code uses `min(50, max)`
- `sweeps` is optional; if omitted, the code uses a default doubling ramp

Behavior:
- with `sweeps`, the code builds a ramp from `min` to `max` over that many sweeps, then repeats `max`
- without `sweeps`, the code uses a doubling ramp such as `50 -> 100 -> 200 -> 400 -> ... -> max`

Examples:

```yaml
dmrg:
  maxdim:
    mode: "warmup"
    min: 50
    max: 800
    sweeps: 5
```

expands to

```yaml
[50, 100, 200, 400, 800]
```

and

```yaml
dmrg:
  maxdim:
    mode: "warmup"
    min: 50
    max: 1200
```

expands to

```yaml
[50, 100, 200, 400, 800, 1200]
```

Constraints:
- all `maxdim` values must be integers `>= 1`

## `cutoff`

### Canonical Generated Form

```yaml
dmrg:
  cutoff:
    mode: "geometric"
    start: 1e-8
    stop: 1e-10
    sweeps: 6
```

Accepted dictionary keys:
- `mode`
- `start`
- `stop`
- `sweeps`

Rules:
- `mode` must be exactly `"geometric"`
- `start` is required
- `stop` is required
- `sweeps` is optional; if omitted, the schedule spans all `nsweeps`

Behavior:
- the schedule is a geometric interpolation from `start` to `stop`
- after the generated part ends, the final value is repeated

Example:

```yaml
dmrg:
  cutoff:
    mode: "geometric"
    start: 1e-8
    stop: 1e-10
    sweeps: 3
```

expands to

```yaml
[1.0e-8, 1.0e-9, 1.0e-10]
```

Constraints:
- all `cutoff` values must be strictly greater than zero

## `noise`

### Canonical Generated Forms

Linear noise ramp:

```yaml
dmrg:
  noise:
    mode: "linear"
    start: 1e-6
    stop: 0.0
    sweeps: 6
```

Geometric noise ramp:

```yaml
dmrg:
  noise:
    mode: "geometric"
    start: 1e-6
    stop: 1e-9
    sweeps: 4
```

Accepted dictionary keys:
- `mode`
- `start`
- `stop`
- `sweeps`

Rules:
- `mode` must be exactly `"linear"` or `"geometric"`
- `start` is required
- `stop` is required
- `sweeps` is optional; if omitted, the schedule spans all `nsweeps`

Behavior:
- `linear` changes by a fixed additive amount each sweep
- `geometric` changes by a fixed multiplicative factor each sweep
- after the generated part ends, the final value is repeated

Examples:

```yaml
dmrg:
  noise:
    mode: "linear"
    start: 1e-6
    stop: 0.0
    sweeps: 4
```

expands to

```yaml
[1.0e-6, 6.67e-7, 3.33e-7, 0.0]
```

and

```yaml
dmrg:
  noise:
    mode: "geometric"
    start: 1e-6
    stop: 1e-9
    sweeps: 4
```

expands to

```yaml
[1.0e-6, 1.0e-7, 1.0e-8, 1.0e-9]
```

Constraints:
- `linear` noise allows `start >= 0` and `stop >= 0`
- `geometric` noise requires `start > 0` and `stop > 0`
- if you want exact zero at the end, use `linear` or an explicit vector

## What Now Breaks

The parser no longer accepts the old alias-heavy forms.

Configs will now fail if they use:
- `mode: "auto"`
- `mode: "automatic"`
- `mode: "warmup"` for `cutoff`
- `mode: "warmup"` for `noise`
- `mode: "logspace"` for `noise`
- `mode: "logspace"` for `cutoff`
- alias keys such as `initial`, `final`, `target`, `from`, `to`, `maximum`
- alias sweep keys such as `warmup_sweeps`, `steps`, `length`
- dict-only convenience keys such as `value`, `values`, or `schedule`

In short:
- `maxdim` dicts must use `mode: "warmup"` with `min`, `max`, `sweeps`
- `cutoff` dicts must use `mode: "geometric"` with `start`, `stop`, `sweeps`
- `noise` dicts must use `mode: "linear"` or `mode: "geometric"` with `start`, `stop`, `sweeps`

## Practical Guidance

Useful defaults when tuning:

1. Increase final `maxdim` if energies or observables look underconverged.
2. Reduce final `maxdim` if the run is too expensive.
3. Use a geometric `cutoff` ramp if you want loose early sweeps and tighter late sweeps.
4. Use linear `noise` if you want the schedule to end at exact zero.
5. Use geometric `noise` only when you want decade-style decay and do not need exact zero.

## Resume Behavior

If a checkpoint is loaded with `resume_mode: "remaining"`:
- the completed prefix of the `maxdim`, `cutoff`, and `noise` schedules is discarded
- only the remaining suffix is applied
- if the checkpoint already reached `nsweeps`, DMRG is skipped entirely

## Update Scheme

`run_dmrg` calls `ITensorMPS.dmrg(H, psi0, sweeps; ...)`, which in the currently installed ITensorMPS implementation performs two-site updates internally.
