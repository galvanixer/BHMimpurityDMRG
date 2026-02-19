# TODO

## Backlog

- [ ] Make `scripts/sample_configurations.jl` independent of external `parameters.yaml` by default.
  - Use `st.params_yaml` from the loaded state for metadata/reproducibility.
  - Keep external `parameters.yaml` optional (only for explicit override use-cases).
  - Preserve current positional argument order:
    1. `parameters.yaml`
    2. `observables.yaml`
    3. `state_path`
    4. `results_path`
