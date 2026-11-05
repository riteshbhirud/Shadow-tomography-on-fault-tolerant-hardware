# [Reproducing results](@id reproducing)

Each result comes from a script in `scripts/`; the header of each script gives its command line. Run them from the repository root.

## `run_all.sh`

```bash
bash scripts/run_all.sh          # everything; days of CPU time
bash scripts/run_all.sh --fast   # the test suite only
```

`run_all.sh` runs the test suite first, then the scripts in this order:

1. the calibration of the noise model against the BB72 memory results of [bravyi2024high](@cite);
2. the surface-code studies;
3. the BB72 studies: memory, circuit distances, gauging calibration, gap decoders, postselection, multi-qubit and small-graph gauging;
4. the surface Y measurement [gidney2024inplace](@cite), the BB72 local settings and their costs;
5. the end-to-end check and the plots;
6. the datasets, calibration tables and logical bases;
7. Task B, then the native-ensemble tables built from it.

Steps marked `~h` in the file take hours. They resume where they stopped when the same command is run again. The stim distance searches need up to about 10 GB of memory each, so run them one at a time.

Task B alone:

```bash
julia --project=scripts/oscar scripts/taskb_run_all.jl
```

The plots use their own environment, `scripts/plots`, which has CairoMakie:

```bash
julia --project=scripts/plots -e 'using Pkg; Pkg.instantiate()'
```

## Where the output goes

Results and datasets are not in the repository. Scripts write to `results/<name>_julia/` (e.g. `results/taskA_julia/`), and the dataset scripts to a `data_julia/` folder next to the stored datasets. They never overwrite stored results.

!!! warning "Quick runs"
    Scripts run with `--quick` write to the same file names as full runs. Move quick outputs away before running a script that reads them.

## Inputs of later steps

Some scripts read the output of an earlier one through `FTShadows.upstream(repo, parts...)`, which returns the file in the `_julia` folder if it exists and the stored file otherwise.

## Determinism

Everything is seeded (see [Seeds](@ref conventions)), so rerunning a script gives the same files up to the date and run-time fields.

!!! note
    Some BB72 scripts decode on several worker processes and keep one Relay-BP decoder per worker between chunks. Per-chunk failure counts then depend on which worker got which chunk. Runs with one worker are exact.
