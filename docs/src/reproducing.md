# [Reproducing results](@id reproducing)

Every reported number comes from a Julia script in `scripts/`, one script per result. Each script starts with a `# Run:` line giving its command line and options. Run the scripts from the repository root. The ones that call Python set `JULIA_PYTHONCALL_EXE` and `JULIA_CONDAPKG_BACKEND` to the project venv unless these are set already.

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

Neither results nor datasets are part of the repository. A study whose stored results are in `results/<name>/` writes to `results/<name>_julia/`, for example `results/taskA_julia/`, `results/calibration_julia/` and `results/taskB_julia/`. The dataset scripts write to a `data_julia/` folder next to the stored datasets, and refuse to write over the stored datasets even when that folder is passed with `--out`. A script run therefore never overwrites a stored result.

!!! warning "Quick runs"
    Scripts run with `--quick` write to the same file names as full runs. Move quick outputs away before running a script that reads them.

## Inputs of later steps

Some scripts read the output of an earlier one. They find it with `FTShadows.upstream(repo, parts...)`: if the `_julia` twin of the file exists, they read the twin, and otherwise the stored file. Only single files are swapped, not folders, because a `_julia` folder can hold part of a dataset. This is what lets `run_all.sh` run from a fresh clone: each step reads the output of the earlier step when there is one.

## Determinism

Everything is seeded (see [Seeds](@ref conventions)). Running a script again gives byte-identical files; only provenance fields such as the date, the run time and the code version change.

!!! note
    Some BB72 scripts decode on several worker processes and keep one Relay-BP decoder per worker between chunks. Per-chunk failure counts then depend on which worker got which chunk. Runs with one worker are exact.
