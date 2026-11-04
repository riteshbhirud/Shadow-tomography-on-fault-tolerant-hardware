# Shadow tomography on fault-tolerant hardware

`FTShadows.jl` simulates fault-tolerant logical measurements for logical shadow tomography. You give it a logical stabilizer state and a measurement setting on the surface code or the [[72,12,6]] bivariate bicycle code (BB72). It builds that code's fault-tolerant measurement circuit, runs it with circuit-level noise, and decodes every shot. For each shot it returns:

- the decoded logical outcomes,
- the decoder's confidence (the "gap"),
- the cost of the setting (rounds, extra qubits, circuit layers).

`src/taskb` is a separate part: it works out which logical Clifford groups BB72, the gross code and La-cross codes implement cheaply, and the shadow norms of those ensembles.

The package started as Python code; this is the Julia version, built on [QuantumClifford.jl](https://github.com/QuantumSavory/QuantumClifford.jl).

## Install

You need Julia 1.12, and Python 3.12 for the tools the package calls through [PythonCall](https://github.com/JuliaPy/PythonCall.jl). [uv](https://docs.astral.sh/uv/) is the easiest way to make the Python environment; plain `venv` + `pip` also works.

```bash
uv venv --python 3.12
uv pip install -r python/requirements.txt      # exact versions behind every reported number
export JULIA_PYTHONCALL_EXE=$PWD/.venv/bin/python
export JULIA_CONDAPKG_BACKEND=Null
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

The two `export` lines tell PythonCall to use this venv instead of installing its own Python. The test suite and `scripts/run_all.sh` set them for you.

Two scripts need their own small Julia environments: Task B needs Oscar (`scripts/oscar`), and the plots need CairoMakie (`scripts/plots`). Install them only if you need them:

```bash
julia --project=scripts/oscar -e 'using Pkg; Pkg.instantiate()'
julia --project=scripts/plots -e 'using Pkg; Pkg.instantiate()'
```

### What runs in Python, and why

Most of the code is Julia. These parts still go through PythonCall:

| Python package | used for | why not Julia |
|---|---|---|
| stim | sampling circuits, detector error models, circuit-distance searches | QuantumClifford has no detector error models, and its Pauli-frame sampler mis-samples X-basis measurements in 0.11.8 |
| PyMatching, Relay-BP, ldpc | decoders | Relay-BP only has a Python binding, and LDPCDecoders.jl's BP-OSD takes one prior for all columns where the noise models need one per column |
| numpy (`default_rng`, `corrcoef`, `.npy`/`.npz`) | seeded random draws, file formats | the seeded streams have to be numpy's to give the same states and permutations |
| midout (Gidney, vendored in `python/gidney_inplace_y`) | in-place Y measurement circuits for the surface code | his own generator (one numpy-2 patch); it reproduces his archived circuits byte for byte |
| matplotlib | the Nov-20 postselection plot | keeps that PNG byte-identical; the other plots use CairoMakie |
| pandas | tests only | checks that the Julia CSV writer writes what pandas writes |

## Usage

```julia
using FTShadows

bell = [("XX", 0), ("ZZ", 0)]  # stabilizers of (|00> + |11>)/sqrt(2), sign bit 0 = +1
r = measure("surface", "local:XY", bell; p=1e-3, shots=1000, seed=1, d=3, memory_rounds=3)
r.outcomes[1:5, :], r.confidence[1:5], r.cost
```

`measure("bb72", "local:YZ", bell; p=1e-3, shots=100, rounds=3)` does the same on BB72, using logical qubits 0 and 1. It takes about a minute. `FTShadows.Circuits.REGISTRY` lists the native operations of each code, and `records(r)` gives one record per shot (outcomes, confidence, cost) as plain dictionaries.

## Layout

| path | contents |
|---|---|
| `src/circuits/` | code definitions, surgery and gauging constructions, the circuit builder |
| `src/decode/` | decoder wrappers and the gap (confidence) computation |
| `src/oracle.jl`, `src/measure.jl` | sampling and decoding; `measure`, the main entry point |
| `src/taskb/` | logical Clifford groups and shadow maps; the group orders use Oscar through `ext/FTShadowsOscarExt.jl` |
| `src/datasets/`, `src/bb72/` | the parts behind the scripts: calibrations, datasets, BB72 studies |
| `python/` | the bridged Python pieces: `requirements.txt` and Gidney's vendored `midout` |
| `scripts/` | one script per reported result; `run_all.sh` runs them all |
| `test/` | the test suite, with its stored reference files in `test/reference/` |

## Tests

```bash
FTSHADOWS_TEST=fast julia --project=. -e 'using Pkg; Pkg.test()'   # about 4 min; this is what CI runs
julia --project=. -e 'using Pkg; Pkg.test()'                        # everything, including the slow items
FTSHADOWS_TEST_FILES=decode,oracle julia --project=. -e 'using Pkg; Pkg.test()'   # only some files
```

Most tests compare against stored reference files in `test/reference/`. A few, tagged `:pyref`, compare against a local copy of the original code in `reference/python` (not tracked) and are skipped without it.

Tests that need `results/` or `handover/` (not in the repository) are skipped when those folders are missing.

## Reproducing results

Results go to `results/` and datasets to `handover/`; neither is tracked.

```bash
bash scripts/run_all.sh          # everything (days of CPU time; long steps resume)
bash scripts/run_all.sh --fast   # tests only
julia --project=scripts/oscar scripts/taskb_run_all.jl   # Task B only
```

The Julia scripts write to `results/<name>_julia/` and `handover/data_julia/`, so they never overwrite the stored reference results. When a script needs the output of an earlier one, it reads the Julia file if it is there and the stored file otherwise. Note that `--quick` runs write to the same file names as full runs.

Everything is seeded. Rerunning a script gives byte-identical files; only provenance fields such as the date and runtime change.

## What is checked and what is not

Every script's output was checked bit for bit against the original Python implementation: the same circuits, shots, decoder outputs and output files for the same seeds (or the same numbers, where only timings differ), some at reduced shot counts.

Checked by tests and scripts in this repository:

- Noiseless runs give the exact Born-rule outcomes for every setting on both codes.
- The noise model reproduces Bravyi et al.'s published BB72 memory error rates within 0.7–1.1×.
- Gidney's in-place Y measurement is reproduced byte for byte. With our padding rounds it has circuit distance d at d = 3, 5, 7.
- Every BB72 gauging combination used by the Bell and GHZ settings has circuit distance 6, with all 12 logical qubits as observables.

Not modelled or not verified:

- State preparation is ideal; noise starts at the first noisy round after it.
- Noise is uniform circuit-level depolarizing only, with no leakage, crosstalk or biased noise.
- BB72 gaps marked `not_found` are unknown: the decoder found no correction for the flipped outcome, so the gap is recorded as infinite but the true value could be anything. This is about half of the noisy BB72 shots at r = 3 and over 80% at r = 6.
- We assume the flip rate q does not depend on the input state. This is checked for single measurements, not for every multi-qubit BB72 setting.

## Licence

MIT (`LICENSE`). Gidney's `midout` in `python/gidney_inplace_y` is CC-BY-4.0 (see its `NOTICE.md`).
