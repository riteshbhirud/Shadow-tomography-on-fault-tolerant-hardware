# Shadow tomography on fault-tolerant hardware

[![docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://riteshbhirud.github.io/Shadow-tomography-on-fault-tolerant-hardware/)

`FTShadows.jl` simulates fault-tolerant logical measurements for logical shadow tomography. You give it a logical stabilizer state and a measurement setting on the surface code or the [[72,12,6]] bivariate bicycle code (BB72). It builds that code's fault-tolerant measurement circuit, runs it with circuit-level noise, and decodes every shot. For each shot it returns:

- the decoded logical outcomes,
- the decoder's confidence (the "gap"),
- the cost of the setting (rounds, extra qubits, circuit layers).

`src/taskb` is a separate part: it works out which logical Clifford groups BB72, the gross code and La-cross codes implement cheaply, and the shadow norms of those ensembles.

It is built on [QuantumClifford.jl](https://github.com/QuantumSavory/QuantumClifford.jl).

## Install

You need Julia 1.12, and Python 3.12 for the tools the package calls through PythonCall.

```bash
uv venv --python 3.12
uv pip install -r python/requirements.txt
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

