# FTShadows.jl

[![docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://riteshbhirud.github.io/Shadow-tomography-on-fault-tolerant-hardware/)

Simulation of fault-tolerant logical measurements for shadow tomography on the surface code and the [[72,12,6]] bivariate bicycle code. Each shot comes with its decoded outcomes, a decoder confidence (the complementary gap), and the cost of the measurement. Built on [QuantumClifford.jl](https://github.com/QuantumSavory/QuantumClifford.jl), with stim, PyMatching and Relay-BP called through PythonCall.

## Installation

Needs Julia 1.12 and Python 3.12.

```bash
uv venv --python 3.12
uv pip install -r python/requirements.txt
export JULIA_PYTHONCALL_EXE=$PWD/.venv/bin/python JULIA_CONDAPKG_BACKEND=Null
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Usage

```julia
using FTShadows

bell = [("XX", 0), ("ZZ", 0)]
r = measure("surface", "local:XY", bell; p=1e-3, shots=1000, d=3)

r.outcomes      # decoded bits, one row per shot
r.confidence    # decoder gap per shot
r.cost          # rounds, extra qubits, circuit layers
```

See the [documentation](https://riteshbhirud.github.io/Shadow-tomography-on-fault-tolerant-hardware/) for the available settings, decoders, noise models and dataset format.

## Tests

```bash
FTSHADOWS_TEST=fast julia --project=. -e 'using Pkg; Pkg.test()'
```
