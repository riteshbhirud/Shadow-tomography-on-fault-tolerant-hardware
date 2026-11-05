# FTShadows.jl

```@meta
DocTestSetup = quote
    using FTShadows
end
```

FTShadows.jl simulates fault-tolerant measurements of logical Pauli operators, for use in logical shadow tomography. Given a logical stabilizer state and a measurement setting on the rotated surface code or on the ``[[72,12,6]]`` bivariate bicycle code (BB72) [bravyi2024high](@cite), it builds the fault-tolerant measurement circuit, samples it with circuit-level noise, and decodes every shot. Each shot comes back with its decoded logical outcomes, a decoder confidence (the "gap") and the cost of the setting. These shot records are the input of classical-shadow estimators [huang2020predicting](@cite). Circuits are sampled with stim ([github.com/quantumlib/stim](https://github.com/quantumlib/stim)), and the package is built on QuantumClifford.jl.

The library consists of a few related components:

- One entry point, [`measure`](@ref), that runs a native setting of a code and returns an [`OracleResult`](@ref). The settings of each code are listed in `FTShadows.Circuits.REGISTRY`; see [Measuring logical settings](@ref manual).

- Decoders that return a per-shot confidence next to the prediction: an exact complementary gap for the surface code ([`MatchingGapDecoder`](@ref)) and an approximate one for BB72 ([`RelayAnchoredGapDecoder`](@ref)). The meaning of each value is given by its kind, see [`CONF_KINDS`](@ref) and [Decoders and the gap](@ref gap).

- Circuit-level noise models with one parameter p, see [Noise models](@ref noise).

- Shot-record files and their loaders, [`load_dataset`](@ref), [`load_records`](@ref) and [`local_view`](@ref), described in [Shot records and datasets](@ref datasets).

- `FTShadows.TaskB`, which computes the logical Clifford groups that BB72, the gross code and La-cross codes implement natively, and the shadow norms of those ensembles. See [Native Clifford ensembles](@ref taskb).

Signs, qubit order and setting ids are collected on the [Conventions](@ref conventions) page.

## Example Usage

Measure logical qubit 0 of a Bell pair in X and logical qubit 1 in Y, on two distance-3 surface-code patches that first idle for 3 noisy rounds. The state is given by its stabilizer generators, as `(Pauli label, sign bit)` tuples.

```@example index
using FTShadows

bell = [("XX", 0), ("ZZ", 0)]   # (|00> + |11>)/√2
r = measure("surface", "local:XY", bell; p = 1e-3, shots = 1000, seed = 1, d = 3, memory_rounds = 3)
r.setting
```

The decoded outcome bits (0 is eigenvalue +1) and the decoder gap of the first shots:

```@example index
r.outcomes[1:4, :]
```

```@example index
r.confidence[1:4]
```

The cost of the setting, the same for every shot:

```@example index
r.cost
```

The product of the two outcomes estimates ``\langle X_0 Y_1 \rangle``, which is 0 for this state.

```@example index
sum((1 .- 2r.outcomes[:, 1]) .* (1 .- 2r.outcomes[:, 2])) / 1000
```

## Installation

FTShadows is not registered. Clone the repository and use its root as the Julia project. Julia 1.12 is needed, and a Python 3.12 environment for the tools called through PythonCall:

```bash
uv venv --python 3.12
uv pip install -r python/requirements.txt
export JULIA_PYTHONCALL_EXE="$PWD/.venv/bin/python" JULIA_CONDAPKG_BACKEND=Null
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

The [Python dependencies](@ref python) page explains these variables and lists what each Python package is used for.

The test suite is run with

```bash
FTSHADOWS_TEST=fast julia --project=. -e 'using Pkg; Pkg.test()'
```

## What is not modelled

!!! warning "Limits of the model"
    - State preparation is ideal; noise starts with the first noisy round after it.
    - The noise is uniform circuit-level depolarizing noise only, with no leakage, crosstalk or biased noise.
    - The surgery and gauging schedules measure X-type and Z-type checks one after the other. They are not optimised, so their flip rates are pessimistic compared with an interleaved schedule.
    - The flip rate of a measured string is assumed not to depend on the input state. This was checked for single measurements, not for every multi-qubit BB72 setting.
    - Task B produces no shot records; its ensembles are evaluated at the logical level only.
