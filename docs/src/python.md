# [Python dependencies](@id python)

Most of FTShadows is Julia. A few tools without a Julia equivalent are called through [PythonCall](https://github.com/JuliaPy/PythonCall.jl). The modules are imported the first time a function needs them, not at `using FTShadows`.

| package | used for |
|:--|:--|
| [stim](https://github.com/quantumlib/stim) | sampling, detector error models, circuit-distance searches |
| PyMatching [higgott2025sparse](@cite) | surface-code decoding and its exact gap |
| Relay-BP (`relay_bp`) [muller2025improved](@cite) | BB72 decoding |
| ldpc [roffe2020decoding](@cite) | BP-OSD, used for cross-checks |
| numpy, scipy | seeded draws (`default_rng`), `.npy`/`.npz` files, sparse matrices for the decoders |
| `midout`, sinter | Gidney's in-place Y measurement circuits [gidney2024inplace](@cite) |
| matplotlib | the postselection figure |
| pandas | tests only |

We use stim rather than QuantumClifford.jl for sampling because detector error models are needed for decoding, and QuantumClifford's Pauli-frame sampler mis-samples X-basis measurements in version 0.11.8.

`midout` is copied from the data archive of [gidney2023inplacedata](@cite) (CC-BY-4.0, see `python/gidney_inplace_y/NOTICE.md`), with one patch for numpy 2. We take his ideal circuit, replace the initialisation by an ideal preparation of the logical state, and apply our own noise model.

## Setting up the environment

`python/requirements.txt` pins the versions used for the results. With [uv](https://docs.astral.sh/uv/) (plain `venv` and `pip` work too), from the repository root:

```bash
uv venv --python 3.12
uv pip install -r python/requirements.txt
export JULIA_PYTHONCALL_EXE="$PWD/.venv/bin/python"
export JULIA_CONDAPKG_BACKEND=Null
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

`JULIA_PYTHONCALL_EXE` points PythonCall to this interpreter, and `JULIA_CONDAPKG_BACKEND=Null` stops CondaPkg from installing its own Python, which would pick different numpy versions (numpy's random streams depend on the version).

The test suite, the scripts in `scripts/` that call Python, and `docs/make.jl` set both variables to `.venv/bin/python` when they are not set. In an interactive session, set them before `using FTShadows`.
