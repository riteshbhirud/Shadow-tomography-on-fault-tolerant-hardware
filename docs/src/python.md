# [Python dependencies](@id python)

Most of FTShadows is Julia. A few tools without a Julia equivalent are called through [PythonCall](https://github.com/JuliaPy/PythonCall.jl). The modules are imported the first time a function needs them, not at `using FTShadows`.

| package | used for | why it is Python |
|:--|:--|:--|
| [stim](https://github.com/quantumlib/stim) | sampling circuits, detector error models, circuit-distance searches | QuantumClifford.jl has no detector error models, and its Pauli-frame sampler mis-samples X-basis measurements in version 0.11.8 |
| PyMatching [higgott2025sparse](@cite) | surface-code decoding and its exact gap | no Julia package with the same matching |
| Relay-BP (`relay_bp`) [muller2025improved](@cite) | BB72 decoding | the decoder is written in Rust and has only a Python binding |
| ldpc [roffe2020decoding](@cite) | BP-OSD, a secondary decoder for cross-checks | LDPCDecoders.jl's BP-OSD takes one prior for all columns; the detector error models need one per column |
| numpy | seeded random draws (`default_rng`), `corrcoef`, `.npy` and `.npz` files | the seeded streams of the dataset generators are numpy's |
| scipy | sparse matrices handed to ldpc and Relay-BP | the input format of those decoders |
| `midout` and sinter | Gidney's in-place Y measurement circuits [gidney2024inplace](@cite) | it is his own generator, vendored in `python/gidney_inplace_y` |
| matplotlib | one plot (the postselection figure) | the other plots use CairoMakie |
| pandas | tests only | checks that the CSV writer writes what pandas writes |

`midout` is copied from the data archive of [gidney2023inplacedata](@cite) (CC-BY-4.0, see `python/gidney_inplace_y/NOTICE.md`), with one patch for numpy 2. With stim 1.16 it reproduces his archived noisy circuits byte for byte. FTShadows takes his ideal circuit, replaces his initialisation by the ideal preparation of the logical state, and applies its own noise model.

## Setting up the environment

`python/requirements.txt` pins every package to the versions behind the reported results. [uv](https://docs.astral.sh/uv/) is the simplest way to make the environment; `venv` and `pip` also work. From the repository root:

```bash
uv venv --python 3.12
uv pip install -r python/requirements.txt
export JULIA_PYTHONCALL_EXE="$PWD/.venv/bin/python"
export JULIA_CONDAPKG_BACKEND=Null
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

`JULIA_PYTHONCALL_EXE` tells PythonCall to use this interpreter, and `JULIA_CONDAPKG_BACKEND=Null` stops CondaPkg from installing a Python of its own. CondaPkg would resolve other numpy and scipy versions than the pinned ones, and some outputs depend on the exact versions, for example numpy's random streams.

The test suite, the scripts in `scripts/` that call Python, and `docs/make.jl` set both variables to `.venv/bin/python` when they are not set. In an interactive session, set them before `using FTShadows`.
