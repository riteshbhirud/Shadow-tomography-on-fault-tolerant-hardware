# [Noise models](@id noise)

```@meta
CurrentModule = FTShadows
DocTestSetup = quote
    using FTShadows
end
```

The noise is circuit-level with one parameter p. A noise model gives the probability of each channel as a multiple of p. `Circuits.NOISE_MODELS` holds the models by name, and `Circuits.DEFAULT_NOISE` is the one [`measure`](@ref) uses.

```@example noise
using FTShadows
using FTShadows: Circuits

Circuits.DEFAULT_NOISE
```

```@example noise
for name in ("bravyi_cal", "sd6", "stim_generated")
    println(rpad(name, 16), Circuits.NOISE_MODELS[name])
end
```

The channels are:

| key | channel |
|:--|:--|
| `reset` | a flip after every reset (`X_ERROR` after a Z-basis reset, `Z_ERROR` after an X-basis reset) |
| `meas` | a flip before every measurement, including the destructive readout of data qubits |
| `gate2` | `DEPOLARIZE2` after every two-qubit gate |
| `idle` | `DEPOLARIZE1` on every idle data qubit in every layer, and after every single-qubit gate (for example the H and `SQRT_X` gates of the Y measurement) |
| `round_data` | `DEPOLARIZE1` on every data qubit once per syndrome round |
| `idle_ancilla` | 1 if waiting ancillas also get the `idle` noise, 0 if not |

- `bravyi_cal` (the default) is the channel set of Bravyi et al. [bravyi2024high](@cite). With their BB72 memory circuit and decoder settings it reproduces their memory logical error rates within a factor 0.7 to 1.1.
- `sd6` adds idle noise on waiting ancillas, so it is more pessimistic.
- `stim_generated` replaces the per-layer idle noise by one data-qubit depolarization per round and is more optimistic. The circuits of the surface `local` settings do not support it.

State preparation is ideal in every model: the logical state is prepared without noise, and noise starts with the first noisy round. The meta of a result names its model even when `p = 0`.

## Choosing a model

[`measure`](@ref) has no noise keyword; [`run_oracle`](@ref) takes `noise=`:

```@example noise
bell = [("XX", 0), ("ZZ", 0)]
exp, labels, decoder = Circuits.build_setting("surface", "surgery:ZZ", bell; d = 3, rounds = 3)
for name in ("bravyi_cal", "stim_generated")
    res = run_oracle(exp, 3e-3, 4000; decoder, seed = 5, setting_labels = labels, noise = name, with_gap = false)
    println(rpad(name, 16), "flip rate ", count(==(1), res.outcomes) / 4000)
end
```

## Flip rates

At the logical level a measured string P is read with a flip rate q, so that ``E[(-1)^b] = (1 - 2q) \langle P \rangle``. q depends on the code, the setting, p and the number of rounds, and not, as far as we have checked, on the input state.

```@example noise
for p in (1e-3, 2e-3, 4e-3)
    res = run_oracle(exp, p, 5000; decoder, seed = 7, setting_labels = labels, with_gap = false)
    println("p = ", p, ": q = ", count(==(1), res.outcomes) / 5000)
end
```

!!! warning "Products of strings"
    The flips of different strings in one setting are not always independent. In surface surgery chains the flips of neighbouring products are correlated. A product of strings must then be corrected with its own measured flip rate, ``E[(-1)^{b_1 \oplus b_2}] = (1 - 2q_{12}) \langle P_1 P_2 \rangle``, not with ``(1 - 2q_1)(1 - 2q_2)``. The patches of a surface `local` setting share no qubits after the ideal preparation and are decoded independently; there the product rule is exact.
