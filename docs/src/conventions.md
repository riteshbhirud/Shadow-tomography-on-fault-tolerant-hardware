# [Conventions](@id conventions)

```@meta
CurrentModule = FTShadows
DocTestSetup = quote
    using FTShadows
end
```

## Qubit order and Pauli strings

- Character i of every Pauli string acts on logical qubit i, counted from 0. There is no reversal. In Julia, logical qubit i is `label[i+1]`.
- Strings are tensor products of I, X, Y, Z with ``Y = iXZ`` on each qubit. Every string is Hermitian, with no extra phase.
- The strings of a setting have one character per logical qubit. `local:XY` on two qubits measures `["XI", "IY"]`.

## Signs

There are three sign conventions.

| where | form | meaning |
|:--|:--|:--|
| input states of [`measure`](@ref) and of the dataset generators | `(label, sign_bit)` | sign bit 0 is +1, 1 is −1 |
| `outcomes` of results and records | bit b | the measured eigenvalue is ``(-1)^b`` |
| `state_for_testing` and `known_answers` in `meta.json` | `"sign"` | ±1 directly, not a bit |

[`local_view`](@ref) converts outcome bits to ±1.

## [Setting ids and parameters](@id setting-ids)

A setting id is `<operation>:<argument>`. [`measure`](@ref) passes its extra keyword arguments to the operation.

| code | setting id | what is measured | parameters (default) | decoder |
|:--|:--|:--|:--|:--|
| `"surface"` | `local:B0B1…`, ``B \in \{X,Y,Z\}`` | patch i in basis ``B_i`` after `memory_rounds` noisy rounds; X and Z by destructive readout, Y in place | `d` (3), `memory_rounds` (0), `boundary_rounds` (see below) | `pymatching` |
| `"surface"` | `surgery:ZZ`, `surgery:XX` | ``Z_iZ_{i+1}`` (or ``X_iX_{i+1}``) for each neighbouring pair, by lattice surgery with conjugate-basis termination | `d` (3), `n_patches` (2), `rounds` (3) | `pymatching` |
| `"bb72"` | `readout:X`, `readout:Z` | all 12 logical qubits in X (or Z), oracle basis | none | `relay` |
| `"bb72"` | `gauge:X0`, `gauge:Y0`, `gauge:Z0` | ``\bar P_0`` by gauging measurement, oracle basis | `rounds` (3) | `relay` |
| `"bb72"` | `local:B0B1…`, at most 3 letters | logical qubits ``0, \dots, m-1`` of the local basis, by readout plus simultaneous gauging | `rounds` (3) | `relay` |

`boundary_rounds` is the number b of padding rounds of Gidney's Y measurement [gidney2024inplace](@cite); the Y measurement then takes ``1 + b + 1`` rounds. The default is 1, 2 and 4 at d = 3, 5 and 7, the smallest values that keep circuit distance d.

!!! warning "Y padding at other distances"
    The default `boundary_rounds` was checked to give circuit distance d only at d = 3, 5 and 7. For any other d the default is `d ÷ 2`, which has not been checked.

## Logical bases

- Surface `local` settings use Gidney's patch: ``\bar X`` is X on the x = 0 column and ``\bar Z`` is Z on the y = 0 row. The surgery settings use their own patch layout, with ``\bar X`` on the top row and ``\bar Z`` on the left column. Logical-level results do not depend on the layout.
- BB72 has two logical bases. The oracle basis, used by `readout` and `gauge`, is one CSS basis in which logical qubit 0 is a pair of weight-6 logicals that overlap on one qubit, so the measured operator equals ``i \bar X_0 \bar Z_0 = \bar Y_0`` exactly. The local basis, used by `local`, has logical qubits 0, 1 and 2 on disjoint supports, so that they can be gauged at the same time. `Circuits.bb_setup()` and `Circuits.local_basis()` return them.

## Confidence

The confidence of a shot is the complementary gap, minimised over the setting's strings. Its kinds are listed on [Decoders and the gap](@ref gap). In dataset files `Inf` is written as `"inf"`, `NaN` as `null`, and finite values are rounded to 6 decimals.

## Costs

See [Costs](@ref costs). For surface `local` settings, rounds and circuit layers are the maximum over the patches, and qubit counts are the sum over the patches.

## Seeds

- `measure(...; seed)` (default 0) seeds stim's sampler. The same circuit, p and seed give the same shots. Relay-BP decoders are rebuilt on every call, so a BB72 call is reproducible as a whole.
- The scripts and generators derive seeds with [`seed_for`](@ref): the SHA-256 of the parts joined by `|`, with the first 4 bytes read as a little-endian `UInt32`. Every configuration and chunk gets its own seed. The parts are written as Python `repr`s, so that strings, numbers and tuples have one fixed text form.

```@example conventions
using FTShadows
seed_for("nov20_draws", "bb72_bell_random_p1e-3_r3")
```

- Random setting draws and logical-noise flags come from numpy's `default_rng`, through the bridge.
- The shots of one setting of a local-setting dataset are simulated in chunks (10,000 shots on the surface code, 100 on BB72), each with `seed_for("nov20", name, bases, flagged, chunk)`.
