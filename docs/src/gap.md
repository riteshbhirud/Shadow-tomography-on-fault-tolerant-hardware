# [Decoders and the gap](@id gap)

```@meta
CurrentModule = FTShadows
DocTestSetup = quote
    using FTShadows
end
```

Every decoder in FTShadows returns, next to its prediction, a per-shot confidence: the complementary gap

```math
\Delta = \ln \frac{P(\text{best correction})}{P(\text{best correction that flips the outcome})},
```

computed under the independent-error model of the circuit's detector error model. For a setting with several strings, the confidence of a shot is the minimum of the gaps of its strings; there is no per-string confidence. Larger values are more trustworthy. Decoder confidences of this kind and their use for error mitigation are studied in [dinca2025error](@cite).

## Decoders

[`DECODERS`](@ref) maps a decoder name to its constructor. Each plug-in picks its decoder (see the [registry](@ref manual)); [`run_oracle`](@ref) accepts any of them through `decoder=`.

| name | type | used for | gap |
|:--|:--|:--|:--|
| `"pymatching"` | [`MatchingGapDecoder`](@ref) | surface code | exact |
| `"relay"` | [`RelayAnchoredGapDecoder`](@ref) | BB72 | approximate |
| `"relay_split"` | [`RelayGapDecoder`](@ref) | comparison studies | approximate |
| `"bposd"` | [`BPOSDGapDecoder`](@ref) | cross-checks | approximate |

- `pymatching` decodes with PyMatching [higgott2025sparse](@cite). The gap comes from two more matchings, one per value of the observable, on a graph whose boundary is split by the observable's parity. It is exact for the decomposed (graphlike) model, and the minimum of the two matchings is checked against the plain matching on every call.
- `relay` takes the prediction and its weight from one plain Relay-BP decode [muller2025improved](@cite). For each string, one more Relay-BP decode with the complementary value of the observable appended to the check matrix gives the complementary weight. The gap is the difference, clipped at 0.
- `relay_split` decodes both values of each observable and predicts the lighter one.
- `bposd` is the same hypothesis-split construction with BP-OSD from ldpc [roffe2020decoding](@cite).

Relay-BP decoders carry a random state that cannot be reseeded, so `run_oracle` builds a new one on every call. A call is reproducible as a whole for a given seed.

## Confidence kinds

The kind of each value is in `confidence_kind`; [`CONF_KINDS`](@ref) lists them.

| kind | value | when |
|:--|:--|:--|
| `exact` | finite or `Inf` | `pymatching`: the exact gap. `Inf` means no correction with the flipped outcome exists. 0 is a tie between the two outcomes |
| `approx` | finite | `relay`: an estimate from two decodes. The error can go either way |
| `not_found` | `Inf` | `relay`: the complementary decode found no correction with the flipped outcome. The gap is unknown |
| `decoder_failed` | 0 | `relay`: the plain decode did not converge |
| `unavailable` | `NaN` | no gap was computed (`with_gap=false`), or none exists (surface surgery with `rounds=1`) |
| `noiseless` | `Inf` | `p = 0`; no decoder runs |

!!! warning "`not_found` is not a large gap"
    A `not_found` shot is stored as `Inf`, but its true gap could be anything; it is not a lower bound. On BB72 at p = 1e-3 about half of the noisy shots are `not_found` at 3 gauging rounds, and over 80% at 6 rounds. Treat these shots as a separate class rather than as the most trustworthy ones.

!!! note "A gap of 0"
    With kind `exact`, 0 is a tie. With kind `approx` it is usually a negative estimate clipped to 0: the complementary decode found a lighter correction than the prediction, so the prediction is probably wrong.

Surface surgery with one merged round has no exact gap: the conjugate-basis termination leaves an odd cycle in the matching graph, and every shot is `unavailable`.

```@example gap
using FTShadows

bell = [("XX", 0), ("ZZ", 0)]
unique(measure("surface", "surgery:ZZ", bell; p = 1e-3, shots = 100, d = 3, rounds = 1).confidence_kind)
```

## Gap and flip rate

The flip rate of the shots kept above a cut falls quickly with the cut. For a `surgery:ZZ` setting on the Bell state, a flip is an outcome of 1:

```@example gap
r = measure("surface", "surgery:ZZ", bell; p = 2e-3, shots = 5000, seed = 7, d = 3, rounds = 3)
extrema(r.confidence)
```

```@example gap
for cut in (0, 2.5, 5, 7.5, 10)
    keep = r.confidence .>= cut
    println("cut ", rpad(cut, 4), " keeps ", rpad(count(keep) / 5000, 7),
            "flip rate ", count(==(1), r.outcomes[keep, 1]) / count(keep))
end
```

The range of the gap depends on the setting. Most shots of a local readout without memory rounds share its largest value:

```@example gap
rl = measure("surface", "local:ZZ", bell; p = 2e-3, shots = 5000, seed = 8, d = 3)
extrema(rl.confidence), count(rl.confidence .== maximum(rl.confidence)) / 5000
```

A cut that suits one setting can therefore remove all shots of another.

## Postselection

Cutting on the gap keeps different settings at different rates. After a cut, the kept shots are distributed over the settings as

```math
\pi'(S) = \frac{\text{kept shots of } S}{\text{kept shots}},
```

not as the nominal ``\pi(S)``. The estimator must then use ``\lambda'_P = \sum_S \pi'(S)\, 1[S \text{ measures } P]`` in place of ``\lambda_P``. Using the nominal ``\lambda_P`` biases every estimate whose settings are kept at different rates.

As an example, take an ensemble of `local:XX` and `surgery:ZZ` with probability 1/2 each, and estimate ``\langle ZZ \rangle`` from the kept shots. Only the surgery setting measures ``ZZ``, so the nominal ``\lambda_{ZZ}`` is 1/2.

```@example gap
loc = measure("surface", "local:XX", bell; p = 2e-3, shots = 5000, seed = 9, d = 3)
sur = r                                              # surgery:ZZ from above
cut = 5
kept_loc, kept_sur = count(loc.confidence .>= cut), count(sur.confidence .>= cut)
kept_loc / 5000, kept_sur / 5000
```

```@example gap
signs = [1 - 2Int(b) for (b, c) in zip(sur.outcomes[:, 1], sur.confidence) if c >= cut]
n_kept = kept_loc + kept_sur
λ_nominal = 1 / 2
λ_kept = kept_sur / n_kept
(nominal = sum(signs) / λ_nominal / n_kept, recomputed = sum(signs) / λ_kept / n_kept)
```

The exact value is +1. With the nominal ``\lambda`` the estimate is too small by the ratio of the keep rates; with the recomputed ``\lambda'`` it is ``1 - 2q`` of the kept surgery shots.

!!! note
    Recomputing ``\lambda'`` assumes that the keep rate of a setting does not depend on the input state. This follows from the detectors being independent of the logical state, but it was not measured directly.
