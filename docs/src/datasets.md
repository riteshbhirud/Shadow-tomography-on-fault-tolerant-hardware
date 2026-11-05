# [Shot records and datasets](@id datasets)

```@meta
CurrentModule = FTShadows
DocTestSetup = quote
    using FTShadows
end
```

A dataset is a set of files with a common prefix. It holds the shots of one ensemble of settings on one prepared state.

| file | content | read with |
|:--|:--|:--|
| `<name>.meta.json` | run-level metadata: code, noise, ensemble, seeds, and for test datasets the prepared state and known answers | `FTShadows.read_meta`, or `load_dataset(prefix).meta` |
| `<name>.records.jsonl` | one shot per line, in draw order | [`load_records`](@ref) |
| `<name>.exact.csv` | every logical Pauli the ensemble can estimate, with its exact value, ``\lambda_P`` and the settings that measure it | [`load_exact`](@ref) |
| `<name>.expected.csv` | local-setting datasets only: exact, after-logical-noise and noise-expected values | any CSV reader |

[`load_dataset`](@ref) reads the first three and returns a named tuple `(; meta, records, exact)`. It does not read `expected.csv`.

The examples on this page use two small noiseless datasets that ship with the test suite.

```@example datasets
using FTShadows

dir = joinpath(pkgdir(FTShadows), "test", "reference", "records")
ds = load_dataset(joinpath(dir, "surface_d3_bell_p0"))
length(ds.records), length(ds.exact)
```

## The record format

One line of `records.jsonl` is one shot:

```@example datasets
println(first(eachline(joinpath(dir, "surface_d3_bell_p0.records.jsonl"))))
```

| field | meaning |
|:--|:--|
| `shot` | index in draw order, from 0 |
| `setting_id` | the setting drawn for this shot, a key of `meta["ensemble"]` |
| `setting` | the commuting logical Pauli strings measured jointly on this shot |
| `outcomes` | one bit per string; bit b means eigenvalue ``(-1)^b`` |
| `confidence` | the decoder gap: a number, `"inf"` or `null` |
| `confidence_kind` | what the confidence means; see [Decoders and the gap](@ref gap) |
| `cost` | the cost of the setting, the same for every shot of it; see [Costs](@ref costs) |

In the file, finite confidences are rounded to 6 decimals, `Inf` is the string `"inf"` and `NaN` is `null`. [`load_records`](@ref) returns `Inf` and `NaN`, and checks that the shots are in draw order.

```@example datasets
ds.records[2]
```

## Reading records

The ensemble of a dataset and the probability of each setting are in its meta file:

```@example datasets
[k => v["probability"] for (k, v) in ds.meta["ensemble"]]
```

A setting with several strings also measures every product of a subset of them. [`subset_index`](@ref) lists these products:

```@example datasets
subset_index(["ZZI", "IZZ"])
```

It refuses settings in which one qubit carries two different letters, since then a product picks up a phase (``XX \cdot ZZ = -YY``). No setting of the current settings does this.

`exact.csv` gives ``\lambda_P``, the probability that a shot measures P:

```@example datasets
ds.exact
```

With these, the standard estimator for an ensemble of Pauli settings [huang2020predicting, bu2022classical](@cite) is a few lines. The single-shot estimate of ``\langle P \rangle`` is ``1[P \text{ measured}] \prod_j (-1)^{b_j} / \lambda_P``:

```@example datasets
function estimate(records, row)
    v = map(records) do r
        idx = get(subset_index(r.setting), row.pauli, nothing)
        idx === nothing ? 0.0 : (-1)^sum(r.outcomes[idx]) / row.lambda
    end
    sum(v) / length(v)
end

[(row.pauli, row.exact, round(estimate(ds.records, row); digits = 3)) for row in ds.exact]
```

The dataset is noiseless, so the deviations are statistical only. With noise, each estimate is ``(1 - 2q)`` times the exact value, with q the flip rate of the string or product ([Noise models](@ref noise)).

## Local settings

For datasets of local Pauli settings, `local:B0B1…`, [`local_view`](@ref) gives one [`LocalShot`](@ref) per shot, with a basis letter and a ``\pm 1`` outcome per logical qubit. It refuses records that are not local settings.

```@example datasets
bb = load_dataset(joinpath(dir, "bb72_bell_fixed_p0"))
shots = local_view(bb.records)
shots[1].bases, shots[1].outcomes, shots[1].gap_kind
```

This dataset has the fixed-setting layout: the shots come in blocks, one per setting, and the probability ``1/3^m`` of each setting is its exact share of the shots.

```@example datasets
bb.meta["sampling"], bb.meta["shots_per_setting"]
```

The Bell fidelity is ``F = (1 + \langle XX \rangle - \langle YY \rangle + \langle ZZ \rangle)/4``, and ``\lambda = 1/9`` for every weight-2 Pauli:

```@example datasets
term(s, P) = s.bases == collect(P) ? 9 * prod(s.outcomes) : 0
v = [(1 + term(s, "XX") - term(s, "YY") + term(s, "ZZ")) / 4 for s in shots]
sum(v) / length(v)
```

The meta file of a local-setting dataset records the known answers:

```@example datasets
ka = bb.meta["known_answers"]
ka["fidelity_ideal"], ka["fidelity_expected_noisy"]
```

`expected.csv` holds, for every Pauli, the value it should have under the calibrated flip rates, with columns `pauli`, `weight`, `exact`, `exact_with_logical_noise`, `expected_noisy`, `expected_noisy_se`, `lambda_P`, `settings` and `q_mean`:

```@example datasets
foreach(println, first(eachline(joinpath(dir, "bb72_bell_fixed_p0.expected.csv")), 3))
```

## Writing records

[`records`](@ref) does not produce the dataset format (see the warning in [Measuring logical settings](@ref manual)). To make files that [`load_records`](@ref) and [`local_view`](@ref) read, write the dataset fields yourself:

```@example datasets
using FTShadows: write_records_raw
using OrderedCollections: OrderedDict

function dataset_records(res, setting_id; first_shot = 0)
    map(1:size(res.outcomes, 1)) do s
        c = res.confidence[s]
        OrderedDict{String,Any}("shot" => first_shot + s - 1, "setting_id" => setting_id,
            "setting" => res.setting, "outcomes" => Int.(res.outcomes[s, :]),
            "confidence" => isnan(c) ? nothing : isinf(c) ? "inf" : round(c; digits = 6),
            "confidence_kind" => res.confidence_kind[s], "cost" => res.cost)
    end
end

bell = [("XX", 0), ("ZZ", 0)]
rxx = measure("surface", "local:XX", bell; p = 1e-3, shots = 300, seed = 11, d = 3)
ryy = measure("surface", "local:YY", bell; p = 1e-3, shots = 300, seed = 12, d = 3)
path = joinpath(mktempdir(), "mine.records.jsonl")
write_records_raw(path, [dataset_records(rxx, "local:XX"); dataset_records(ryy, "local:YY"; first_shot = 300)])

mine = local_view(load_records(path))
[B => sum(prod(s.outcomes) for s in mine if String(s.bases) == B) / 300 for B in ("XX", "YY")]
```

Use a different seed for every call: two calls with the same seed and the same circuit give the same noise realisation.
