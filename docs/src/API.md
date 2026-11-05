# [Full API](@id API)

## FTShadows

Measurement, decoders, results and shot-record loaders.

```@autodocs
Modules = [FTShadows]
Private = false
```

## FTShadows.Circuits

Code constructions, logical bases, surgery and gauging experiments, the circuit builder and the settings registry.

```@autodocs
Modules = [FTShadows.Circuits]
Private = false
```

## FTShadows.TaskB

Logical Clifford groups of BB72, the gross code and La-cross codes, and the shadow maps of their native ensembles. [`group_order`](@ref FTShadows.TaskB.group_order) and [`reduce_generators`](@ref FTShadows.TaskB.reduce_generators) need `using Oscar`.

```@autodocs
Modules = [FTShadows.TaskB]
Private = false
```

## Other documented names

Names that are not exported but are used in this manual.

```@docs
FTShadows.write_records_raw
FTShadows.read_records_raw
FTShadows.read_meta
FTShadows.pyjson
FTShadows.Dem
FTShadows.Circuits.min_weight_codeword
FTShadows.TaskB.ENSEMBLES
FTShadows.TaskB.logical_weights
FTShadows.TaskB.lacross
```
