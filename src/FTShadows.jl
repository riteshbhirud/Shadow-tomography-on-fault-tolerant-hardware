"""
FTShadows: fault-tolerant logical measurements for logical shadow tomography,
built on QuantumClifford.jl. Python tools without a Julia equivalent (stim,
PyMatching, Relay-BP, ldpc, Gidney's `midout`) are called through PythonCall.
"""
module FTShadows

using QuantumClifford
using PythonCall

include("bridge.jl")
include("records.jl")
include("circuits/Circuits.jl")
include("seeds.jl")
include("stats.jl")
include("decode/Decode.jl")
include("oracle.jl")
include("measure.jl")
include("taskb/TaskB.jl")

export seed_for, wilson, two_proportion_test, holm
export MatchingGapDecoder, BPOSDGapDecoder, RelayGapDecoder, RelayAnchoredGapDecoder, decode, CONF_KINDS, DECODERS
export OracleCircuit, OracleResult, PrebuiltExperiment, oracle_circuit, load_prebuilt, run_oracle, records
export measure
export Record, ExactRow, LocalShot, load_dataset, load_records, load_exact, subset_index, local_view

end
