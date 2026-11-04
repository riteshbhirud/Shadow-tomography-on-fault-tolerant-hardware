"""
The circuit side of the measurement oracle: codes, logical bases, surgery and gauging
constructions, the epoch circuit builder with its stabilizer tracker, and the native-operation
catalog. Circuits are a small IR emitted as stim text.

Qubit ids are 0-based, as in the stim text. Python's iteration
orders (insertion-ordered dicts, networkx adjacency, small-int sets, `random.Random(0)`) decide
which detectors and gate orders get emitted, so they are reproduced exactly.
"""
module Circuits

using OrderedCollections: OrderedDict
using Printf: @sprintf
using PythonCall
using QuantumClifford: PauliOperator, Stabilizer, expect, apply!, nqubits, sHadamard, sPhase, sCNOT, sX, sZ
import Graphs
import JSON

using ..FTShadows: pymodule

export PauliDict, Check, Epoch, Experiment, BuildResult, build, stim_text, NOISE_MODELS, DEFAULT_NOISE,
       BBCode, ordered_checks, bb_setup, bb_gauging, local_basis, local_experiment, compile_setting,
       single_patch_readout, multi_patch_surgery, multi_patch_readout, setting_experiment,
       surface_local_experiment, ComposedExperiment, memory_circuit, logical_to_physical, physical_state,
       exact_expectation, bb72_readout_experiment, random_stabilizer_state, code_distance_isd,
       code_distance_stim, css_basis, memory_experiment, stim_search, code_capacity_circuit, REGISTRY, build_setting

include("bits.jl")
include("pyorder.jl")
include("gf2.jl")
include("pauli.jl")
include("ir.jl")
include("tracker.jl")
include("builder.jl")
include("logical.jl")
include("surface.jl")
include("bb.jl")
include("pygraph.jl")
include("gauging.jl")
include("bb_surgery.jl")
include("bb_local.jl")
include("bb_memory.jl")
include("fragments.jl")
include("bravyi.jl")
include("distance.jl")
include("plugins.jl")

end
