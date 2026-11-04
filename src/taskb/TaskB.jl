"""
Task B: the logical Clifford groups that BB72, the gross code and La-cross implement cheaply, and the shadow
maps of the resulting native ensembles. Conventions: Pauli rows (x|z), Cliffords act from the right, v -> v F.

Group orders and the exact generator reduction need `using Oscar` (extension `FTShadowsOscarExt`).
"""
module TaskB

using Base64: base64decode
using LinearAlgebra: I, dot
using OrderedCollections: OrderedDict
using Random: MersenneTwister
using ..FTShadows: pyjson, pyfloat, pymodule
using Printf: Printf
using PythonCall: Py, pyconvert
using QuantumClifford
import JSON

export StabilizerCode, make_bb_code, make_lacross_code, logical_action, build_ensembles, orbits,
       channel_for_group, group_order, reduce_generators, shadow_map

include("gf2.jl")
include("codes.jl")
include("physops.jl")
include("ensembles.jl")
include("groups.jl")
include("shadow.jl")
include("designed_cover.jl")
include("basis_independence.jl")
include("lacross.jl")
include("build_groups.jl")
include("bruteforce.jl")
include("outputs.jl")

end
