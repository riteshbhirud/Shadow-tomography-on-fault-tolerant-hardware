# Task A studies of BB72 gauging: the reusable parts behind scripts/taskA_bb72_*.jl and taskA_bb_*.jl.
# Outputs keep the reference's file formats.

using OrderedCollections: OrderedDict
using JSON

include("helpers.jl")
include("distances.jl")
include("sampling.jl")
