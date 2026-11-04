# Decoders returning (predicted observable flips, per-shot confidence, confidence kind).
#
# The confidence is the complementary gap ln[P(best correction) / P(best correction with the other value of
# the observable)] under the DEM's independent-error model; with several observables, the minimum over them.
# MatchingGapDecoder computes it exactly (graphlike DEMs, surface code). The BP-type decoders (BB codes)
# decode under each hypothesis with the observable rows appended to H, so their gaps are approximate.
# stim, PyMatching, Relay-BP and ldpc are bridged; everything around them is here.

include("numerics.jl")
include("dem.jl")
include("matching.jl")
include("hypergraph.jl")
