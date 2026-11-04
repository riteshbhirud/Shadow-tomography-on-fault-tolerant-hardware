# Export the logical bases used by the handover as explicit physical supports.
# Qubit index q = block*l*m + a*m + b for monomial x^a y^b (block 0 = L, 1 = R; Bravyi et al.'s convention).
# Run: julia --project=. scripts/make_bases.jl [OUTDIR]   (default handover/data_julia/bases)
const REPO = dirname(@__DIR__)
get!(ENV, "JULIA_PYTHONCALL_EXE", joinpath(REPO, ".venv", "bin", "python"))
get!(ENV, "JULIA_CONDAPKG_BACKEND", "Null")
using FTShadows

out = isempty(ARGS) ? joinpath(REPO, "handover", "data_julia", "bases") : ARGS[1]
abspath(out) == abspath(joinpath(REPO, "handover", "bases")) && error("refusing to write into handover/bases")
mkpath(out)
for (fname, rec) in FTShadows.handover_bases()
    FTShadows.write_meta(joinpath(out, fname), rec)
    if occursin("taskB", fname)
        println(split(rec["code"], ':')[1], " weights X ", length.(rec["xbar_support"]), " Z ", length.(rec["zbar_support"]))
    end
end
