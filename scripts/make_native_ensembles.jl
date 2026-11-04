# handover/native_ensembles/<code>_lambda_w3.csv: exact lambda_P and 1/lambda_P for weight <= 3.
# Reads results/taskB/lambda_mix_*.npy.
# The gross 'aut_plus_all_folds' ensemble (unverified second fold family) is left out.
# Run: julia --project=. scripts/make_native_ensembles.jl [OUTDIR]   (default handover/data_julia/native_ensembles)
const REPO = dirname(@__DIR__)
using FTShadows

out = isempty(ARGS) ? joinpath(REPO, "handover", "data_julia", "native_ensembles") : ARGS[1]
abspath(out) == abspath(joinpath(REPO, "handover", "native_ensembles")) && error("refusing to write into handover/native_ensembles")
mkpath(out)
tb(f) = FTShadows.upstream(REPO, "results", "taskB", f)
for code in ("BB72", "gross")
    ens = [e for e in FTShadows.NATIVE_ENSEMBLES if !((code, e) in FTShadows.NATIVE_EXCLUDE)]
    lam = Dict(e => FTShadows.read_npy(tb("lambda_mix_$(code)_$e.npy")) for e in ens)
    t = FTShadows.native_ensemble_table(code, lam)
    FTShadows.write_csv_table(joinpath(out, "$(code)_lambda_w3.csv"), t)
    println(code, " ", length(t.rows), " rows")
end
