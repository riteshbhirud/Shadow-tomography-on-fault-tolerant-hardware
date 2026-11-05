# Decoder-gap postselection on every noisy Nov-20 dataset: kept fraction, flip rates and fidelity estimate per gap cut.
# Run: julia --project=. scripts/nov20_postselection.jl
# Output: results/postselection_julia/
const REPO = dirname(@__DIR__)
get!(ENV, "JULIA_PYTHONCALL_EXE", joinpath(REPO, ".venv", "bin", "python"))
get!(ENV, "JULIA_CONDAPKG_BACKEND", "Null")
using FTShadows, PythonCall
using Printf: @sprintf

data = joinpath(REPO, "handover", "data", "nov20")
out = joinpath(REPO, "results", "postselection_julia")
mkpath(out)
prefixes = sort([joinpath(data, f[1:end-length(".meta.json")]) for f in readdir(data)
                 if endswith(f, ".meta.json") && FTShadows.read_meta(joinpath(data, f))["noise"]["p"] > 0])

mpl = pyimport("matplotlib")
mpl.use("Agg")
plt = pyimport("matplotlib.pyplot")
arr(x) = pyimport("numpy").array(x)  # matplotlib gets numpy arrays
fig, axes = plt.subplots(1, 2; figsize = (11, 4.2))
for prefix in prefixes
    name = basename(prefix)
    meta, rows, frows = FTShadows.postselection_analyse(prefix)
    isempty(rows) || FTShadows.pandas_csv(joinpath(out, "$name.csv"), rows)
    FTShadows.pandas_csv(joinpath(out, "$(name)_fidelity.csv"), frows)
    println(name)
    for r in frows
        println(join((@sprintf("%.4f", r[c]) for c in ("cut", "fraction_kept", "F", "se")), "  "))
    end
    flush(stdout)
    if occursin("random", name) && !occursin("logical", name)
        # kept shots and flips summed per cut, cuts ascending
        cuts = sort(unique(r["cut"] for r in rows))
        fk = [sum(r["shots_kept"] for r in rows if r["cut"] == c) for c in cuts]
        fl = [sum(r["flips"] for r in rows if r["cut"] == c) for c in cuts]
        tot = sum(r["shots_kept"] for r in rows if r["cut"] == minimum(cuts))
        lab = replace(name, "_random" => "")
        axes[0].plot(arr(fk ./ tot), arr(fl ./ fk), "o-"; ms = 4, lw = 1.5, label = lab)
        axes[1].errorbar(arr([r["fraction_kept"] for r in frows]), arr([r["F"] - r["F_expected_all_shots"] for r in frows]);
                         yerr = arr([r["se"] for r in frows]), fmt = "o-", ms = 4, lw = 1.5, capsize = 2, label = lab)
    end
end
axes[0].set_xlabel("fraction of shots kept")
axes[0].set_ylabel("flip rate of kept stabilizer outcomes")
axes[0].set_yscale("log")
axes[1].axhline(0; color = "0.5", lw = 1)
axes[1].set_xlabel("fraction of shots kept")
axes[1].set_ylabel("F estimate - noise-expected F (all shots)")
for ax in axes
    ax.invert_xaxis()
    ax.grid(; alpha = 0.3)
end
axes[1].legend(; fontsize = 6, ncol = 1, frameon = false)
fig.suptitle("Postselection on decoder gap, random-Pauli datasets (p = 1e-3)")
fig.tight_layout()
fig.savefig(joinpath(out, "postselection_summary.png"); dpi = 150)
