# Task B shadow-map figures (mixed readout): learnable fraction and median squared norm per weight, norm histogram of the largest ensemble.
# Run: julia --project=scripts/plots scripts/plots/taskb_shadow_plots.jl DIR [BB72] [gross] [--dump FILE]
# Reads shadow_<code>.json and lambda_mix_*.npy from DIR and writes the PNGs there.
using Printf
import NPZ
include(joinpath(@__DIR__, "common.jl"))

# matplotlib's default colour cycle
const TAB10 = Makie.to_colormap(:tab10)

# logical weight of every Pauli index v = x | z << k
logical_weights(k) = [count_ones((v & (1 << k - 1)) | (v >> k)) for v in 0:(1<<2k)-1]

function shadow_plots(dir, name)
    res = JSON.parsefile(joinpath(dir, "shadow_$name.json"))
    k = res["k"]
    ens_names = collect(keys(res["ensembles"]))
    fig = Figure(size=(1200, 450))
    ax1 = styled_axis(fig[1, 1]; xlabel="logical weight w", ylabel="fraction of weight-w logical Paulis learnable")
    ax2 = styled_axis(fig[1, 2]; xlabel="logical weight w", ylabel="median squared shadow norm (learnable ones)")
    ax1.xticks = ax2.xticks = 2:2:k
    ws = collect(1:k)
    for (i, e) in enumerate(ens_names)
        rows = res["ensembles"][e]["readouts"]["mix"]["learnable_by_weight_all"]
        col = TAB10[mod1(i, 10)]
        fr = Float64[r["frac_learnable"] for r in rows]
        lab = "$e (|G|=$(@sprintf("%.3g", Float64(res["ensembles"][e]["order"]))))"
        rec!("plot"; label=lab, x=ws, y=fr)
        # a zero fraction cannot be shown on the log axis (matplotlib drops it too)
        scatterlines!(ax1, ws, [f > 0 ? f : NaN for f in fr]; color=col, label=lab)
        med = Float64[r["norm2_median"] === nothing ? NaN : r["norm2_median"] for r in rows]
        rec!("semilogy"; label=e, x=ws, y=med)
        scatterlines!(ax2, ws, med; color=col, label=e)
    end
    rec!("semilogy"; label="local random Pauli 3^w", x=ws, y=3.0 .^ ws)
    lines!(ax2, ws, 3.0 .^ ws; color=:black, linestyle=:dash, label="local random Pauli 3^w")
    rec!("axhline"; label="global Clifford 2^k+1", y=[2.0^k + 1])
    hlines!(ax2, 2^k + 1; color=:gray, linestyle=:dot, label="global Clifford 2^k+1")
    axislegend(ax1; labelsize=9, position=:lb)
    axislegend(ax2; labelsize=9, position=:lt)
    Label(fig[0, :], "$name: native ensemble 'random g in G, Pauli frame, all-Z or all-X readout (50/50)'"; fontsize=15)
    savefig(fig, joinpath(dir, "shadow_$(name)_by_weight.png"); px_per_unit=1.3)

    big = ens_names[end]
    lam = NPZ.npzread(joinpath(dir, "lambda_mix_$(name)_$big.npy"))
    wts = logical_weights(k)
    fig = Figure(size=(700, 450))
    ax = styled_axis(fig[1, 1]; title="$name, ensemble '$big', mixed readout", xlabel="log10 squared shadow norm",
                     ylabel="count", ylog=false)
    for w in 1:4
        lw = lam[wts.==w]
        lw = lw[lw.>0]
        isempty(lw) && continue
        x = log10.(1.0f0 ./ lw)   # float32, as numpy computes it
        rec!("hist"; label="w=$w (learnable $(length(lw)))", x=x, bins=60)
        stephist!(ax, x; bins=60, label="w=$w (learnable $(length(lw)))", color=TAB10[w])
    end
    axislegend(ax; position=:lt)
    savefig(fig, joinpath(dir, "shadow_$(name)_norm_hist_$big.png"); px_per_unit=1.3)
end

let args = filter(a -> !startswith(a, "--"), ARGS)
    i = findfirst(==("--dump"), ARGS)
    i === nothing || filter!(!=(ARGS[i+1]), args)
    dir = args[1]
    for name in (length(args) > 1 ? args[2:end] : ["BB72", "gross"])
        shadow_plots(dir, name)
    end
    dump_plotlog()
end
