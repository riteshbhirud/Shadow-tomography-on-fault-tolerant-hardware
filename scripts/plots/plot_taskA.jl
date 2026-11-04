# Task A plots: surface q(r), BB72 q(r) and BB72 postselection curves.
# Run: julia --project=scripts/plots scripts/plots/plot_taskA.jl [--in DIR] [--out DIR] [--dump FILE]
# Reads results/taskA/*.csv and *.npz; writes the PNGs to results/taskA_julia/, not next to the data.
import CSV, NPZ
include(joinpath(@__DIR__, "common.jl"))

const REPO = dirname(dirname(@__DIR__))
argval(flag, default) = (i = findfirst(==(flag), ARGS); i === nothing ? default : ARGS[i+1])
const IN = argval("--in", joinpath(REPO, "results", "taskA"))
const OUT = argval("--out", joinpath(REPO, "results", "taskA_julia"))

const REF_SURF = Dict(5 => Dict(1 => 8.3e-2, 2 => 3.4e-2, 3 => 3.3e-3, 4 => 1.2e-3, 5 => 4.3e-4),
                      7 => Dict(1 => 1.2e-1, 3 => 4.5e-3, 5 => 3.3e-4))
const REF_SURF_PLATEAU = Dict(5 => 2.9e-4, 7 => 4.0e-5)
const REF_BB = Dict("X" => [1.13e-1, 2.84e-2, 3.20e-3, 5.3e-4], "Z" => [9.7e-2, 2.27e-2, 1.97e-3, 6.7e-4],
                    "Y" => [1.87e-1, 5.37e-2, 5.91e-3, 1.90e-3])

# rows of a CSV as named tuples, every column read as text and parsed where needed (the state column is "+1")
readrows(path) = collect(CSV.File(path; types=String))
num(x) = parse(Float64, x)

# q(r) with 95% Wilson error bars, colour and marker i
function errbar!(ax, x, y, lo, hi, i, label)
    rec!("errorbar"; label, x, y, yerr_lo=y .- lo, yerr_hi=hi .- y)
    errorbars!(ax, x, y, y .- lo, hi .- y; color=SERIES[i], whiskerwidth=4)
    scatterlines!(ax, x, y; color=SERIES[i], marker=MARK[i], markersize=10, linewidth=2, label)
end

function surface()
    rows = readrows(joinpath(IN, "surface_qr.csv"))
    fig = Figure(size=(1100, 480))
    axes = [styled_axis(fig[1, j]; title="Surface code d=$d, p=1e-3 (bravyi_cal): q(r)", xlabel="merged rounds r")
            for (j, d) in enumerate((5, 7))]
    linkyaxes!(axes...)
    foreach(ax -> ax.xticks = 2:2:12, axes)   # whole rounds only
    sel = [("ZZ", "conjugate", "bravyi_cal", "ZZ, conjugate"), ("XX", "conjugate", "bravyi_cal", "XX, conjugate"),
           ("ZZ", "same", "bravyi_cal", "ZZ, same-basis (optimistic)"),
           ("ZZ", "conjugate", "stim_generated", "ZZ, conjugate (stim-generated noise)"),
           ("ZZ", "conjugate", "sd6", "ZZ, conjugate (sd6, pessimistic)")]
    for (ax, d) in zip(axes, (5, 7))
        for (i, (k, t, nz, lab)) in enumerate(sel)
            g = [r for r in rows if parse(Int, r.d) == d && r.kind == k && r.term == t && r.noise == nz]
            sort!(g; by=r -> parse(Int, r.r))
            isempty(g) && continue
            errbar!(ax, num.(getproperty.(g, :r)), num.(getproperty.(g, :q)), num.(getproperty.(g, :ci_lo)),
                    num.(getproperty.(g, :ci_hi)), i, lab)
        end
        rr = sort(collect(keys(REF_SURF[d])))
        rec!("plot"; label="earlier reference runs", x=Float64.(rr), y=[REF_SURF[d][r] for r in rr])
        scatter!(ax, rr, [REF_SURF[d][r] for r in rr]; marker=:circle, color=:transparent, strokecolor=REFGRAY,
                 strokewidth=1.5, markersize=13, label="earlier reference runs")
        rec!("axhline"; label="reference plateau", y=[REF_SURF_PLATEAU[d]])
        hlines!(ax, REF_SURF_PLATEAU[d]; color=REFGRAY, linestyle=:dot, linewidth=1, label="reference plateau")
    end
    axes[1].ylabel = "q(r)  (95% Wilson CI)"
    Legend(fig[2, :], axes[1]; orientation=:horizontal, nbanks=3, framevisible=false, labelsize=11)
    savefig(fig, joinpath(OUT, "surface_qr.png"); px_per_unit=1.5)
end

function bb72()
    f = joinpath(IN, "bb72_qr.csv")
    isfile(f) || return
    rows = readrows(f)
    fig = Figure(size=(650, 420))
    ax = styled_axis(fig[1, 1]; title="BB72 gauging (common basis), p=1e-3 bravyi_cal, Relay-BP: q(r)",
                     xlabel="merged rounds r", ylabel="q(r)  (95% Wilson CI)", )
    ax.xticks = 1:6
    for (i, P) in enumerate(("X", "Z", "Y"))
        g = sort([r for r in rows if r.P == P && r.state in ("+1", "1")]; by=r -> parse(Int, r.r))
        ex = parse(Int, g[1].cost_extra_qubits)
        errbar!(ax, num.(getproperty.(g, :r)), num.(getproperty.(g, :q)), num.(getproperty.(g, :ci_lo)),
                num.(getproperty.(g, :ci_hi)), i, "$P  (+$ex qubits)")
        rec!("plot"; x=collect(1.0:4.0), y=REF_BB[P])
        scatter!(ax, 1:4, REF_BB[P]; marker=MARK[i], color=:transparent, strokecolor=(SERIES[i], 0.8),
                 strokewidth=1.2, markersize=14)
    end
    rec!("plot"; label="hollow = earlier reference runs", x=Float64[], y=Float64[])
    scatter!(ax, [NaN], [NaN]; marker=:circle, color=:transparent, strokecolor=REFGRAY, strokewidth=1.2,
             label="hollow = earlier reference runs")
    axislegend(ax; framevisible=false, labelsize=11)
    savefig(fig, joinpath(OUT, "bb72_qr.png"); px_per_unit=1.5)

    fig = Figure(size=(650, 420))
    ax = styled_axis(fig[1, 1]; xlabel="fraction of shots kept (highest gap first)", ylabel="q among kept shots",
                     title="BB72: confidence postselection, plain-anchored Relay gap (3e4 shots each;\n" *
                           "curves end where no kept shot is an error; infinite-gap shots are tied)")
    ax.titlesize = 12
    for (i, P) in enumerate(("X", "Z", "Y")), r in (2, 3)   # colour = P (as in the q(r) plot), dash = r
        fz = joinpath(IN, "bb72_post_$(P)_r$r.npz")
        isfile(fz) || continue
        z = NPZ.npzread(fz)
        # np.nan_to_num(gap, nan=0.0): NaN -> 0, ±inf -> the largest finite values; then a stable sort, highest first
        gap = [isnan(x) ? 0.0 : isinf(x) ? copysign(floatmax(Float64), x) : x for x in z["gap"]]
        o = sortperm(-gap; alg=MergeSort)
        err = Int.(z["err"][o])
        n = collect(1:length(err))
        keep = n ./ length(err)
        q = cumsum(err) ./ n
        m = (keep .>= 0.05) .& (cumsum(err) .> 0)   # no point where zero errors were kept
        rec!("plot"; label="$P, r=$r", x=keep[m], y=q[m])
        lines!(ax, keep[m], q[m]; color=SERIES[i], linewidth=2, linestyle=r == 3 ? :solid : :dash, label="$P, r=$r")
    end
    axislegend(ax; framevisible=false, labelsize=11, nbanks=2, position=:lt)
    savefig(fig, joinpath(OUT, "bb72_postselection.png"); px_per_unit=1.5)
end

mkpath(OUT)
surface()
bb72()
dump_plotlog()
println("plots written")
