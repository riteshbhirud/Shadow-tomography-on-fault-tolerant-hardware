# Shared plotting helpers: the series colours and a log of the data passed to each plot call (`--dump FILE` writes it as JSON).
using CairoMakie, JSON

const SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]  # fixed categorical order
const MARK = [:circle, :rect, :utriangle, :diamond, :dtriangle]
const GRID = "#e4e3df"
const REFGRAY = "#6b6a64"

const PLOTLOG = Any[]
_arr(x) = [isnan(v) ? nothing : Float64(v) for v in x]

"""Record one plot call: kind as in matplotlib (plot, semilogy, errorbar, hist, axhline), label, data."""
function rec!(kind; label=nothing, kw...)
    d = Dict{String,Any}("kind" => kind, "label" => label)
    for (k, v) in kw
        d[string(k)] = v isa AbstractVector ? _arr(v) : v
    end
    push!(PLOTLOG, d)
end

"""Save a figure and add a savefig marker to the plot log."""
function savefig(fig, path; px_per_unit=1)
    push!(PLOTLOG, Dict{String,Any}("kind" => "savefig", "name" => basename(path)))
    save(path, fig; px_per_unit)
end

dump_plotlog() = (i = findfirst(==("--dump"), ARGS); i === nothing || write(ARGS[i+1], JSON.json(PLOTLOG)))

# matplotlib-like axis: log y, light grid, no top/right spines
function styled_axis(pos; title="", xlabel="", ylabel="", ylog=true)
    ax = Axis(pos; title, titlealign=:left, titlesize=13, titlefont=:regular, xlabel, ylabel, xgridcolor=GRID, ygridcolor=GRID,
              topspinevisible=false, rightspinevisible=false, yscale=ylog ? log10 : identity)
    return ax
end
