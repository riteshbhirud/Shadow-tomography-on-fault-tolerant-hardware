# Compiled cost of every local Pauli setting on 2 and 3 logical qubits:
# surface code d = 3, 5, 7 (memory rounds r, then per-patch readout, Y by Gidney's measurement) and BB72
# (readout plus simultaneous gauging for r rounds). Costs come from the built circuits; flip rates are
# joined from the calibration tables.

"""Rows of a CSV file as dicts from column name to typed cell (see `CSVTable`)."""
csv_rows(path::AbstractString) = (t = read_csv_table(path); [Dict(zip(t.header, r)) for r in t.rows])

# Python's builtin max: the first element unless a later one compares greater (so NaN handling follows `>`).
pymax(xs) = foldl((a, b) -> b > a ? b : a, xs)

const COST_COLUMNS = ["code", "distance", "m", "setting", "r", "implementation", "readout", "gauged", "n_gauged",
    "needs_surgery", "merged_rounds", "data_qubits", "initial_data_qubits", "ancilla_qubits",
    "initial_ancilla_qubits", "extra_qubits", "circuit_layers", "q_max_string"]

"""
    local_setting_costs(surface_q, bb72_q) -> Vector{OrderedDict}

The rows of `results/taskA/local_setting_costs.csv`, from the calibration tables `surface_local_q.csv` and
`bb72_local_q.csv` (paths). The BB72 flip rates are joined in: the largest string q at the setting's rounds (3 for readout-only settings).
"""
function local_setting_costs(surface_q::AbstractString, bb72_q::AbstractString)
    sq = csv_rows(surface_q)
    bq = [r for r in csv_rows(bb72_q) if r["kind"] == "string"]
    rows = OrderedDict{String,Any}[]
    seen = Set{Tuple}()
    function add!(row)
        key = (row["code"], row["distance"], row["m"], row["setting"], row["r"])
        key in seen || (push!(seen, key); push!(rows, row))
    end
    for m in (2, 3), bases in product_strings("XYZ", m)
        for d in (3, 5, 7), r in (3, 6)
            exp, _ = Circuits.surface_local_experiment(d, bases, []; memory_rounds = r)
            c = Circuits.build(exp, 1e-3).cost
            qs = Dict(only(x["basis"]) => Float64(x["q"]) for x in sq if x["d"] == d && x["memory_rounds"] == r)
            nY = count(==('Y'), bases)
            impl = "per-patch destructive readout" *
                   (nY > 0 ? "; Y on $nY patch(es) by in-place Y measurement (b = $(Circuits.Y_DEFAULT_ROUNDS[d]))" : "")
            add!(OrderedDict{String,Any}("code" => "surface", "distance" => d, "m" => m, "setting" => bases, "r" => r,
                "implementation" => impl, "readout" => "per patch", "gauged" => "", "n_gauged" => 0,
                "needs_surgery" => false, c..., "q_max_string" => pymax(get(qs, b, NaN) for b in bases)))
        end
        for r in (3, 6)
            exp, _, comp = Circuits.local_experiment(bases, [], r)
            c = Circuits.build(exp, 1e-3).cost
            g = comp.gauged
            R = string(comp.readout)
            impl = isempty(g) ? "destructive readout in $R" :
                   "simultaneous gauging of " * join(("$(P)bar_$i" for (i, P) in g), ", ") * " for $r rounds, then readout in $R"
            rr = isempty(g) ? 0 : r
            qrows = [x["q"] for x in bq if x["setting_id"] == "local:$bases" && x["rounds_parameter"] == (rr == 0 ? 3 : rr)]
            add!(OrderedDict{String,Any}("code" => "bb72", "distance" => 6, "m" => m, "setting" => bases, "r" => rr,
                "implementation" => impl, "readout" => R, "gauged" => join(("$P$i" for (i, P) in g), ";"),
                "n_gauged" => length(g), "needs_surgery" => !isempty(g), c...,
                "q_max_string" => isempty(qrows) ? NaN : Float64(maximum(qrows))))
        end
    end
    return rows
end

_cell(x::Bool) = x ? "True" : "False"
_cell(x::AbstractFloat) = isnan(x) ? missing : Float64(x)
_cell(x) = x

"""
    write_setting_costs(path, rows)

Write the cost rows as pandas' `to_csv(index=False)` does.
"""
function write_setting_costs(path::AbstractString, rows)
    write_csv_table(path, CSVTable(COST_COLUMNS, [Any[_cell(r[c]) for c in COST_COLUMNS] for r in rows], "\n"))
end
