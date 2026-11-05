# Cost of preparing a logical state vs measuring it with our local Pauli settings.
# Unit: physical qubit x round. Measurement volumes come from our circuits' cost records, preparation volumes from
# the surface GHZ surgery circuit or from the cited numbers in scripts/data/b2_cost_sources.json.

const PREP_VS_MEASURE_COLUMNS = ["task", "architecture", "n_logical", "prep_volume", "measure_volume",
    "ratio_measure_over_total", "prep_source", "measure_source", "assumptions"]

"""Space-time volume of a cost record: all qubits times (merged rounds + 1)."""
cost_volume(c) = (c["data_qubits"] + c["ancilla_qubits"]) * (c["merged_rounds"] + 1)

"""Python's `repr` of a small dict of str -> int, e.g. `{'X': 9, 'Y': 76, 'Z': 9}`."""
pyrepr_dict(d) = "{" * join(("'$k': $v" for (k, v) in d), ", ") * "}"

"""
    surface_patch_measure(d) -> (mean volume, per-basis volumes)

Expected volume of measuring one patch in a uniformly random basis (X, Y, Z), measured right after prep.
"""
function surface_patch_measure(d::Integer)
    v = OrderedDict(b => cost_volume(Circuits.build(Circuits.surface_local_experiment(d, string(b), [])[1], 1e-3).cost)
                    for b in "XYZ")
    return sum(values(v)) / 3, v
end

"""Mean volume of the BB72 settings on m qubits at r rounds (readout-only settings count at any r)."""
function bb72_measure(costs_csv::AbstractString, m::Integer, r::Integer)
    t = [x for x in csv_rows(costs_csv) if x["code"] == "bb72" && x["m"] == m && (x["r"] == r || x["n_gauged"] == 0)]
    v = [Float64((x["data_qubits"] + x["ancilla_qubits"]) * (x["merged_rounds"] + 1)) for x in t]
    return np_mean(v), length(t)
end

"""
    prep_vs_measure_rows(costs_csv, sources_json) -> Vector{OrderedDict}

The rows of `results/b2/prep_vs_measure.csv`. `costs_csv` is `local_setting_costs.csv`, `sources_json` the
cited numbers (`scripts/data/b2_cost_sources.json`). Missing values are `nothing`.
"""
function prep_vs_measure_rows(costs_csv::AbstractString, sources_json::AbstractString)
    rows = OrderedDict{String,Any}[]
    row(task, arch, n, prep, meas, ratio, ps, ms, as) = push!(rows, OrderedDict{String,Any}(zip(PREP_VS_MEASURE_COLUMNS,
        (task, arch, n, prep, meas, ratio, ps, ms, as))))
    for d in (3, 5, 7)
        pm, per = surface_patch_measure(d)
        for n in (3, 6, 12)
            prep = cost_volume(Circuits.build(Circuits.multi_patch_surgery(d, n, "ZZ", [], d), 1e-3).cost)
            meas = n * pm
            row("GHZ_n (prep: |+>^n, then all Z_iZ_{i+1} by one lattice-surgery step of d rounds)", "surface d=$d", n,
                prep, meas, meas / (prep + meas),
                "our circuit: surface.multi_patch_surgery(d, n, 'ZZ', rounds=d), cost record",
                "our circuits: per-patch X/Y/Z volumes $(pyrepr_dict(per))",
                "patch initialisation and Pauli-frame fixing not counted (prep is a lower " *
                "bound); surgery schedule not optimised (schedule factor)")
        end
    end
    for r in (3, 6)
        meas, nset = bb72_measure(costs_csv, 3, r)
        row("GHZ_3 on logical qubits 0-2 of one block", "BB72 (r = $r)", 3, nothing, meas, nothing,
            "not modelled: no published in-block GHZ / CNOT cost for BB72 (B2 sources: 'not " *
            "found'); Tour de gross gives gross-code costs in gate-layer timesteps only",
            "our compiled settings, mean over $nset settings (local_setting_costs.csv)",
            "unused logical qubits idle in the same block")
    end
    src = OrderedDict(e["quantity"] => e for e in JSON.parsefile(sources_json))
    cult3 = src["expected qubit-rounds per kept T, d1=3"]
    cult5 = src["expected qubit-rounds per kept T, d1=5"]
    lit = first(e for (q, e) in src if startswith(q, "(15-to-1)_{17,7,7}"))
    pm11, per11 = surface_patch_measure(11)
    for (name, prep, e) in (("cultivation d1=3 (p_out 3e-6)", cult3["value"], cult3),
                            ("cultivation d1=5 (p_out 2e-9)", cult5["value"], cult5),
                            ("15-to-1 (17,7,7) distillation (p_out 4.5e-8)", lit["value"]["qubitcycles"], lit))
        row("single |T> state: $name", "surface d=11 (output patch)", 1, prep, pm11, pm11 / (prep + pm11),
            "$(e["source"]), $(e["location"])", "our circuits: $(pyrepr_dict(per11))",
            "p = 1e-3; cultivation volume includes retries and is 'slightly optimistic' " *
            "(ignores packing); cultivation ends in a grafted code ~ d = 11 surface code " *
            "(2409.17595 Fig. 16); Y measurement at d = 11 uses b = d//2 (cost only, distance " *
            "not verified at d = 11)")
    end
    fl = first(e for (q, e) in src if startswith(q, "11x11 TFIM"))["value"]
    d = fl["d"]
    prep = fl["physical_qubits"] * fl["L_FLASQ"] * d  # logical timestep = d cycles (FLASQ Sec. II B)
    pm13, per13 = surface_patch_measure(d)
    meas = 121 * pm13  # the 121 spins of the 11 x 11 lattice
    for (steps, label) in ((20, "20 Trotter steps"), (1, "1 Trotter step (1/20 of the run, our division)"))
        p_ = prep * steps / 20
        row("11x11 transverse-field Ising, 2nd-order Trotter: $label", "surface d=$d", 121, p_, meas, meas / (p_ + meas),
            "arXiv:2511.08508 (FLASQ) Table II / Table VII: 62,720 physical qubits, 55,995 " *
            "logical timesteps at d = 13, p = 1e-3",
            "our circuits, 121 patches at d = 13: $(pyrepr_dict(per13))",
            "all 62,720 qubits held for the whole run (FLASQ's own spacetime accounting); " *
            "only the 121 system qubits are measured; Y at d = 13 uses b = d//2 (cost only)")
    end
    return rows
end

"""
    write_prep_vs_measure(path, rows)

Write the rows as CSV in pandas format. The volume and ratio columns have missing values, so they are
float columns (`260.0`, missing empty).
"""
function write_prep_vs_measure(path::AbstractString, rows)
    floatcols = ("prep_volume", "measure_volume", "ratio_measure_over_total")
    cell(c, v) = v === nothing ? missing : (c in floatcols ? Float64(v) : v)
    write_csv_table(path, CSVTable(PREP_VS_MEASURE_COLUMNS,
                                   [Any[cell(c, r[c]) for c in PREP_VS_MEASURE_COLUMNS] for r in rows], "\n"))
end
