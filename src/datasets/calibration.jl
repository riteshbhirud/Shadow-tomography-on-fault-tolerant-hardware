# Flip rates behind the noise-expected values of the Nov-20 datasets: the surface-code single-patch
# calibration and its rows of calibration_local.csv. The BB72 rows take CPU-hours and are not regenerated here.

const SURFACE_CAL_P = 1e-3
const SURFACE_CAL_CHUNK = 100_000
const SURFACE_CAL_TARGET = 300
const SURFACE_CAL_MAX = Dict(3 => 2_000_000, 5 => 10_000_000)

"""
    surface_cal_job(job) -> (config, errors per target, product errors, cost)

One chunk of `SURFACE_CAL_CHUNK` shots of `job = (d, bases, r, chunk)` on +1 eigenstates, plain matching.
"""
function surface_cal_job(job)
    d, bases, r, chunk = job
    st = [(join(j == i ? b : 'I' for j in 1:length(bases)), 0) for (i, b) in enumerate(bases)]
    exp, _ = Circuits.surface_local_experiment(d, bases, st; memory_rounds = r)
    res = run_oracle(exp, SURFACE_CAL_P, SURFACE_CAL_CHUNK; seed = seed_for("surface_local_q", d, bases, r, chunk),
                     with_gap = false)
    o = Int.(res.outcomes)
    return (d, bases, r), vec(sum(o; dims = 1)), count(isodd, vec(sum(o; dims = 2))), res.cost
end

"""
    surface_calibration(cfgs; pmap=map) -> Dict(config => (errors, shots, cost))

Run chunks two at a time per configuration `(d, bases, r)` until every target has
`SURFACE_CAL_TARGET` errors or the shot cap for d is reached.
"""
function surface_calibration(cfgs; pmap = map)
    st = Dict(c => Dict{String,Any}("err" => nothing, "prod" => 0, "shots" => 0, "chunk" => 0, "cost" => nothing) for c in cfgs)
    active = Set(cfgs)
    while !isempty(active)
        jobs = Tuple[]
        for c in sort!(collect(active))
            for _ in 1:2
                push!(jobs, (c..., st[c]["chunk"]))
                st[c]["chunk"] += 1
            end
        end
        for (c, e, pr, cost) in pmap(surface_cal_job, jobs)
            x = st[c]
            x["err"] = x["err"] === nothing ? e : x["err"] .+ e
            x["prod"] += pr
            x["shots"] += SURFACE_CAL_CHUNK
            x["cost"] = cost
        end
        filter!(c -> !(minimum(st[c]["err"]) >= SURFACE_CAL_TARGET || st[c]["shots"] >= SURFACE_CAL_MAX[c[1]]), active)
    end
    return st
end

"""
    calibration_local_rows(st) -> Vector{OrderedDict}

calibration_local.csv rows (code surface, kind single_patch) of single-patch results, sorted by
(d, basis, memory rounds).
"""
function calibration_local_rows(st)
    rows = OrderedDict{String,Any}[]
    for (d, b, r) in sort!(collect(keys(st)))
        x = st[(d, b, r)]
        e = x["err"][1]
        lo, hi = wilson(e, x["shots"])
        push!(rows, OrderedDict{String,Any}("code" => "surface", "distance" => d, "setting_id" => "patch:$b", "r" => r,
            "target" => b, "kind" => "single_patch", "shots" => x["shots"], "flips" => e, "q" => e / x["shots"],
            "ci_lo" => lo, "ci_hi" => hi, "noise" => "bravyi_cal", "p" => SURFACE_CAL_P))
    end
    return rows
end
