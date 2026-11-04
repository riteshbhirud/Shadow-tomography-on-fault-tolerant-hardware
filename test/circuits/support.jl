# Shared by the circuit test items: rebuild a configuration of
# test/reference/circuits/configs.json (exported from the original implementation).
@testmodule CircRef begin

using FTShadows
using FTShadows.Circuits
using FTShadows.Circuits: PauliDict, Experiment, Epoch, BuildResult, build, stim_text, normalized_text, bb_setup,
                          bb_gauging, local_experiment, multi_patch_surgery, single_patch_readout, setting_experiment,
                          multi_patch_readout, bb72_readout_experiment, surface_local_experiment, memory_circuit,
                          code_capacity_circuit, multi_gauging_experiment, ordered_checks, exact_expectation, check_dicts,
                          logical_basis, logical_basis_with, dict_to_vec, vec_to_dict, pycircuit
using JSON, SHA, PythonCall

const CIRC_REF = joinpath(pkgdir(FTShadows), "test", "reference", "circuits")

circ_configs() = JSON.parsefile(joinpath(CIRC_REF, "configs.json"))
circ_codes() = JSON.parsefile(joinpath(CIRC_REF, "codes.json"))
circ_states() = JSON.parsefile(joinpath(CIRC_REF, "states.json"))

st(state) = [(String(l), Int(s)) for (l, s) in state]
pd(pairs) = PauliDict(Int(q) => only(s) for (q, s) in pairs)

function make_experiment(spec)
    k = spec["kind"]
    if k == "single_patch_readout"
        return single_patch_readout(spec["d"], only(spec["basis"]), st(spec["state"]), spec["rounds"])
    elseif k == "multi_patch_surgery"
        return multi_patch_surgery(spec["d"], spec["n_patches"], spec["type"], st(spec["state"]), spec["rounds"],
                                   spec["termination"])
    elseif k == "setting_experiment"
        return setting_experiment(spec["setting"], spec["d"], spec["n_patches"], st(spec["state"]), spec["rounds"])[1]
    elseif k == "multi_patch_readout"
        return multi_patch_readout(spec["d"], spec["bases"], st(spec["state"]), spec["rounds"])
    elseif k == "bb72_readout"
        return bb72_readout_experiment(only(spec["basis"]), st(spec["state"]))[1]
    elseif k == "bb_gauging"
        gkw = Dict{Symbol,Any}()
        if haskey(spec, "graph_kw")
            gkw[:edges] = Tuple(Tuple(e) for e in spec["graph_kw"]["edges"])
        end
        e = bb_gauging(only(spec["P"]), spec["rounds"]; state_kind=get(spec, "state_kind", "+1"),
                       termination=get(spec, "termination", "conjugate"), extra_targets=get(spec, "extra_targets", false),
                       order_seed=get(spec, "order_seed", nothing), separate_gauss=get(spec, "separate_gauss", false),
                       use_defaults=get(spec, "use_defaults", true), flux_keep=get(spec, "flux_keep", nothing),
                       state_labels=get(spec, "state_labels", nothing) === nothing ? nothing : st(spec["state_labels"]),
                       gkw...)
        get(spec, "tracker_readout_fallback", false) && (e.tracker_readout_fallback = true)
        return e
    elseif k == "bb_local"
        return local_experiment(spec["bases"], st(spec["state"]), spec["rounds"]; spectators=get(spec, "spectators", false))[1]
    elseif k == "surface_local"
        return surface_local_experiment(spec["d"], spec["bases"], st(spec["state"]);
                                        boundary_rounds=get(spec, "boundary_rounds", nothing),
                                        memory_rounds=spec["memory_rounds"])[1]
    elseif k == "multi_gauging"
        code = bb_setup().code
        Ls = [pd(L) for L in spec["Ls"]]
        pairs = logical_basis_with(code, [dict_to_vec(L, code.n) for L in Ls])
        basis = [(vec_to_dict(a), vec_to_dict(b)) for (a, b) in pairs]
        return multi_gauging_experiment(code.n, ordered_checks(code), Ls, spec["rounds"], only(spec["readout_basis"]);
                                        state=[(L, s) for (L, s) in zip(Ls, spec["signs"])], logical_basis=basis,
                                        extra_ops=[pd(op) for op in get(spec, "extra_ops", [])], name="bb72_multi_XX")
    elseif k == "memory_exp"
        s = bb_setup()
        code, basis = s.code, s.basis
        bk = only(spec["basis"])
        idx = bk == 'Z' ? 2 : 1
        targets = [pr[idx] for pr in basis if all(==(bk), values(pr[idx]))]
        return Experiment(; n_data=code.n, active0=0:code.n-1, prep_checks=check_dicts(code),
                          logical_state=[(t, 0) for t in targets], epochs=[Epoch(ordered_checks(code), spec["rounds"])],
                          readout=Dict(q => bk for q in 0:code.n-1), targets, logical_basis=basis)
    end
    error("unknown kind $k")
end

"""(circuit, build result, experiment) of a configuration (the last two `nothing` for bare circuits)."""
function build_config(spec)
    k = spec["kind"]
    if k == "bravyi_memory"
        return memory_circuit(spec["code"], only(spec["basis"]), spec["cycles"], spec["p"]), nothing, nothing
    elseif k == "code_capacity"
        e = make_experiment(spec["of"])
        ch = [c.pauli for c in e.epochs[1].checks]
        n = e.n_data
        S = Matrix{Bool}(permutedims(reduce(hcat, [dict_to_vec(c, n) for c in ch])))
        obs = [vec_to_dict(pr[spec["half"]+1]) for pr in logical_basis(S)]
        return code_capacity_circuit(n, ch, obs), nothing, nothing
    end
    e = make_experiment(spec)
    b = build(e, spec["p"]; eigen_targets=spec["eigen"], noise=spec["noise"])
    return b.circuit, b, e
end

sha_hex(s) = bytes2hex(sha256(s))

"""
    compare_config(rec; normalize=true)

`(ok, message)` for one exported configuration: the Julia text, normalized by stim, must hash to
the exported reference text, and the build record must agree.
"""
function compare_config(rec; normalize::Bool=true)
    spec = rec["spec"]
    if haskey(rec, "error")
        try
            build_config(spec)
        catch err
            return true, "raised as in the reference: $(sprint(showerror, err))"
        end
        return false, "the reference raised ($(rec["error"])) but Julia built a circuit"
    end
    c, b, e = build_config(spec)
    text = stim_text(c)
    raw_ok = sha_hex(text) == rec["sha256"]
    norm = normalize && !raw_ok ? normalized_text(text) : text
    sha_hex(norm) == rec["sha256"] || return false, "text differs (raw hash match: $raw_ok)"
    if b !== nothing
        b.n_meas == rec["n_meas"] || return false, "n_meas"
        b.n_prep == rec["n_prep"] || return false, "n_prep"
        b.detector_count == rec["detector_count"] || return false, "detector_count"
        b.target_records == rec["target_records"] || return false, "target_records"
        all(b.cost[k] == v for (k, v) in rec["cost"]) && length(b.cost) == length(rec["cost"]) || return false, "cost"
        length(b.epoch_info) == length(rec["epoch_info"]) || return false, "epoch_info"
        for (a, r) in zip(b.epoch_info, rec["epoch_info"])
            all(a[k] == v for (k, v) in r) || return false, "epoch_info"
        end
    end
    if haskey(rec, "graph_info")
        gi = e.graph_info
        (length(gi) == length(rec["graph_info"]) && all(gi[k] == v for (k, v) in rec["graph_info"])) ||
            return false, "graph_info"
    end
    return true, raw_ok ? "identical text" : "identical after normalization"
end

"""
    noiseless_outcomes(b, shots; seed=0)

Sample a p = 0 build with stim (bridge) and return the outcome bits (shots × targets), the parity
of each target's records; fails if any detector fires (oracle.run's noiseless path).
"""
function noiseless_outcomes(b, shots::Integer; seed::Integer=0)
    c = pycircuit(b.circuit)
    m = c.compile_sampler(; seed).sample(shots)
    dets = c.compile_m2d_converter().convert(; measurements=m, separate_observables=true)[0]
    any(pyconvert(Matrix{Bool}, dets)) && error("noiseless circuit fired a detector")
    M = pyconvert(Matrix{Bool}, m)
    return Matrix{Bool}(hcat([reduce(.⊻, (M[:, r+1] for r in recs); init=falses(shots)) for recs in b.target_records]...))
end

"""Number of subset products of the setting `bases` whose noiseless outcomes contradict the exact expectation."""
function check_products(state, bases, outcomes)
    x = 1 .- 2 .* Int.(outcomes)
    m = length(bases)
    bad = 0
    for mask in 1:2^m-1
        sub = [i for i in 1:m if (mask >> (i - 1)) & 1 == 1]
        P = join(i in sub ? bases[i] : 'I' for i in 1:m)
        e = exact_expectation(state, P)
        v = vec(prod(x[:, sub]; dims=2))
        bad += e in (1, -1) ? any(v .!= e) : length(unique(v)) != 2
    end
    return bad
end

"""Logical Z on qubit q of a state (signs of generators with X or Y there flip)."""
flipped(state, q) = [(l, s ⊻ Int(l[q+1] in "XY")) for (l, s) in state]

const BELL = [("XX", 0), ("ZZ", 0)]
const GHZ = [("XXX", 0), ("ZZI", 0), ("IZZ", 0)]

export noiseless_outcomes, check_products, flipped, BELL, GHZ, exact_expectation
export CIRC_REF, circ_configs, circ_codes, circ_states, st, pd, make_experiment, build_config, compare_config, sha_hex
export PauliDict, Experiment, Epoch, build, stim_text, normalized_text, bb_setup, bb_gauging, local_experiment,
       multi_patch_surgery, single_patch_readout, setting_experiment, multi_patch_readout, bb72_readout_experiment,
       surface_local_experiment, memory_circuit, code_capacity_circuit, ordered_checks, check_dicts, logical_basis,
       logical_basis_with, dict_to_vec, vec_to_dict, pycircuit, JSON

end
