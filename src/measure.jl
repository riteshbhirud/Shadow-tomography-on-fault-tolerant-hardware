# Glue between the circuit builder and the oracle: `oracle_circuit` for the builder's experiment types,
# and `measure`, the single entry point for plug-in settings.

for T in (Circuits.Experiment, Circuits.ComposedExperiment)
    @eval function oracle_circuit(exp::$T, p; eigen_targets::Bool = false, noise::AbstractString = "bravyi_cal")
        b = Circuits.build(exp, p; eigen_targets, noise)
        return OracleCircuit(exp.name, Circuits.stim_text(b.circuit), b.n_prep, b.target_records, b.detector_count, b.cost)
    end
end

"""
    measure(code, setting_id, state, p, shots; seed=0, with_gap=true, params...) -> OracleResult
    measure(code, setting_id, state; p, shots, seed=0, with_gap=true, params...)

Simulate and decode one native setting of a code plug-in on the logical state `[(label, sign bit)]`:
the plug-in builds the experiment and picks the decoder, then
[`run_oracle`](@ref) samples and decodes it. `params` are the operation's parameters, e.g. `d` and
`memory_rounds` (surface) or `rounds` (BB72 gauging).
"""
function measure(code::AbstractString, setting_id::AbstractString, state, p::Real, shots::Integer;
                 seed::Integer = 0, with_gap::Bool = true, params...)
    exp, labels, dec = Circuits.build_setting(code, setting_id, state; params...)
    return run_oracle(exp, p, shots; decoder = dec, seed, setting_labels = labels, with_gap)
end

measure(code::AbstractString, setting_id::AbstractString, state; p::Real, shots::Integer, kw...) =
    measure(code, setting_id, state, p, shots; kw...)
