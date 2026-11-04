# BB memory experiments for the syndrome-schedule search and the memory sanity check.

"""
    css_basis(code)

A CSS logical basis `[(Xbar_i, Zbar_i)]` (Pauli dicts) with X-type Xbar and Z-type Zbar: the null
spaces of H_Z and H_X reduced modulo the stabilizers, then paired by Z' = M^{-T} Z.
"""
function css_basis(code::BBCode)
    n = code.n
    NZ, NX = nullspace_gf2(code.HZ), nullspace_gf2(code.HX)
    function indep(vs)
        cur = RowSpace(code.S)
        return [Vector{Bool}(v) for v in vs if addrow!(cur, v)]
    end
    Xs = indep(vcat(NZ[i, :], falses(n)) for i in axes(NZ, 1))
    Zs = indep(vcat(falses(n), NX[i, :]) for i in axes(NX, 1))
    @assert length(Xs) == length(Zs) == code.k
    return LogicalBasis([(vec_to_dict(a), vec_to_dict(b)) for (a, b) in _pair_rest([], Xs, Zs)])
end

"""
    memory_experiment(code, checks, basis_kind, rounds, logical_basis)

Memory in `basis_kind` ('X' or 'Z'): ideal preparation, `rounds` noisy rounds of `checks`, readout
of every qubit in that basis; every logical of the basis that is purely of that type is a +1 target.
"""
function memory_experiment(code::BBCode, checks, basis_kind::Char, rounds::Integer, logical_basis)
    idx = basis_kind == 'Z' ? 2 : 1
    targets = [pr[idx] for pr in logical_basis if all(==(basis_kind), values(pr[idx]))]
    return Experiment(; n_data=code.n, active0=0:code.n-1, prep_checks=check_dicts(code),
                      logical_state=[(t, 0) for t in targets], epochs=[Epoch(checks, rounds)],
                      readout=OrderedDict(q => basis_kind for q in 0:code.n-1), targets, logical_basis)
end
