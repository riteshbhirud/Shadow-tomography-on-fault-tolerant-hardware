# Where the scripts read their inputs. Scripts write to results/X_julia/ and handover/data_julia/,
# next to the stored outputs. A later script reads each input from the Julia folder if an earlier
# script has written it there, otherwise from the stored folder.

"""
    upstream(repo, parts...) -> path

The input file `joinpath(repo, parts...)`, taken from its Julia twin when
that file exists: `results/X/…` -> `results/X_julia/…`, `handover/data/…` and `handover/…` ->
`handover/data_julia/…`. Otherwise the stored path. Only single files are swapped: a Julia folder can
hold part of a dataset.
"""
function upstream(repo::AbstractString, parts::AbstractString...)
    twin = if length(parts) >= 2 && parts[1] == "results"
        joinpath(repo, "results", parts[2] * "_julia", parts[3:end]...)
    elseif length(parts) >= 2 && parts[1] == "handover"
        joinpath(repo, "handover", "data_julia", (parts[2] == "data" ? parts[3:end] : parts[2:end])...)
    end
    return twin !== nothing && isfile(twin) ? twin : joinpath(repo, parts...)
end
