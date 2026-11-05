# Build the documentation into docs/build/. From the repository root:
#   export JULIA_PYTHONCALL_EXE="$PWD/.venv/bin/python" JULIA_CONDAPKG_BACKEND=Null
#   julia --project=docs -e 'using Pkg; Pkg.instantiate()'
#   julia --project=docs docs/make.jl
# The examples call stim and the decoders through PythonCall, so the venv of python/requirements.txt is needed.
const REPO = dirname(@__DIR__)
get!(ENV, "JULIA_PYTHONCALL_EXE", joinpath(REPO, ".venv", "bin", "python"))
get!(ENV, "JULIA_CONDAPKG_BACKEND", "Null")

using Documenter
using DocumenterCitations
using FTShadows

ENV["LINES"] = 80     # for forcing `displaysize(io)` to be big enough
ENV["COLUMNS"] = 100

bib = CitationBibliography(joinpath(@__DIR__, "src", "references.bib"); style = :authoryear)

makedocs(;
    plugins = [bib],
    doctest = false,
    clean = true,
    sitename = "FTShadows.jl",
    authors = "The FTShadows authors",
    modules = [FTShadows, FTShadows.Circuits, FTShadows.TaskB],
    repo = Remotes.GitHub("riteshbhirud", "Shadow-tomography-on-fault-tolerant-hardware"),
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://riteshbhirud.github.io/Shadow-tomography-on-fault-tolerant-hardware/",
        size_threshold_ignore = ["API.md", "manual.md"],
    ),
    checkdocs = :exports,
    warnonly = [:missing_docs],
    pages = [
        "FTShadows.jl" => "index.md",
        "Manual" => [
            "Measuring logical settings" => "manual.md",
            "Decoders and the gap" => "gap.md",
            "Noise models" => "noise.md",
            "Shot records and datasets" => "datasets.md",
            "Conventions" => "conventions.md",
        ],
        "Native Clifford ensembles (Task B)" => "taskb.md",
        "Python dependencies" => "python.md",
        "API" => "API.md",
        "References" => "references.md",
    ],
)
