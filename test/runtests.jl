# The bridged Python tools come from the project venv built from python/requirements.txt unless the
# caller configured PythonCall.
const REPO = dirname(@__DIR__)
get!(ENV, "JULIA_PYTHONCALL_EXE", joinpath(REPO, ".venv", "bin", "python"))
get!(ENV, "JULIA_CONDAPKG_BACKEND", "Null")

using TestItemRunner

# FTSHADOWS_TEST=fast skips items tagged :slow (the CI subset).
# FTSHADOWS_TEST_FILES=decode,oracle runs only the files whose names contain one of the words.
const FAST = get(ENV, "FTSHADOWS_TEST", "") == "fast"
const FILES = filter(!isempty, split(get(ENV, "FTSHADOWS_TEST_FILES", ""), ","))
keep(ti) = !(FAST && :slow in ti.tags) && (isempty(FILES) || any(f -> occursin(f, basename(ti.filename)), FILES))
@run_package_tests filter = keep
