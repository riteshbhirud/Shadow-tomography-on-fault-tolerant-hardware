# Python modules used through PythonCall. They are imported lazily, so
# loading FTShadows does not start Python until a bridged function is called.

const _PYMODULES = Dict{Symbol,Py}()

"""
    pymodule(name)

The Python module `name` (e.g. `:stim`, `:pymatching`, `:relay_bp`), imported
once per session from the interpreter PythonCall is configured to use.
"""
pymodule(name::Symbol) = get!(() -> pyimport(String(name)), _PYMODULES, name)
