# Dataset tables (setting costs, calibration tables, the v1 and Nov-20 datasets, the basis and native-ensemble
# exports) and the calibration, verification and postselection studies behind results/. The scripts in
# scripts/ run them and write to new folders.

using OrderedCollections: OrderedDict
using Printf: @sprintf

include("numerics.jl")
include("setting_costs.jl")
include("calibration.jl")
include("local_datasets.jl")
include("handover_tables.jl")
include("v1_datasets.jl")
include("prep_vs_measure.jl")
include("pandas_csv.jl")
include("surface_studies.jl")
include("bravyi_calibration.jl")
include("dataset_settings.jl")
include("postselection.jl")
