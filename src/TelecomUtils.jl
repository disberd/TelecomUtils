module TelecomUtils

using DocStringExtensions
using PlutoDevMacros
using StaticArrays
using LinearAlgebra

# Create the filtering functions for include calls from PlutoDevMacros
# exclude_exprs = vcat(PlutoDevMacros.default_exprlist,[Symbol("@test"), :(using InteractiveUtils)])
# pluto_mapexpr = include_mapexpr(exclude_exprs)
# Constants
const c₀ = 299_792_458 # Speed of light [m/s]
const k_B = 1.38064852e-23 # Boltzmann constant [m² kg / (s² K)]

export wgs84_ellipsoid

include("../notebooks/utils.jl")
include("../notebooks/satview_basics.jl")
include("../notebooks/satview_transformations.jl")
include("../notebooks/satview_struct.jl")


end
