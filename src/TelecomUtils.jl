module TelecomUtils

using DocStringExtensions
using PlutoDevMacros
using StaticArrays
using LinearAlgebra

# Create the filtering functions for include calls from PlutoDevMacros
exclude_exprs = vcat(PlutoDevMacros.default_exprlist,[Symbol("@test"), :(using InteractiveUtils)])
pluto_mapexpr = include_mapexpr(exclude_exprs)
# Constants
const c₀ = 299_792_458 # Speed of light [m/s]
const k_B = 1.38064852e-23 # Boltzmann constant [m² kg / (s² K)]
# Initialize the vector that contains the matrix to compute the beam coloring. We limit ourselves at 500 colors to start
const F_reuse_matrix = (square = SMatrix{2,2,Float64,4}[], triangular = SMatrix{2,2,Float64,4}[])
include("./utils.jl")
include(pluto_mapexpr,"./SatView.jl")

export generate_regular_lattice, generate_hex_lattice, generate_square_lattice, generate_rect_lattice, db2lin, lin2db, f2λ, λ2f
export generate_colors
end
