module TelecomUtils

using DocStringExtensions

# Constants
c_0 = 299_792_458 # Speed of light [m/s]
k_B = 1.38064852e-23 # Boltzmann constant [m² kg / (s² K)]
include("./utils.jl")

export c_0, k_B
export generate_regular_lattice, db2lin, lin2db, f2λ, λ2f
end
