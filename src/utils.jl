"""
$(TYPEDSIGNATURES)
Generate a regular lattice with following parameters
- `dx` : element spacing on the x axis
- `dy` : element spacing on the y axis
- `ds` : displacement along x between rows of elements
- `f_cond`: function of x and y returning 1 if element at position x-y must be kept and 0 otherwise
"""
function generate_regular_lattice(dx::T, dy::T, ds::T, f_cond::Function = (x, y)->true;dx0 = T(0), dy0 = T(0), M::Int = 70,N::Int = M) where T<:Real

	# Function to generate x position as function of row,column number m,n
	x(m, n) = m * dx + n * ds + dx0
	# Function to generate y position as function of row,column number m,n
	y(n) = n * dy + dy0
	# Generate the elements. For each row, shift the columns to always have the search domain around x=0
	# out = [(x(m, n), y(n)) for n in -N:N for m in range(-M - round(n * ds / dx);length = 2 * M + 1) if f_cond(x(m, n), y(n))]
	out = [(x(m - round(Int,n * ds / dx), n), y(n)) for n in -N:N,m in -M:M if f_cond(x(m - round(Int,n * ds / dx), n), y(n))]
	return out
end
generate_regular_lattice(dx::Real,dy::Real,ds::Real,args...;kwargs...) = generate_regular_lattice(promote(dx,dy,ds)...,args...;kwargs...)

# Get the conversion from linear to db and viceversa
"""
$(TYPEDSIGNATURES)
Convert a number from linear to dB
"""
lin2db(x::Real) = 10log10(x)
"""
$(TYPEDSIGNATURES)
Convert a number from dB to linear
"""
db2lin(x::Real) = 10^(x/10)

# Convert between frequency and wavelength
"""
$(TYPEDSIGNATURES)
Get the wavelength (in m) starting from the frequency (in Hz)
"""
f2λ(f::Real) = c_0/f
"""
$(TYPEDSIGNATURES)
Get the frequency (in Hz) starting from the wavelength (in m) 
"""
λ2f(λ::Real) = c_0/λ