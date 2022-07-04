### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ c0a30957-4c7b-4d7b-bfa9-c2fb691a077b
# ╠═╡ skip_as_script = true
#=╠═╡
begin
	using Revise
	using PlutoUtils
	using PlotlyBase
	using BenchmarkTools
end
  ╠═╡ =#

# ╔═╡ 74975885-9a4e-4857-8135-9e4f69061caf
begin
	using DocStringExtensions
	using StaticArrays
	using LinearAlgebra
	using Dictionaries
	using SplitApplyCombine
end

# ╔═╡ 379613ec-0973-4000-ae8c-d7c33ddca18e
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
# Packages
"""
  ╠═╡ =#

# ╔═╡ 736b0cf6-bec2-4226-8ef4-70f6a865d34a
# ╠═╡ skip_as_script = true
#=╠═╡
ToC()
  ╠═╡ =#

# ╔═╡ f8243a65-9f5e-464e-bb06-0bb4f5131b8b
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
# Exports
"""
  ╠═╡ =#

# ╔═╡ 7f645e69-3334-44db-9ba1-9f2d3e0127a2
const ColorOrderDict{I} = Dictionary{SVector{2, I}, I}

# ╔═╡ 8660a7c4-eb78-4e7c-966b-d759df7f3dfa
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
# Lattice Functions
"""
  ╠═╡ =#

# ╔═╡ 71163795-9695-4f11-acc2-6e3838c8a158
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## generate\_regular_lattice
"""
  ╠═╡ =#

# ╔═╡ da97848f-a7ff-4f2d-b98d-e8bf1ccc3038
function lattice_generator(dx::T, dy::T, ds::T;x0::T = T(0), y0::T = T(0), M::Int = 70,N::Int = M) where T<:Real
	# Function to generate x position as function of row,column number m,n
	x(m, n) = m * dx + n * ds + x0
	# Function to generate y position as function of row,column number m,n
	y(n) = n * dy + y0
	# Generate the elements. For each row, shift the columns to always have the search domain around x=0
	gen = (SVector(x(m - round(Int,n * ds / dx), n), y(n)) for n in -N:N,m in -M:M)
	return gen
end

# ╔═╡ f8a53711-e07f-4b6b-84ea-803679496571
"""
Generate a regular lattice of points
$(TYPEDSIGNATURES)

# Arguments
- `dx` → element spacing on the x axis
- `dy` → element spacing on the y axis
- `ds` → displacement along x between rows of elements
- `f_cond::Function` → function of two arguments (`f(x,y) = ...`) returning `true` if element at position `x,y` must be kept and `false` otherwise

# Keyord Arguments
- `x0 = 0` → x coordinate of the origin of the lattice
- `y0 = 0` → y coordinate of the origin of the lattice
- `M::Int = 70` → Number of elements to generate per row of points before appliying the filtering function `f_cond`
- `N::Int = M` → Number of rows of points to generate before appliying the filtering function `f_cond`
"""
function generate_regular_lattice(dx::T, dy::T, ds::T, f_cond::Function = (x, y) -> true;x0::T = T(0), y0::T = T(0), M::Int = 70,N::Int = M) where T<:Real
	gen = lattice_generator(dx, dy, ds; x0, y0, M, N)
	return [x for x ∈ gen if f_cond(x...)]
end

# ╔═╡ ce22b91e-6bba-4312-a89b-1a78f84034d3
function regular_lattice_nelements(dx::T, dy::T, ds::T, f_cond::Function = (x, y) -> true;x0::T = T(0), y0::T = T(0), M::Int = 70,N::Int = M) where T<:Real
	gen = lattice_generator(dx, dy, ds; x0, y0, M, N)
	sum(x -> f_cond(x...), gen)
end

# ╔═╡ 0eb5d19d-b535-4cda-a89b-26ba295e2711
generate_regular_lattice(dx::Real,dy::Real,ds::Real,args...;kwargs...) = generate_regular_lattice(promote(dx,dy,ds)...,args...;kwargs...)

# ╔═╡ f0834b38-8efe-4e77-b0f9-47e5b7595191
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## generate\_rect_lattice
"""
  ╠═╡ =#

# ╔═╡ 9e5db472-f96a-4acb-96ae-024b5c73a93d
"""
    generate_rect_lattice(spacing_x::Real,spacing_y::Real[,f_cond];kwargs...)
# Summary
Generate a rectangular lattice of points (with different spacing among x and y directions)
# Arguments
- `spacing_x` → spacing between points on the x axis
- `spacing_y` → spacing between points on the y axis

See [`generate_regular_lattice`](@ref) for a description of `f_cond` and of  the keyword arguments
"""
generate_rect_lattice(spacing_x::Real,spacing_y::Real,args...;kwargs...) = generate_regular_lattice(spacing_x,spacing_y,0,args...;kwargs...)

# ╔═╡ c41f9f41-d8bd-4001-98cb-2ab788404b1b
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## generate\_square\_lattice
"""
  ╠═╡ =#

# ╔═╡ 8787134f-9d14-4329-8dda-72557e3175b8
"""
`generate_square_lattice(spacing::Real[,f_cond];kwargs...)`
# Summary
Generate a square lattice of points (with equal spacing among x and y directions)
# Arguments
- `spacing` → spacing between points on both x and y axis

See [`generate_regular_lattice`](@ref) for a description of `f_cond` and of  the keyword arguments
"""
generate_square_lattice(spacing::Real,args...;kwargs...) = generate_regular_lattice(spacing,spacing,0,args...;kwargs...)

# ╔═╡ 2ba73b51-ecb9-4632-9f39-bdaeb8c5bd34
# ╠═╡ skip_as_script = true
#=╠═╡
@benchmark generate_square_lattice(1, (x,y) -> x^2 + y^2 < 100)
  ╠═╡ =#

# ╔═╡ f589306c-919d-468a-a0fd-9367acc36a7b
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## generate\_hex\_lattice
"""
  ╠═╡ =#

# ╔═╡ 6ca05079-4c0d-4c45-8486-a4291310189d
"""
`generate_hex_lattice(spacing::Real[,f_cond];kwargs...)`
# Summary
Generate a hexagonal lattice of points (with equal distance between neighboring points).
The hexagonal lattice generated by this function has distance between points on the same
column √3 times greater than the distance between points on the same row.
# Arguments
- `spacing` → spacing between points

See [`generate_regular_lattice`](@ref) for a description of `f_cond` and of  the keyword arguments
"""
generate_hex_lattice(spacing::Real,args...;kwargs...) = generate_regular_lattice(spacing .* (1,√3/2,.5)...,args...;kwargs...)

# ╔═╡ 243621c4-245a-4267-9bb4-568e673450fa
# ╠═╡ skip_as_script = true
#=╠═╡
generate_hex_lattice(1;M=20, x0 = .5) |> x -> scatter(x; mode="markers") |> Plot
  ╠═╡ =#

# ╔═╡ 059045e0-9acc-438d-b3f5-602f8d5892f7
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
# Misc Functions
"""
  ╠═╡ =#

# ╔═╡ 2f1f02c5-3ea5-40c1-8fae-704d150036e6
# Get the conversion from linear to db and viceversa
"""
$(TYPEDSIGNATURES)
Convert a number from linear to dB
"""
lin2db(x::Real) = 10log10(x)

# ╔═╡ 5be397fe-a531-423c-8be0-5d31df79dd2f
"""
$(TYPEDSIGNATURES)
Convert a number from dB to linear
"""
db2lin(x::Real) = 10^(x/10)

# ╔═╡ b2e80c33-bbfe-43ca-8795-c9d8d6fa52a9
# Convert between frequency and wavelength
"""
$(TYPEDSIGNATURES)
Get the wavelength (in m) starting from the frequency (in Hz)
"""
f2λ(f::Real) = c₀/f

# ╔═╡ 9165c4d4-69b5-456c-813c-4725feeb5b52
"""
$(TYPEDSIGNATURES)
Get the frequency (in Hz) starting from the wavelength (in m) 
"""
λ2f(λ::Real) = c₀/λ

# ╔═╡ a5ca5a8a-8497-41e2-9af0-92db5db9ce73
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
# Generate Colors
"""
  ╠═╡ =#

# ╔═╡ 0ddf072d-009d-42f2-9a8f-f69fbab750c6
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
The idea here is to find the optimal color ordering that starting from the color associated with the lattice point in (0,0) and populate the available colors in order to always selected the *un-picked* color that has the highest distance from the currently active colors (*picked* lattice points).

To do so, the idea is to verify the distane not just with the unique lattice points, but also with the neighbor lattice points on the direction of the center of the lattice parallelpiped.
To find the direction of the center of the parallelepiped, we subtract half of the lattice generating vectors from each point and compute the resulting point projection on the base of the parallelepipeid. The sign of the coefficients gives 4 possible quadrant which identify where to extend the lattice points for the computation of the closest active beam. 
"""
  ╠═╡ =#

# ╔═╡ 5db224d4-7379-4c8f-bcee-9cf00011d286
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## plot\_color\_basis
"""
  ╠═╡ =#

# ╔═╡ 5f62c7d2-ebfb-43e6-b9e3-06ca78c99390
"""
$TYPEDSIGNATURES
Plot the two arrows of the two basis vectors identified by matrix F (one vector per column).
Plotly anotations are used for creating the arrows
"""
function plot_basis_vector(F)
	common_params = (
		ax = 0,
		ay = 0,
		axref = "x",
		ayref = "y",
		arrowhead = 1,
		arrowwidth = 2,
		xanchor = "right", # Needed to have the tail match the origin
		yanchor = "top", # Needed to have the tail match the origin
	)
	# Plot the first vector arrow
	a1 = attr(;
		x = F[1,1],
		y = F[2,1],
		common_params...
	)
	# Plot the first vector arrow
	a2 = attr(;
		x = F[1,2],
		y = F[2,2],
		common_params...
	)
	Plot(Layout(annotations = [a1, a2]))
end

# ╔═╡ 6107a3dc-26dd-4a0d-aeff-5eca2cd1dcd4
# ╠═╡ skip_as_script = true
#=╠═╡
@bind NN Slider(4:20)
  ╠═╡ =#

# ╔═╡ 0d33b162-cc06-4c3b-8ade-5b40106dec0e
# ╠═╡ skip_as_script = true
#=╠═╡
NN
  ╠═╡ =#

# ╔═╡ 063db114-1b95-47d8-8f8b-26eaff8f9574
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## adjugate
"""
  ╠═╡ =#

# ╔═╡ 8d4eb806-8af1-4127-8f72-6dd68f810eb5
"""
	adjugate(A::SMatrix{2,2,<:Number,4})
Compute the [adjugate matrix](https://en.wikipedia.org/wiki/Adjugate_matrix) for a 2x2 `SMatrix`.

Used to compute the integer vector modulo in 2D space.

See: [`mod`](@ref)
"""
adjugate(A::T) where T <: SMatrix{2, 2, <:Number, 4} = T(A[4], -A[2], -A[3], A[1])

# ╔═╡ d175246c-552a-4f0f-8415-2339e9af833c
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Base.mod
"""
  ╠═╡ =#

# ╔═╡ 087730c6-cf72-4133-8064-d5619ea4b188
"""
# Integer Vector Modulo
	mod(m::SVector{2,<:Integer}, M::SMatrix{2,2,<:Integer,4})
Compute the integer modulo vector operation of an integer 2D vector with respect to an integer 2x2 matrix. Taken from `[1]`, explots the adjugate matrix instead of the inverse to avoid rounding problems

See: [`adjugate`](@ref)

References: 
- [[1]](https://doi.org/10.1109/TSP.2020.3023584) L. Xiao, X. -G. Xia and Y. -P. Wang, "Exact and Robust Reconstructions of Integer Vectors Based on Multidimensional Chinese Remainder Theorem (MD-CRT)," in IEEE Transactions on Signal Processing, vol. 68, pp. 5349-5364, 2020, doi: 10.1109/TSP.2020.3023584.
"""
function Base.mod(m::SVector{2,<:Integer}, M::SMatrix{2,2,<:Integer,4})
	M * mod.(adjugate(M)m, det(M))/det(M) |> SVector{2,Int}
end;

# ╔═╡ 745bb1c1-a312-4dbd-a29f-0836b7dbe8a7
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## get\_deterministic\_color\_order
"""
  ╠═╡ =#

# ╔═╡ 76275b93-5668-4ac9-a6a2-2f4b30ca8ab3
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## get\_color\_illumination\_order
"""
  ╠═╡ =#

# ╔═╡ 9ccf08d4-2ffb-4f36-a632-3b0ed4017d92
# ╠═╡ skip_as_script = true
#=╠═╡
let
	a = rand(SVector{2,Float64},100)
	@benchmark map(x -> Int.(sign.(x)), $a)
end
  ╠═╡ =#

# ╔═╡ f700edf6-0c44-4507-bbf5-4c5dc02fa74c
# ╠═╡ skip_as_script = true
#=╠═╡
NNN = 16
  ╠═╡ =#

# ╔═╡ 019c9ab2-ef2b-4677-a5aa-d0363fffef72
# ╠═╡ skip_as_script = true
#=╠═╡
n5 = 17
  ╠═╡ =#

# ╔═╡ cffd2b81-66c1-4f50-948c-e38ec011105d
# ╠═╡ skip_as_script = true
#=╠═╡
N4 = 4
  ╠═╡ =#

# ╔═╡ 7785e6b1-881e-4088-82fe-3dad106b07be
# ╠═╡ skip_as_script = true
#=╠═╡
n4 = 4
  ╠═╡ =#

# ╔═╡ 9cd5609b-1e58-4d39-87ab-b3d7542de691
# ╠═╡ skip_as_script = true
#=╠═╡
n6, nc = 8, 13
  ╠═╡ =#

# ╔═╡ 89888da7-4912-4557-a375-0150f80ee703
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## new\_color\_order
"""
  ╠═╡ =#

# ╔═╡ f942e39c-0d93-4bc8-8ff2-112fed566014
"""
$TYPEDSIGNATURES
Find the new optimal color ordering when starting from color `n` rather than from 1.
# Arguments
- `n` → New starting color.
- `color_order_dict` → original optimal color order for the current coloring scheme obtained via [`get_color_illumination_order`](@ref).
- `F` → lattice generating matrix for the current coloring scheme, obtained using [`generate_F_reuse_matrix`](@ref).
"""
function new_color_order(n::Int, color_order_dict::ColorOrderDict, F)
	@assert n <= length(color_order_dict) "The provided color is higher than the length of the `color_order_dict`"
	# Find the point lattice point associated to the current color
	p_current = filter(x -> x == n, color_order_dict).indices.values[1]
	idict = sort(color_order_dict) |> x->Dictionary(x.values, x.indices.values)
	# Find the new color
	new_order = map(idict) do p
		v = mod(p + p_current, F)
		color_order_dict[v]
	end
	return new_order.values
end

# ╔═╡ df0602a2-d278-4fae-9957-30c85b55598a
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## get\_projection
"""
  ╠═╡ =#

# ╔═╡ a86b2ff6-1967-424f-86eb-40fb4288c8b5
"""
$TYPEDSIGNATURES
Compute the generic projection coefficients (coordinates) of `v1` on the base identified by matrix `B`, composed of the two generating vectors `v2` and `v3`

Return the projection coefficients as a StaticVector of 2 elements.

References:
- [`https://math.stackexchange.com/questions/148199/equation-for-non-orthogonal-projection-of-a-point-onto-two-vectors-representing`](https://math.stackexchange.com/questions/148199/equation-for-non-orthogonal-projection-of-a-point-onto-two-vectors-representing)
"""
function get_projection(v1, B)
	# extract the basis generating vectors
	# Compute the matrix of the linear system of equations
	A = B'B
	# Compute the target vector 
	b = B'v1
	# solve the system
	x = A\b
end

# ╔═╡ 41817b1f-4c5e-4a55-93b7-520d6b71ea9c
function get_color_illumination_order(F, unique_points)
	N = length(unique_points) # Number of colors
	# We preallocate the array of the distances
	check_vec = Vector{eltype(unique_points)}(undef, N) # +3 because the origin is replicated 3 times
	# Preallocate the output to be modified
	out = Dictionary(unique_points, ones(Int,N))
	pts = deepcopy(unique_points)
	# Put the expanded origin points
	check_vec[1] = popat!(pts,1)
	@views @inbounds for n = 2:N
		# Find distance between all remaining points and currently populated points
		maxdist = 0.0
		id = 0
		for (i,p) ∈ enumerate(pts)
			mindist = Inf
			# Find the projection of on the generator axis
			h = get_projection(p, F)
			α, β = -Int(sign(h[1])), -Int(sign(h[2]))
			# Here we don't only check p but also the relevant closest lattice points
			for (nn,v) ∈ enumerate((SA[0,0], α*F[:,1], β*F[:,2], F * SA[α, β]))
				pp = v + p
				for c ∈ check_vec[1:n-1]					
					dist = norm(pp-c)
					mindist = dist < mindist ? dist : mindist
				end
			end
			# Find the index of the point that has the maximum minimum distance
			maxdist, id = mindist > maxdist ? (mindist, i) : (maxdist, id)
		end
		# Put that vector in the check_vec and change it's order in the dict
		pt = popat!(pts, id)
		check_vec[n] = pt
		set!(out, pt, n)
	end
	out
end

# ╔═╡ 3f2a31d4-0fa8-40fa-9dc4-bd6a26d2ddc9
# Initialize the vector that contains the matrix to compute the beam coloring. We limit ourselves at 500 colors to start
const F_reuse_matrix = (square = SMatrix{2,2,Int,4}[], triangular = SMatrix{2,2,Int,4}[])

# ╔═╡ dbef7000-7c39-49f7-b24e-f0fb436eb54e
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## compute\_F\_cell
"""
  ╠═╡ =#

# ╔═╡ 9183aa40-9dc1-4237-8bf5-de42c93b149f
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## \_coloring\_inner\_bruteforce!
"""
  ╠═╡ =#

# ╔═╡ 6b9beb62-dc7e-4b8b-9b7c-8fee5b1da98f
function _coloring_inner_bruteforce!(T_mat, rot_mat, grid_max)
    max_colours = length(T_mat)
    check_vec = fill((typemax(Int),typemax(Int)),max_colours)
    @inline norm2(x) = sum(abs2.(x))
    @inbounds for x1 = 0:grid_max, y1 = grid_max:-1:-grid_max, x2 = 0:grid_max, y2 = 0:grid_max
        mat = @SMatrix [x1 y1;x2 y2]
        # Compute the determinant
        t_det = Int(det(mat))
        # Skip points which are not likely to give useful results
        if (t_det < 1)  || (t_det > max_colours) || (t_det > grid_max^2) || (maximum(abs.(mat)) > ceil(sqrt(t_det) + 3))
            continue
        end
        # Compute the angle between the basis vectors identified by [x1,y1] and [x2,y2]
        angle = abs(atan(y1,x1) - atan(y2,x2))
        # Skip cases where the angle between vectors is either too acute or too obtuse
        if abs(π/2-angle) > π/4
            continue
        end
        # Create temp variables for computation of the minimum distance
        dmat = mat*rot_mat
        # Compute frobenius norm and minimum squared distance for the candidate lattice generating matrix
        frobe = round(Int,norm2(dmat))
        # display(frobe)
        # Minimum squared distance is either the modulo of one of the vectors or the modulo of sum or difference of them
        dmin = round(Int,minimum((norm2(dmat[1,:]), norm2(dmat[2,:]), norm2(sum(dmat;dims=1)), norm2(diff(dmat;dims=1)))))
        # Check if the current realization is better than the saved one
        if isless((-dmin,frobe),check_vec[t_det])
            # Update the check_vec
            check_vec[t_det] = (-dmin,frobe)
            # Update the vector containing the generating matrices
            T_mat[t_det] = round.(Int, mat')
        end
    end
end

# ╔═╡ 14cb2a0b-2ea8-471b-987f-1647f1516992
## Here we have the functions for the coloring computation
function compute_F_cell(max_colours::Int;grid_max::Int=ceil(Int,sqrt(max_colours))+5)
    #=
    This function is used to compute all the possible 2x2 lattice generating matrices for possible coloring schemes up to 'max_colours' colors
    Computation is done with a brute-force approach, generating all possible 2x1 vectors with maximum elements up to grid_max.

    The very heuristic automatic number of the grid_max is valid up to 5000 colors
    =#

    @assert max_colours <= 5000 "A number of colors greater than 5000 is currently not supported"
    
    # Find the current length of the pre-computed vector
    current_length = length(F_reuse_matrix.square)
    if current_length >= max_colours
        # We already computed the function for the required number of colors
        return
    end
    n_missing = max_colours - current_length
    append!(F_reuse_matrix.square,Vector{SMatrix{2,2,Int,4}}(undef,n_missing))
    append!(F_reuse_matrix.triangular,Vector{SMatrix{2,2,Int,4}}(undef,n_missing))
    # Compute the matrix for the square lattice
    _coloring_inner_bruteforce!(F_reuse_matrix.square,I,grid_max)
    # Compute the matrix for the triangular lattice
    _coloring_inner_bruteforce!(F_reuse_matrix.triangular,@SMatrix([1 0;cosd(60) sind(60)]),grid_max)
end

# ╔═╡ 1155e836-99d0-4cc8-83d0-355a6ab6fcc0
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## generate\_F\_reuse\_matrix
"""
  ╠═╡ =#

# ╔═╡ 1d44cf1c-11a5-4366-94f3-85b695c6ca12
function generate_F_reuse_matrix(N_Colours::Int=4; lattice_type::Symbol=:triangular, max_colours::Int=max(10,N_Colours))
	@assert lattice_type ∈ (:triangular, :square) "The lattice type has to be either :triangular or :square"
    compute_F_cell(max_colours)
    return getproperty(F_reuse_matrix,lattice_type)[N_Colours]
end

# ╔═╡ a7037a49-c2cb-48b6-bbf9-ca8150afcfbe
# ╠═╡ skip_as_script = true
#=╠═╡
FF = generate_F_reuse_matrix(:triangular, NN)
  ╠═╡ =#

# ╔═╡ 259ea307-2862-42f3-866a-be8f4eb83cf3
# ╠═╡ skip_as_script = true
#=╠═╡
plot_basis_vector(FF)
  ╠═╡ =#

# ╔═╡ e6033b3e-51fa-4eeb-87ff-e5b5359aebc5
"""
$TYPEDSIGNATURES
Assign a given set of unique of N integer 2D vectors, assign a deterministic and consistent order for associating these vectors to indices (1:N) used for coloring.

The order is obtained by sorting the vectors first by the norm and then by the angle it makes with respect to the X-axis (first dimension)
"""
function get_deterministic_color_order(N::Int; lattice_type = :triangular)
	T = SVector{2,Int}
	# Get the F_matrix
	F = generate_F_reuse_matrix(N; lattice_type)
	# Find the range over which to generate the potential lattice points
	M = [F*SA[1,1] F]
	xstart, xend = extrema(M[1,:])
	ystart, yend = extrema(M[2,:])
	# pts = ((x,y) for x ∈ xstart:xend, y ∈ ystart:yend)
	# Create the vector of possible values
	V = zeros(T, N)
	# Start from n = 1 as 0,0 is for sure in the set
	n = 1
	@inbounds for x ∈ xstart:xend, y ∈ ystart:yend
		p = T(x,y)
		# Get the vector modulo the color generating matrix
		v = mod(p,F)
		if v ∉ V
			n += 1
			V[n] = v
			n == N && break
		end
	end
	# Check if we found all the vector
	@assert n == N "Not all unique vectors were found!"
	# We now sort this vector
	sort!(V; by = v -> (norm(v), mod(atan(reverse(v)...), 2π)))		
end
	

# ╔═╡ 272234bf-ea90-4d24-b5f4-1cd6cdaea20b
# ╠═╡ skip_as_script = true
#=╠═╡
get_deterministic_color_order(150)
  ╠═╡ =#

# ╔═╡ 5da650b5-9309-418b-92d9-3388b5385e9b
# ╠═╡ skip_as_script = true
#=╠═╡
fffvec = get_deterministic_color_order(NNN; lattice_type = :triangular)
  ╠═╡ =#

# ╔═╡ 6579518a-2130-4e25-b4eb-966e60203f17
"""
$TYPEDSIGNATURES
Get the illumination order to achieve maximum possible distance between already illuminated colors and next illuminated color.

Returns a `Dictionary` associating each unique color vector obtained from [`get_deterministic_color_order`](@ref) to the ordering color. 
"""
function get_color_illumination_order(N; lattice_type = :triangular)
	F = generate_F_reuse_matrix(N;lattice_type)
	vec = get_deterministic_color_order(N; lattice_type)
	get_color_illumination_order(F, vec), F
end

# ╔═╡ 0f605f2c-029b-4af8-8e52-2c7bc4c7626a
# ╠═╡ skip_as_script = true
#=╠═╡
get_color_illumination_order(4; lattice_type = :square)
  ╠═╡ =#

# ╔═╡ ca09052d-1fce-422d-9205-ae4e87dc4db4
# ╠═╡ skip_as_script = true
#=╠═╡
@benchmark get_color_illumination_order(50; lattice_type = :square)
  ╠═╡ =#

# ╔═╡ f97bb119-0248-4782-8cf9-cca61a666dbf
# ╠═╡ skip_as_script = true
#=╠═╡
let
	lattice_type = :triangular
	cdict, F = get_color_illumination_order(64; lattice_type)
	order = sort(cdict) |> x -> Dictionary(x.values, x.indices.values)
	pts = map(i -> order[i],1:n5)
	outline = scatter(
		[SA[0,0],F[:,1], F*SA[1,1], F[:,2], SA[0,0]];
		mode = "lines",
		line_color = "black",
	)
	data = [
		outline,
		scatter(pts; mode = "markers"),
	]
	Plot(data)
end
  ╠═╡ =#

# ╔═╡ de86f516-b65a-445e-8988-4cfc9dafd000
# ╠═╡ skip_as_script = true
#=╠═╡
FFF = generate_F_reuse_matrix(NNN; lattice_type = :triangular)
  ╠═╡ =#

# ╔═╡ 7f4b9a79-ce48-487f-be37-d884ee9db0e9
# ╠═╡ skip_as_script = true
#=╠═╡
dio = get_color_illumination_order(FFF, fffvec)
  ╠═╡ =#

# ╔═╡ 39820d4a-d05f-4dbc-8d6d-7e20b579688c
# ╠═╡ skip_as_script = true
#=╠═╡
dio.indices.values
  ╠═╡ =#

# ╔═╡ 7bbafd36-6130-4f00-ac86-63b2bc467b77
# ╠═╡ skip_as_script = true
#=╠═╡
idio = sort(dio) |> x->Dictionary(x.values, x.indices.values)
  ╠═╡ =#

# ╔═╡ 2a45ea3b-ea3e-4d5e-912d-36b0055cf307
# ╠═╡ skip_as_script = true
#=╠═╡
new_color_order(1, dio, FFF)
  ╠═╡ =#

# ╔═╡ 759d6fd6-9649-4bea-b088-1cae1647ad3d
# ╠═╡ skip_as_script = true
#=╠═╡
map(x -> mod(x + SA[2, -2], FFF), idio)
  ╠═╡ =#

# ╔═╡ 39948c5f-d472-48ce-8385-1c67b8f1580a
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## get\_L
"""
  ╠═╡ =#

# ╔═╡ f5f5a045-1c9b-40b8-b326-7d0e1336412c
"""
$TYPEDSIGNATURES
Get the lattice generating matrix as a function of the `lattice_type` and of the lattice `spacing`.

`lattice_type` has to be either `:triangular` or `:square`
"""
function get_L(spacing; lattice_type = :triangular)
	@assert lattice_type ∈ (:triangular, :square) "`lattice_type` has to be either `:triangular` or `:square`"
	if lattice_type === :triangular
		return SA[1 1/2;0 √3/2] * spacing
	else
		return SA[1 0;0 1] * spacing
	end
end

# ╔═╡ cd9ee96f-91a9-43cf-980c-6253a0c6018f
md"""
## generate\_colors
"""

# ╔═╡ 16614e6a-2ef0-4b36-913d-6fd19440b60b
asdf = let
	N = 4
	CD, F = get_color_illumination_order(N; lattice_type = :triangular)
	spacing = .5
	BC = generate_hex_lattice(spacing; M = 100)
	L = spacing * SA[1 1/2; 0 √3/2]
	(BC, CD, L, F)
end

# ╔═╡ d3e49592-cee0-4414-9b92-ce2fb66529df
"""
$TYPEDSIGNATURES
Provide a the coloring breakdown for a given set of lattice points.
# Arguments
- `BeamCenters` → Array of SVector{2,<:Number} points identifying the beam lattice centers
- `colordict` → Dictionary of the optimal color order for the given coloring scheme, obtained using [`get_color_illumination_order`](@ref)
- `L` → Beam Lattice generating matrix, function of the beam spacing and lattice type and obtained using [`get_L`](@ref)
- `F` → Color Lattice generating matrix, function of the number of colors and beam lattice type and obtained using [`generate_F_reuse_matrix`](@ref)
# Keyword Arguments
- `first_color_idx` → Provide the index of the point in BeamCenters that should be associated with color 1.
- `first_color_coord` → Provide the coordinates of the point that should be associated with color 1, if both coord and idx are provided, the idx is ignored. Defaults to the giving color 1 to the point in (0,0).

# Note
The colors are ordered in such a way that if illuminating the colors in sequential order, the next color is always the one that has the maximum minimum distance with points of the previoulsy illuminated colors (so the one most likely to create lower interference).
"""
function generate_colors(BeamCenters, colordict, L, F; first_color_idx = nothing, first_color_coord = first_color_idx === nothing ? zero(eltype(BeamCenters)) : BeamCenters[first_color_idx])
	map(BeamCenters) do p
		T = SVector{2, Int}
		# Find the integer lattice representation
		l = inv(L)*(p - first_color_coord)
		# For some reason doing broadcast round would allocate
		v = SA[round(Int,l[1]), round(Int,l[2])]
		# Get the resulting color based on the modulo value of l w.r.t F
		color = colordict[mod(v,F)]
	end
end;

# ╔═╡ 7e68054e-4268-424e-b413-ef18baf832ac
"""
    generate_colors(BeamCenters::AbstractVector,N_Colours::Int=4;lattice_type::Symbol=:triangular)   

Provide a the coloring breakdown for a given set of lattice points.
# Arguments
- `BeamCenters` → Vector of Tuple or StaticVectors expressing the U-V coordinates of each point in the lattice for which the coloring is being computed
- `N_Colours` → Number of colors to divide the lattice in. Defaults to `4`

# keyword Arguments
- `lattice_type` → Symbol that can either be `:triangular` or `:square`, idenifying the type of lattice. Defaults to `:triangular`
"""
function generate_colors(BeamCenters::AbstractVector,N_Colours::Int=4;first_color_coord=nothing,first_color_idx=nothing,precision_digits::Int=7,lattice_type::Symbol=:triangular)
    #=
    **************************************************************************
       Generate frequency colouring file
    **************************************************************************

     References:
     [1]   "On the frequency allocation for mobile radio telephone systems", C. de
           Almeida; R. Palazzo, Proceedings of 6th International Symposium on
           Personal, Indoor and Mobile Radio Communications, Year: 1995, Volume: 1, Pages: 96 - 99 vol.1
     [2]   P. Angeletti, "Simple implementation of vectorial modulo operation based
           on fundamental parallelepiped," in Electronics Letters, vol. 48, no. 3, pp. 159-160,
           February 2 2012. doi: 10.1049/el.2011.3667
     [3]   L. Xiao, X.  -G. Xia and Y. -P. Wang, "Exact and Robust
           Reconstructions of Integer Vectors Based on Multidimensional Chinese
           Remainder Theorem (MD-CRT)," in IEEE Transactions on Signal Processing,
           vol. 68, pp. 5349-5364, 2020, doi: 10.1109/TSP.2020.3023584

     Authors: Alberto Mengali, 2020, European Space Agency

    Input Arguments:
    BeamCenters:         A vector containing the U,V coordinates of the beam centers as tuples or staticvectors
    N_Colours:           An integer number specifying the number of colors to generate in the association
    first_color_coord:   A tuple or static vector containing the U,V coordinates of the beam that will have the first color
    first_color_idx:     The beam index of the beam containing the first color, either this variable of first_order_coord can be specified, not together
    precision_digits:    The number of digits to be used in the flooring operation
    =#

    # Check if either the first color coordinates or first color idx are given
    if first_color_coord !== nothing
        if first_color_idx !== nothing
            @warn "Both first_color_idx and first_color_coord were given, disregarding the idx variable"
        end
        first_color_idx = findfirst(x -> x == first_color_coord, BeamCenters)
        if first_color_idx === nothing
          error("The provided `first_color_coord` $(first_color_coord) is not part of the beam centers vector")
        end
    end
    if lattice_type ∉ (:triangular, :square)
        @error "The specified lattice type ($lattice_type) is not recognized, it should be either :triangular or :square"
    end

    n_beams = length(BeamCenters)

    # If only 1 color is requested, simply return the trivial result
    if N_Colours == 1
        Colours = ones(n_beams)
        nbeams_per_color = n_beams
        idxs = fill(fill(true,n_beams))
    end

    # Find the minimum distance in U and V
    minU_dist = minimum(diff(sort(first.(BeamCenters)) |> x -> unique!(y -> round(y;digits=precision_digits),x)))
    minV_dist = minimum(diff(sort(last.(BeamCenters)) |> x -> unique!(y -> round(y;digits=precision_digits),x)))
    if lattice_type === :triangular
        beamU_spacing = 2minU_dist
        beamV_spacing = 2minV_dist
        # Matrix to normalize the grid points in u-v into a integer grid with y(v) axis not pointing north but north-east with 60° inclination
        D = @SMatrix [1 -1/2;0 1]
    elseif lattice_type === :square
        beamU_spacing = minU_dist
        beamV_spacing = minV_dist
        D = @SMatrix [1 0;0 1]
    end
    # Get the coloring generating matrix
    F_reuse_matrix = generate_F_reuse_matrix(N_Colours;lattice_type)
    # Create the set that will contain the unique results
    unique_colors_vecs = SVector{2,Int}[]
    # Initialize the colors vector
    Colors = similar(BeamCenters, Int)
    @inbounds for (n,p₀) ∈ enumerate(BeamCenters)
        # Compute the integer beam centers
        p₀_normalized = p₀ ./ SA_F64[beamU_spacing, minV_dist]
        # Transform the beam centers from u-v coordinates in radians into integer indexes
        beam_index_vector = round.(Int,D*p₀_normalized);
        # Find the values of the beam indexes modulo F_reuse_matrix
        unique_beam_index = round.(Int,beam_index_vector .- (F_reuse_matrix*floor.(round.(inv(F_reuse_matrix)*beam_index_vector,digits=precision_digits))))
        # Check if this color has already been assigned/registered
        idx = findfirst(x -> x == unique_beam_index,unique_colors_vecs)
        if idx isa Nothing
            push!(unique_colors_vecs, unique_beam_index)
            cidx = length(unique_colors_vecs)
        else
            cidx = idx
        end
        Colors[n] = cidx
    end

    if first_color_idx !== nothing
      # Reshuffle the colors in order to have first_color_idx as color number 1
      ref_col = Colors[first_color_idx]
      @inbounds @simd for n ∈ eachindex(Colors)
        c = Colors[n]
        Colors[n] = mod(c - ref_col,N_Colours) + 1
      end
    end

    
    return Colors
end;

# ╔═╡ d8584d03-8eb4-4864-b646-a6a0656a2e12
begin
	export generate_regular_lattice, generate_square_lattice, generate_hex_lattice, generate_rect_lattice
	export generate_colors
	export f2λ, λ2f, db2lin, lin2db
	export get_color_illumination_order, get_L, new_color_order
end

# ╔═╡ 66417309-8d6d-47e9-b5ef-7fcc0cde9194
function plot_illumination(n, N, spacing; lattice_type = :triangular)
	cdict, F = get_color_illumination_order(N; lattice_type)
	# Get the max colors
	plotly_colors = PlutoUtils.color_order_64
	# Generate the lattice of the beams
	p = generate_hex_lattice(spacing; M = round(Int, 5/spacing))
	# Find the lattice generating matrix
	L = get_L(spacing; lattice_type)
	# Find the colors
	colors = generate_colors(p, cdict, L, F)
	order = group(i -> colors[i], 1:length(p))
	data = @views map(1:n) do i
		pts = p[order[i]]
		scatter(pts;mode = "markers", marker_color = plotly_colors[i], showlegend = false)
	end
	Plot(data)
end	

# ╔═╡ ba3660ba-e274-4e11-9fd1-75aafbb0a776
# ╠═╡ skip_as_script = true
#=╠═╡
let
	n = n4
	Nmax = N4 .^ [1 2 3]
	spacing = [1 .5 .25]
	P = map(Nmax, spacing) do N,sp
		plot_illumination(min(n, N), N, sp)
	end
	vcat(P...)
end
  ╠═╡ =#

# ╔═╡ 47ffd006-29b2-4774-a791-11b3f88a7aba
# ╠═╡ skip_as_script = true
#=╠═╡
let n = n6, spacing = .5
	N = 16
	lattice_type = :triangular
	cdict, F = get_color_illumination_order(N; lattice_type)
	# Get the max colors
	plotly_colors = PlutoUtils.color_order_64
	# Generate the lattice of the beams
	p = generate_hex_lattice(spacing; M = round(Int, 5/spacing))
	# Find the lattice generating matrix
	L = get_L(spacing; lattice_type)
	# Find the colors
	colors = generate_colors(p, cdict, L, F)
	order = group(i -> colors[i], 1:length(p))
	corder = new_color_order(nc, cdict, F)
	data = @views map(corder[1:n]) do i
		pts = p[order[i]]
		scatter(pts;mode = "markers", marker_color = plotly_colors[i], showlegend = false)
	end
	Plot(data)
end	
  ╠═╡ =#

# ╔═╡ c2a099fb-446c-4a39-b505-973223c38a27
# ╠═╡ skip_as_script = true
#=╠═╡
let
	N = 4
	lattice_type = :triangular
	lat = lattice_type === :square ? generate_square_lattice : generate_hex_lattice
	CD, F = get_color_illumination_order(N; lattice_type)
	spacing = 1
	BC = lat(spacing; M = 10)
	L = get_L(spacing; lattice_type)
	colors = generate_colors(BC, CD, L, F; first_color_coord = SA[-1/2,√3/2])
	# colors = generate_colors(BC, N; lattice_type)
	data = scatter(BC; mode = "markers", marker = attr(
		color = colors, 
		colorscale = predefined_colorscale(N),
	))
	Plot(data)
end
  ╠═╡ =#

# ╔═╡ 8160086a-6349-447c-87ae-880b02fa97f5
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
# Tests
"""
  ╠═╡ =#

# ╔═╡ 3ea0415c-af14-430c-bf7c-2c7d71b7a333
# ╠═╡ skip_as_script = true
#=╠═╡
let
	N_colors = 6
	lat = generate_hex_lattice(1; M = 10)
	colors = generate_colors(lat,N_colors)	
	data = scatter(lat;mode="markers", marker_color = colors)
	Plot(data)
end
  ╠═╡ =#

# ╔═╡ 65b74f4a-e2af-4caf-8719-5a59c6349bb9
# ╠═╡ skip_as_script = true
#=╠═╡
let
	N_colors = 4
	lat = generate_square_lattice(1; M = 10)
	colors = generate_colors(lat,N_colors; lattice_type = :square)	
	data = scatter(lat;mode="markers", marker_color = colors)
	Plot(data)
end
  ╠═╡ =#

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
Dictionaries = "85a47980-9c8c-11e8-2b9f-f7ca1fa99fb4"
DocStringExtensions = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlotlyBase = "a03496cd-edff-5a9b-9e67-9cda94a718b5"
PlutoUtils = "ed5d0301-4775-4676-b788-cf71e66ff8ed"
Revise = "295af30f-e4ad-537b-8983-00126c2a3abe"
SplitApplyCombine = "03a91e81-4c3e-53e1-a0a4-9c0c8f19dd66"
StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[compat]
BenchmarkTools = "~1.3.1"
Dictionaries = "~0.3.15"
DocStringExtensions = "~0.8.6"
PlotlyBase = "~0.8.18"
PlutoUtils = "~0.5.9"
Revise = "~3.3.3"
SplitApplyCombine = "~1.2.0"
StaticArrays = "~1.2.13"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.0-rc1"
manifest_format = "2.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "4c10eee4af024676200bc7752e536f858c6b8f93"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.3.1"

[[deps.Chain]]
git-tree-sha1 = "8c4920235f6c561e401dfe569beb8b924adad003"
uuid = "8be319e6-bccf-4806-a6f7-6fae938471bc"
version = "0.5.0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "2dd813e5f2f7eec2d1268c57cf2373d3ee91fcea"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.1"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "1e315e3f4b0b7ce40feded39c73049692126cf53"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.3"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "6d4fa04343a7fc9f9cb9cff9558929f3d2752717"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.0.9"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "1fd869cc3875b57347f7027521f561cf46d1fcd8"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.19.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "d08c20eef1f2cbc6e60fd3612ac4340b89fea322"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.9"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "924cdca592bc16f14d2f7006754a621735280b74"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.1.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "fb5f5316dd3fd4c5e7c30a24d50643b73e37cd40"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.10.0"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.Dictionaries]]
deps = ["Indexing", "Random"]
git-tree-sha1 = "7669d53b75e9f9e2fa32d5215cb2af348b2c13e2"
uuid = "85a47980-9c8c-11e8-2b9f-f7ca1fa99fb4"
version = "0.3.21"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.Glob]]
git-tree-sha1 = "4df9f7e06108728ebf00a0a11edee4b29a482bb2"
uuid = "c27321d9-0574-5035-807b-f59d2c89b15c"
version = "1.3.0"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.Indexing]]
git-tree-sha1 = "ce1566720fd6b19ff3411404d4b977acd4814f9f"
uuid = "313cdc1a-70c2-5d6a-ae34-0150d3930a38"
version = "1.1.1"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "b3364212fb5d870f724876ffcd34dd8ec6d98918"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.7"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "52617c41d2761cc05ed81fe779804d3b7f14fff7"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.13"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "09e4b894ce6a976c354a69041a04748180d43637"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.15"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoweredCodeUtils]]
deps = ["JuliaInterpreter"]
git-tree-sha1 = "dedbebe234e06e1ddad435f5c6f4b85cd8ce55f7"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "2.2.2"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "0044b23da09b5608b4ecacb4e5e6c6332f833a7e"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.3.2"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlotlyBase]]
deps = ["ColorSchemes", "Dates", "DelimitedFiles", "DocStringExtensions", "JSON", "LaTeXStrings", "Logging", "Parameters", "Pkg", "REPL", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "180d744848ba316a3d0fdf4dbd34b77c7242963a"
uuid = "a03496cd-edff-5a9b-9e67-9cda94a718b5"
version = "0.8.18"

[[deps.PlutoDevMacros]]
deps = ["MacroTools", "Requires"]
git-tree-sha1 = "994167def8f46d3be21783a76705228430e29632"
uuid = "a0499f29-c39b-4c5c-807c-88074221b949"
version = "0.4.5"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "8d1f54886b9037091edf146b517989fc4a09efec"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.39"

[[deps.PlutoUtils]]
deps = ["AbstractPlutoDingetjes", "Chain", "Colors", "DocStringExtensions", "Glob", "HypertextLiteral", "OrderedCollections", "PlutoDevMacros", "PlutoUI", "PrettyTables", "Reexport", "Requires", "StaticArrays", "UUIDs"]
git-tree-sha1 = "3f8dfe27dbb980ad5e83ecd641ded8eed91f3265"
uuid = "ed5d0301-4775-4676-b788-cf71e66ff8ed"
version = "0.5.9"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "dfb54c4e414caa595a1f2ed759b160f5a3ddcba5"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.3.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Revise]]
deps = ["CodeTracking", "Distributed", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "Pkg", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "4d4239e93531ac3e7ca7e339f15978d0b5149d03"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.3.3"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "a9e798cae4867e3a41cae2dd9eb60c047f1212db"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.6"

[[deps.SplitApplyCombine]]
deps = ["Dictionaries", "Indexing"]
git-tree-sha1 = "48f393b0231516850e39f6c756970e7ca8b77045"
uuid = "03a91e81-4c3e-53e1-a0a4-9c0c8f19dd66"
version = "1.2.2"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "3c76dde64d03699e074ac02eb2e8ba8254d428da"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.2.13"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "5ce79ce186cc678bbb5c5681ca3379d1ddae11a1"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.7.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╠═379613ec-0973-4000-ae8c-d7c33ddca18e
# ╠═c0a30957-4c7b-4d7b-bfa9-c2fb691a077b
# ╠═74975885-9a4e-4857-8135-9e4f69061caf
# ╠═736b0cf6-bec2-4226-8ef4-70f6a865d34a
# ╠═f8243a65-9f5e-464e-bb06-0bb4f5131b8b
# ╠═d8584d03-8eb4-4864-b646-a6a0656a2e12
# ╠═7f645e69-3334-44db-9ba1-9f2d3e0127a2
# ╠═8660a7c4-eb78-4e7c-966b-d759df7f3dfa
# ╠═71163795-9695-4f11-acc2-6e3838c8a158
# ╠═f8a53711-e07f-4b6b-84ea-803679496571
# ╠═ce22b91e-6bba-4312-a89b-1a78f84034d3
# ╠═da97848f-a7ff-4f2d-b98d-e8bf1ccc3038
# ╠═0eb5d19d-b535-4cda-a89b-26ba295e2711
# ╠═2ba73b51-ecb9-4632-9f39-bdaeb8c5bd34
# ╠═f0834b38-8efe-4e77-b0f9-47e5b7595191
# ╠═9e5db472-f96a-4acb-96ae-024b5c73a93d
# ╠═c41f9f41-d8bd-4001-98cb-2ab788404b1b
# ╠═8787134f-9d14-4329-8dda-72557e3175b8
# ╠═f589306c-919d-468a-a0fd-9367acc36a7b
# ╠═6ca05079-4c0d-4c45-8486-a4291310189d
# ╠═243621c4-245a-4267-9bb4-568e673450fa
# ╠═059045e0-9acc-438d-b3f5-602f8d5892f7
# ╠═2f1f02c5-3ea5-40c1-8fae-704d150036e6
# ╠═5be397fe-a531-423c-8be0-5d31df79dd2f
# ╠═b2e80c33-bbfe-43ca-8795-c9d8d6fa52a9
# ╠═9165c4d4-69b5-456c-813c-4725feeb5b52
# ╠═a5ca5a8a-8497-41e2-9af0-92db5db9ce73
# ╠═0ddf072d-009d-42f2-9a8f-f69fbab750c6
# ╠═5db224d4-7379-4c8f-bcee-9cf00011d286
# ╠═5f62c7d2-ebfb-43e6-b9e3-06ca78c99390
# ╠═6107a3dc-26dd-4a0d-aeff-5eca2cd1dcd4
# ╠═0d33b162-cc06-4c3b-8ade-5b40106dec0e
# ╠═a7037a49-c2cb-48b6-bbf9-ca8150afcfbe
# ╠═259ea307-2862-42f3-866a-be8f4eb83cf3
# ╠═063db114-1b95-47d8-8f8b-26eaff8f9574
# ╠═8d4eb806-8af1-4127-8f72-6dd68f810eb5
# ╠═d175246c-552a-4f0f-8415-2339e9af833c
# ╠═087730c6-cf72-4133-8064-d5619ea4b188
# ╠═745bb1c1-a312-4dbd-a29f-0836b7dbe8a7
# ╠═e6033b3e-51fa-4eeb-87ff-e5b5359aebc5
# ╠═272234bf-ea90-4d24-b5f4-1cd6cdaea20b
# ╠═76275b93-5668-4ac9-a6a2-2f4b30ca8ab3
# ╠═6579518a-2130-4e25-b4eb-966e60203f17
# ╠═9ccf08d4-2ffb-4f36-a632-3b0ed4017d92
# ╠═41817b1f-4c5e-4a55-93b7-520d6b71ea9c
# ╠═0f605f2c-029b-4af8-8e52-2c7bc4c7626a
# ╠═ca09052d-1fce-422d-9205-ae4e87dc4db4
# ╠═f700edf6-0c44-4507-bbf5-4c5dc02fa74c
# ╠═de86f516-b65a-445e-8988-4cfc9dafd000
# ╠═5da650b5-9309-418b-92d9-3388b5385e9b
# ╠═7f4b9a79-ce48-487f-be37-d884ee9db0e9
# ╠═39820d4a-d05f-4dbc-8d6d-7e20b579688c
# ╠═2a45ea3b-ea3e-4d5e-912d-36b0055cf307
# ╠═7bbafd36-6130-4f00-ac86-63b2bc467b77
# ╠═759d6fd6-9649-4bea-b088-1cae1647ad3d
# ╠═019c9ab2-ef2b-4677-a5aa-d0363fffef72
# ╠═f97bb119-0248-4782-8cf9-cca61a666dbf
# ╠═cffd2b81-66c1-4f50-948c-e38ec011105d
# ╠═7785e6b1-881e-4088-82fe-3dad106b07be
# ╠═ba3660ba-e274-4e11-9fd1-75aafbb0a776
# ╠═66417309-8d6d-47e9-b5ef-7fcc0cde9194
# ╠═9cd5609b-1e58-4d39-87ab-b3d7542de691
# ╠═47ffd006-29b2-4774-a791-11b3f88a7aba
# ╠═89888da7-4912-4557-a375-0150f80ee703
# ╠═f942e39c-0d93-4bc8-8ff2-112fed566014
# ╠═df0602a2-d278-4fae-9957-30c85b55598a
# ╠═a86b2ff6-1967-424f-86eb-40fb4288c8b5
# ╠═3f2a31d4-0fa8-40fa-9dc4-bd6a26d2ddc9
# ╠═dbef7000-7c39-49f7-b24e-f0fb436eb54e
# ╠═14cb2a0b-2ea8-471b-987f-1647f1516992
# ╠═9183aa40-9dc1-4237-8bf5-de42c93b149f
# ╠═6b9beb62-dc7e-4b8b-9b7c-8fee5b1da98f
# ╠═1155e836-99d0-4cc8-83d0-355a6ab6fcc0
# ╠═1d44cf1c-11a5-4366-94f3-85b695c6ca12
# ╠═39948c5f-d472-48ce-8385-1c67b8f1580a
# ╠═f5f5a045-1c9b-40b8-b326-7d0e1336412c
# ╠═cd9ee96f-91a9-43cf-980c-6253a0c6018f
# ╠═16614e6a-2ef0-4b36-913d-6fd19440b60b
# ╠═d3e49592-cee0-4414-9b92-ce2fb66529df
# ╠═c2a099fb-446c-4a39-b505-973223c38a27
# ╠═7e68054e-4268-424e-b413-ef18baf832ac
# ╠═8160086a-6349-447c-87ae-880b02fa97f5
# ╠═3ea0415c-af14-430c-bf7c-2c7d71b7a333
# ╠═65b74f4a-e2af-4caf-8719-5a59c6349bb9
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
