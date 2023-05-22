### A Pluto.jl notebook ###
# v0.19.25

#> custom_attrs = ["hide-enabled"]

using Markdown
using InteractiveUtils

# ╔═╡ 36f00194-59ac-4e1a-a746-f41c9057e972
begin
	using PlutoDevMacros
end

# ╔═╡ f43c934c-84c8-4c3d-b4d9-2b716753d89c
# ╠═╡ skip_as_script = true
#=╠═╡
begin
	using BenchmarkTools
	using PlutoTest
	using PlutoExtras
end
  ╠═╡ =#

# ╔═╡ f41cdadb-808d-4714-983a-b871151ff1c0
@fromparent begin
	import *
	using >.Rotations
	using >.StaticArrays
	using >.LinearAlgebra
	using >.CoordinateTransformations
	import >.SatelliteToolbox: geodetic_to_ecef, ecef_to_geodetic, wgs84_ellipsoid
	import >.Unitful: °, rad, km, m
end

# ╔═╡ d852d113-2be1-4580-92dd-bf4082d0df11
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
# Packages
"""
  ╠═╡ =#

# ╔═╡ f41cdadb-808d-4714-983a-b871151ff32f
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
# Exports
"""
  ╠═╡ =#

# ╔═╡ 059edd4a-b3b7-4db2-9ecd-ca8a36021d2e
# ╠═╡ skip_as_script = true
#=╠═╡
ExtendedTableOfContents()
  ╠═╡ =#

# ╔═╡ e43a64ba-d776-42dd-97be-2be24a2769a7
# ╠═╡ skip_as_script = true
#=╠═╡
initialize_eqref()
  ╠═╡ =#

# ╔═╡ 91045805-53e1-457a-b7d1-db5e6df5af19
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
# Load previous notebook
"""
  ╠═╡ =#

# ╔═╡ f5577c80-ffdd-44ae-bc05-2baed9de552d
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
# Helper Functions
"""
  ╠═╡ =#

# ╔═╡ b2c827b1-2177-4b81-bdea-ea89242152ea
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Rotation Matrix
"""
  ╠═╡ =#

# ╔═╡ 3cc3b232-01e8-4064-8a2a-abe14aa6e5c0
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### User-Centric
"""
  ╠═╡ =#

# ╔═╡ 00d31f8c-dd75-4d8f-83b6-d8e976b040d0
# Generic definition, the @inline here was necessary to avoid allocations, see https://discourse.julialang.org/t/dispatch-on-value-allocating/26337/11
@inline _rotation_matrix(s::Symbol,lat,lon)::RotMatrix3{Float64} = _rotation_matrix(Val(s),lat,lon)

# ╔═╡ 46730818-1bb8-4c79-8b6f-f8cf0188c918
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### Satellite-Centric
"""
  ╠═╡ =#

# ╔═╡ 965e7534-cc27-4657-b3cf-5a5b36be2a9c
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Code Generation
"""
  ╠═╡ =#

# ╔═╡ 5113cbdb-6c07-4258-9d19-2d2a6b596fcd
_origin_transformation_docstring(srcname,dstname) = """
Convert a point from $srcname coordinates to $dstname ones

# Fields
- `origin::SVector{3,Float64}`: ECEF coordinates of the reference CRS origin
- `R::RotMatrix3{Float64}`: Rotation matrix to align the source to the destination CRS axes
- `ellipsoid::Ellipsoid{Float64}`: Reference ellipsoid used for the transformation between ECEF and other coordinates
"""

# ╔═╡ 40363971-4729-435a-b3ae-515ac30634b0
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### Basic Transformation
"""
  ╠═╡ =#

# ╔═╡ 99a35555-52e7-4e45-b265-d3868da813a8
function _basic_origin_transformation(srcname,dstname,parent,docstring=_origin_transformation_docstring(srcname,dstname))
		name = Symbol(dstname,:from,srcname)
		expr = quote
		@doc $docstring
		struct $(name) <: $(parent)
			"ECEF coordinates of the CRS origin"
			origin::SVector{3,Float64}
			"Rotation matrix for the tropocentric transformation"
			R::RotMatrix3{Float64}
			"Reference ellipsoid used in the transformation"
			ellipsoid::Ellipsoid{Float64}
		end
		function $(name)(lla::LLA;ellipsoid = wgs84_ellipsoid)
			origin = ECEFfromLLA(ellipsoid)(lla)
			R = _rotation_matrix($(Meta.quot(name)),lla.lat,lla.lon)
			$(name)(origin,R,ellipsoid)
		end
		function $(name)(ecef::StaticVector{3,Float64};ellipsoid = wgs84_ellipsoid)
			lla = LLAfromECEF(ellipsoid)(ecef)
			R = _rotation_matrix($(Meta.quot(name)),lla.lat,lla.lon)
			$(name)(ecef,R,ellipsoid)
		end
	end
	return expr |> Base.remove_linenums!
end

# ╔═╡ 2222735a-0a6b-43cd-81ea-24f9288ffa59
function _full_origin_transformation(name1,name2,parent)
	block = Expr(:block)
	
	fwdname = Symbol(name1,:from,name2)
	rvsname = Symbol(name2,:from,name1)
	
	# Do the forward direction
	expr = _basic_origin_transformation(name1,name2,parent)
	push!(block.args,expr.args...)
	
	# Do the reverse direction
	expr = _basic_origin_transformation(name2,name1,parent)
	push!(block.args,expr.args...)
	
	# Do the inversions
	expr = quote
		Base.inv(t::$fwdname) = $(rvsname)(t.origin,inv(t.R),t.ellipsoid)
		Base.inv(t::$rvsname) = $(fwdname)(t.origin,inv(t.R),t.ellipsoid)
	end |> Base.remove_linenums!
	push!(block.args,expr.args...)
	
	return block
end

# ╔═╡ 9cdffed4-bc8e-4c3f-8d3e-196b687815f6
# Boilerplate code for generating a UserCentric Transformation
macro user_transformation(name1,name2)
	block = _full_origin_transformation(name1,name2,:UserCentricTransformation)	
	esc(block)
end

# ╔═╡ 4d406319-9640-47c5-915c-0e291c30bd15
# Boilerplate code for generating a SatCentric Transformation
macro sat_transformation(name1,name2)
	block = _full_origin_transformation(name1,name2,:SatCentricTransformation)	
	esc(block)
end

# ╔═╡ 2af0a4dd-bc00-4563-91ed-7ba1caf6a0d6
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
# Generic Transformations
"""
  ╠═╡ =#

# ╔═╡ 3a537233-2f08-451c-9cb4-dcd3723cd6c8
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
The transformations here do not depend on the specific position of the user or the satellite.
"""
  ╠═╡ =#

# ╔═╡ a01e22f4-02cd-4427-8762-a3db2dd15112
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## ECEF <-> LLA 
"""
  ╠═╡ =#

# ╔═╡ 7bcb5a5d-6ee3-4cde-8d3c-97699c765fd3
begin
	struct ECEFfromLLA <: CoordinateTransformations.Transformation
		ellipsoid::Ellipsoid{Float64}
	end
	struct LLAfromECEF <: CoordinateTransformations.Transformation
		ellipsoid::Ellipsoid{Float64}
	end
	
	Base.inv(t::LLAfromECEF) = ECEFfromLLA(t.ellipsoid)
	Base.inv(t::ECEFfromLLA) = LLAfromECEF(t.ellipsoid)
	
	LLAfromECEF() = LLAfromECEF(wgs84_ellipsoid)
	ECEFfromLLA() = ECEFfromLLA(wgs84_ellipsoid)
	
	function (trans::LLAfromECEF)(ecef::StaticVector{3})
		el = trans.ellipsoid
		lat,lon,alt = ecef_to_geodetic(ecef;ellipsoid=el)
		return LLA(lat,lon,alt)
	end
	
	function (trans::ECEFfromLLA)(lat::Number,lon::Number,alt::Number) 
		el = trans.ellipsoid
		ecef = geodetic_to_ecef(lat,lon,alt;ellipsoid=el)
		return ecef
	end
	(trans::ECEFfromLLA)(lla::LLA) = trans(lla.lat,lla.lon,lla.alt) 
end

# ╔═╡ e19edbc3-c268-4294-b551-f6dd6964316a
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## ENU <-> ERA
"""
  ╠═╡ =#

# ╔═╡ fcb9caa8-2ee3-469a-8bb7-d462ab4162bd
begin
	# The transformation between ERA and tropocentric is simply a transformation between spherical and cartesian coordinates. While one needs the user location to compute ENU or ERA, the conversion between the two systems (assuming the referene UT to be the same) is indpendent on the specific user location.
struct ERAfromENU <: CoordinateTransformations.Transformation end
struct ENUfromERA <: CoordinateTransformations.Transformation end
	
Base.inv(::ERAfromENU) = ENUfromERA()
Base.inv(::ENUfromERA) = ERAfromENU()
	
function (::ERAfromENU)(enu::StaticVector{3,T}) where T
	x,y,z = enu
	# If the up coordinate is negative, the target point is not valid, so we return NaN
	z < 0 && return ERA(NaN,NaN,NaN)
	r = hypot(x, y, z)
	θ = r == 0 ? 0 : acos(z/r)
	ϕ = r == 0 ? 0 : atan(y,-x) # -x because we want the angle measured from West to North
	ERA((π/2 - θ) * rad,r * m, ϕ * rad)
end
function (::ENUfromERA)(era::ERA)
	θ = π/2 - era.el
	r = era.r
	φ = era.az
	sθ,cθ = sincos(θ)
	sφ,cφ = sincos(φ)
	x = -r * sθ * cφ  # -r because we want the angle measured from West to North
	y = r * sθ * sφ 
	z = r * cθ
	# Return the ECEF coordinates
	return SVector(x,y,z)
end
end

# ╔═╡ 0ac44137-9d7f-4746-868e-ae09b628f5e0
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## UV <-> ThetaPhi
"""
  ╠═╡ =#

# ╔═╡ ef2c3b39-5487-42ec-a006-20e0794ed21e
begin

"""
	UVfromThetaPhi <: CoordinateTransformations.Transformation
Convert a 2-D point representing θ,φ pointing angles [in rad] in the uv coordinates. Theta (θ) and Phi (φ) follow the ISO convention for spherical coordinates so represent the polar and azimuth angle respectively.
	
See also [`ThetaPhifromUV`](@ref), [`XYZfromUV`](@ref), [`UVfromXYZ`](@ref)
"""
struct UVfromThetaPhi <: CoordinateTransformations.Transformation end
"""
	ThetaPhifromUV <: CoordinateTransformations.Transformation
Convert a 2-D point representing uv pointing coordinates into the same pointing represented in polar angle θ and azimuth angle φ [expressed in rad].
	
See also [`UVfromThetaPhi`](@ref), [`XYZfromUV`](@ref), [`UVfromXYZ`](@ref)
"""
struct ThetaPhifromUV <: CoordinateTransformations.Transformation end
	
Base.inv(::UVfromThetaPhi) = ThetaPhifromUV()
Base.inv(::ThetaPhifromUV) = UVfromThetaPhi()
	
function (::ThetaPhifromUV)(uv)
	u,v = uv
	# We use the ISO Physics convention, so θ is the polar angle
	θ = asin(sqrt(u^2 + v^2))
	# φ is the angle from the U (W in WND) axis going towards V (N in WND)
	φ = atan(v,u)
	return SVector(θ,φ)
end
(t::ThetaPhifromUV)(u,v) = t((u,v))

function (::UVfromThetaPhi)(θφ)
	_check_radians(θφ)
	θ,φ = θφ
	v, u = sin(θ) .* sincos(φ)
	return SVector(u,v)
end
(t::UVfromThetaPhi)(θ,φ) = t((θ,φ))
end

# ╔═╡ 42c1d4da-2ac8-4b44-92af-8c5f0a3958e9
md"""
## UV <-> XYZ
"""

# ╔═╡ 31089d3a-e122-4f4a-bf6a-33bd6a7bff3f
begin

"""
	UVfromXYZ <: CoordinateTransformations.Transformation
	(::UVfromXYZ)(xyz::Point3D)
Convert a 3-D point representing a point in cartesian coordinates for a generic CRS `XYZ` to an pointing direction in U-V, which simply corresponds to the X and Y components (respectively) of the normalized starting point in `XYZ`.

When applying this transformation to a `Point2D` a `SVector{3, Float64}` is returned.
	
```julia
xyz2uv = UVfromXYZ()
uv = xyz2uv((0,0,100))
uv == [0,0]
```

See also [`XYZfromUV`](@ref), [`UVfromThetaPhi`](@ref), [`ThetaPhifromUV`](@ref)
"""
struct UVfromXYZ <: CoordinateTransformations.Transformation end
"""
	XYZfromUV <: CoordinateTransformations.Transformation
	(::XYZfromUV)(uv::Point2D)
Convert a 2-D point representing uv pointing coordinates and a distance `r` into the corresponding 3-D point in the 3-D CRS `XYZ`.

When applying this transformation to a `Point3D` a `SVector{2, Float64}` is returned

```julia
uv2xyz = XYZfromUV()
xyz = uv2xyz((0,0), 100)
xyz == [0,0,100]
```

See also [`UVfromXYZ`](@ref), [`UVfromThetaPhi`](@ref), [`ThetaPhifromUV`](@ref)
"""
struct XYZfromUV <: CoordinateTransformations.Transformation end
	
Base.inv(::UVfromXYZ) = XYZfromUV()
Base.inv(::XYZfromUV) = UVfromXYZ()
	
function (::XYZfromUV)(uv::Point2D, r)
	u,v = uv
	w = sqrt(1 - (u^2 + v^2))
	return SVector(u,v,w) * r
end
function (::UVfromXYZ)(xyz::Point3D)
	uvw = normalize(xyz)
	u,v,w = uvw
	return SVector(u,v)
end
end

# ╔═╡ f91fbe7d-137f-4e05-a7c7-0486db54e39e
begin
	"""
	_rotation_matrix(::Union{Val{:ENUfromECEF},Val{:ERAfromECEF}},lat,lon)
	
	Compute the rotation matrix to compute the tropocentric coordinates with tropocentric origin in the point located at geodetic coordinates `lat` and `lon` expressed in radians or Unitful Angles (both `rad` and `°`)
	"""
@inline	function _rotation_matrix(::Union{Val{:ENUfromECEF},Val{:ERAfromECEF}},lat,lon)::RotMatrix3{Float64}
		# Precompute the sines and cosines
		sλ, cλ = sincos(lon)
		sφ, cφ = sincos(lat)
		
		# Generate the rotation matrix as a StaticArray
		# Rotation matrix ECEF -> ENU [https://gssc.esa.int/navipedia/index.php/Transformations_between_ECEF_and_ENU_coordinates]
		return SA_F64[
			-sλ      cλ      0
			-cλ*sφ  -sλ*sφ   cφ
			 cλ*cφ   sλ*cφ   sφ
			] |> RotMatrix
	end
	_rotation_matrix(::Union{Val{:ECEFfromENU},Val{:ECEFfromERA}},lat,lon)::RotMatrix3{Float64} = inv(_rotation_matrix(Val(:ENUfromECEF),lat,lon))
end

# ╔═╡ 17d1271f-713d-4a85-b6ef-43e2632b74cf
begin
	# Define the relevant rotation matrix
		function _rotation_matrix(::Union{Val{:ECEFfromUV},Val{:ECEFfromWND},Val{:LLAfromUV}},lat,lon)::RotMatrix3{Float64}
		# Precompute the sines and cosines
		sλ, cλ = sincos(lon)
		sφ, cφ = sincos(lat)
		
		# Generate the rotation matrix as a StaticArray
		# Rotation matrix WND -> ECEF [mod NED -> ECEF matrix from "Global Positioning Systems, Inertial Navigation, and Integration, 2nd Ed.", Mohinder, p.472]
		return SA_F64[
			 sλ -cλ*sφ -cλ*cφ
			-cλ -sλ*sφ -sλ*cφ
			 0   cφ    -sφ
			] |> RotMatrix
	end
	_rotation_matrix(::Union{Val{:UVfromECEF},Val{:WNDfromECEF},Val{:UVfromLLA}},lat,lon)::RotMatrix3{Float64} = inv(_rotation_matrix(Val(:ECEFfromUV),lat,lon))
end

# ╔═╡ 41896117-5597-40e0-b6a1-27bba86398f5
#=╠═╡
@benchmark _rotation_matrix(Val(:UVfromECEF), x, y) setup=(x=.2;y=.8)
  ╠═╡ =#

# ╔═╡ f582bd71-774b-4745-adb0-5c2bbd00d515
#=╠═╡
@benchmark Ref(_rotation_matrix(:UVfromECEF, x, y))[] setup=(x=.2;y=.8)
  ╠═╡ =#

# ╔═╡ 436a9f46-4a05-41e1-b95a-62deb6337a8d
#=╠═╡
@benchmark Ref(_rotation_matrix(:ECEFfromUV, x, y))[] setup=(x=.2;y=.8)
  ╠═╡ =#

# ╔═╡ fc816e38-ac19-40d7-a2ab-925b97b48910
#=╠═╡
@benchmark map((x,y) -> _rotation_matrix(Val(:UVfromECEF), x, y), x,y) setup=(x = rand(1000);y=rand(1000))
  ╠═╡ =#

# ╔═╡ 1b0bb6f6-648d-46c8-b45b-85fbac0b2ed9
# ╠═╡ skip_as_script = true
#=╠═╡
let
	xyz = SA_F64[0,0,100]
	UVfromXYZ()(xyz)
end
  ╠═╡ =#

# ╔═╡ 88ba47dd-5845-4c36-83bc-d02c3cabcd63
# ╠═╡ skip_as_script = true
#=╠═╡
let
	uv = (0,0)
	XYZfromUV()(uv, 100)
end
  ╠═╡ =#

# ╔═╡ 017d2193-9a59-4186-b897-04232d61e02a
md"""
## Add Angular Offset
"""

# ╔═╡ 28a20e60-c694-408d-8aee-5aa35c498878
# ╠═╡ skip_as_script = true
#=╠═╡
todeg(x) = @. $Tuple(rad2deg(x) * °)
  ╠═╡ =#

# ╔═╡ c122e5eb-bc43-429c-b663-bc3574f2d029
# ╠═╡ skip_as_script = true
#=╠═╡
function test_φ(x, y)
	x̂ = rem2pi(to_radians(x), RoundNearest)
	ŷ = rem2pi(to_radians(y), RoundNearest)
	result = abs(x̂) ≈ abs(ŷ) ≈ π || Base.isapprox(x̂, ŷ; atol=1e-10, rtol=1e-5)
	result || @info "Phi" x y x̂ ŷ Base.isapprox(x̂, ŷ; atol=1e-10, rtol=1e-5)
	result
end
  ╠═╡ =#

# ╔═╡ 50184617-bc1f-49cd-ae9d-712e15250398
begin
# Compute the rotation matrix to find offset points following the procedure in
# this stackexchnge answer:
# https://math.stackexchange.com/questions/4343044/rotate-vector-by-a-random-little-amount
function test_offset_rotation(θ, φ)
	# Precompute the sines and cosines
	sθ, cθ = sincos(θ)
	sφ, cφ = sincos(φ)

	x̂ = SA_F64[sφ^2 + cθ * cφ^2, -sφ*cφ + cθ*sφ*cφ, -sθ*cφ]
	ŷ = SA_F64[x̂[2], x̂[1], x̂[3]]
	ẑ = SA_F64[sθ*cφ, sθ*sφ, cθ]
	
	
	_R = hcat(x̂, ŷ, ẑ) # φ̂ has to change sign to maintain the right-rule axis order
	# We have to create a rotation matrix around Z that is equivalent to π/2 - φ

	return _R
end
test_offset_rotation(θφ) = test_offset_rotation(θφ...)
end

# ╔═╡ f52f8229-4320-4b17-ab8c-cfc430d4fa1b
begin
# Compute the rotation matrix to find offset points following the procedure in
# this stackexchnge answer:
# https://math.stackexchange.com/questions/4343044/rotate-vector-by-a-random-little-amount
function angle_offset_rotation(θ, φ)
	# Precompute the sines and cosines
	sθ, cθ = sincos(θ)
	sφ, cφ = sincos(φ)
	
	# Compute the versors of the spherical to cartesian transformation as per
	# [Wikipedia](https://en.wikipedia.org/wiki/Spherical_coordinate_system#Integration_and_differentiation_in_spherical_coordinates)
	r̂ = SA_F64[sθ*cφ, sθ*sφ, cθ]
	θ̂ = SA_F64[cθ*cφ, cθ*sφ, -sθ]
	φ̂ = SA_F64[-sφ, cφ, 0]

	# The standard basis for spherical coordinates is r̂, θ̂, φ̂. We instead
	# want a basis that has r̂ as third vector (e.g. ẑ in normal cartesian
	# coordinates), and we want to rotate the other two vectors in a way that
	# the second vector is pointing towards Positive ŷ (i.e. similar to how ENU
	# has the second direction pointing towards North). 
	# To achieve this, we have to reorder the versor and then perform a matrix
	# rotation around the new third axis (which is r̂) by an angle that depends
	# on φ.
	# See
	# ![image](https://upload.wikimedia.org/wikipedia/commons/thumb/a/a2/Kugelkoord-lokale-Basis-s.svg/360px-Kugelkoord-lokale-Basis-s.svg.png)
	# for a reference figure of the original spherical versors.
	_R = hcat(-φ̂, θ̂, r̂) # φ̂ has to change sign to maintain the right-rule axis order
	# We have to create a rotation matrix around Z that is equivalent to π/2 - φ

	# We use the RotZ Matrix definition in
	# https://en.wikipedia.org/wiki/Rotation_matrix#Basic_rotations, but
	# remembering that: 
	# cos(π/2 - φ) = sin(φ)
	# sin(π/2 - φ) = cos(φ)
	__R = SA_F64[
		sφ -cφ 0
		cφ sφ 0
		0 0 1
	]

	# return __R, _R, __R*_R, _R*__R
	return _R*__R
end
angle_offset_rotation(θφ) = angle_offset_rotation(θφ...)
end

# ╔═╡ 9ee2004f-f50a-44c1-8bdb-290471a2617e
begin
	test_rot(θ, φ) = RotZYZ(π + φ, -θ, π - φ)
	test_rot(tp) = test_rot(tp...)
end

# ╔═╡ caad1fed-fea3-4571-a295-75aaa9929862
# ╠═╡ disabled = true
# ╠═╡ skip_as_script = true
#=╠═╡
let
	tp = (10°, 35°)
	a = @benchmark angle_offset_rotation($tp)
	b = @benchmark test_rot($tp)
	a,b
end
  ╠═╡ =#

# ╔═╡ 9cb24af1-0023-4a5b-9b01-9e1aa3ba2e6e
#=╠═╡
let
	tp = [(rand() * 180°, rand()*360°) for _ in 1:100]
	a = map(angle_offset_rotation, tp)
	b = map(test_rot, tp)
	b1 = @benchmark map(angle_offset_rotation, $tp)
	b2 = @benchmark map(test_rot, $tp)
	valid = @test all(a .≈ b)
	b1, b2, valid
end
  ╠═╡ =#

# ╔═╡ 92d05476-cbb4-44eb-bba6-055aed1a1a17
#=╠═╡
let
	tp = [(rand() * 180°, rand()*360°) for _ in 1:1000]
	tp2 = (rand() * 180°, rand()*360°)
	function f1(tp, tp2)
		θ, φ = to_radians(tp2)
		sθ, cθ = sincos(θ)
		sφ, cφ = sincos(φ)
		R = angle_offset_rotation(tp)
		perturbation = SA_F64[sθ*cφ, sθ*sφ, cθ]
		p3_out = R * perturbation
	end
	function f2(tp, tp2)
		θ, φ = to_radians(tp2)
		sθ, cθ = sincos(θ)
		sφ, cφ = sincos(φ)
		R = test_rot(tp)
		perturbation = SA_F64[sθ*cφ, sθ*sφ, cθ]
		p3_out = R * perturbation
	end
	a = map(y -> f1(y,tp2), tp)
	b = map(y -> f2(y,tp2), tp)
	b2 = @benchmark map(y -> $f2(y,$tp2), $tp)
	b1 = @benchmark map(y -> $f1(y,$tp2), $tp)
	valid = @test all(a .≈ b)
	b1, b2, valid
end
  ╠═╡ =#

# ╔═╡ 1a723bf6-686c-4c05-a876-463010285757
#=╠═╡
let
	x = (30°, 90°)
	@benchmark angle_offset_rotation($x)
end
  ╠═╡ =#

# ╔═╡ b79949d3-de02-4608-8248-89ad72e85fb9
begin
function _add_angular_offset(p₀, offset_angles, input_type::PointingType, output_type::PointingType)
	θ, φ = to_radians(offset_angles)
	sθ, cθ = sincos(θ)
	sφ, cφ = sincos(φ)
	
	θφ_in = if input_type isa ThetaPhi
		to_radians(p₀)
	else
		ThetaPhifromUV()(p₀)
	end
	R = angle_offset_rotation(θφ_in)
	perturbation = SA_F64[sθ*cφ, sθ*sφ, cθ]
	p3_out = R * perturbation
	u,v,w = p3_out
	@assert w >= 0 "The resulting point has a θ > π/2, so it is located behind the viewer"
	out = if output_type isa ThetaPhi
		ThetaPhifromUV()(u,v)
	else
		SVector{2,Float64}(u,v)
	end
end
end

# ╔═╡ b02ee0d4-3bfa-4131-90cc-2bbfef7ef586
begin
"""
	p = add_angular_offset(p₀, offset_angles; input_type = :UV, output_type = :UV)
	p = add_angular_offset(p₀, θ::ValidAngle, φ::ValidAngle = 0.0; input_type = :UV, output_type = :UV)

Compute the resulting pointing direction `p` obtained by adding an angular offset
specified by angles θ, φ (following the ISO/Physics convention for spherical
coordinates) [rad] to the starting position identified by p₀.

The input starting position `p₀` can be specified as any iterable of 2 elements
and is interpreted (depending on the value of the `input_type` kwarg) as a set
of 2D coordinates specified either as U-V or Theta-Phi [rad]. 

The output pointing `p` is always provided as a `SVector{2, Float64}` and similarly
to `p₀` can be given as either UV or ThetaPhi coordinates depending on the value
of the `output_type` kwarg.

For both `input_type` and `output_type`, the following Symbols are supported:
- `:thetaphi`, `:ThetaPhi` and `:θφ` can be used to represent pointing in ThetaPhi
- `:UV` and `:uv` can be used to represent pointing in UV

The offset angles can be provided either as 2nd and 3rd argument to the function
(2nd method) or as an iterable containing θ and φ as first and second element
respectively (1st method).

This function performs the inverse operation of [`get_angular_offset`](@ref) so the following code should return true
```julia
uv1 = (.3,.4)
uv2 = (-.2,.5)
offset = get_angular_offset(uv1, uv2; input_type = :uv, output_type = :thetaphi)
p = add_angular_offset(uv1, offset; input_type = :uv, output_type = :uv)
all(p .≈ uv2)
```

See also: [`get_angular_offset`](@ref), [`get_angular_distance`](@ref)
"""
function add_angular_offset(p₀, offset_angles; input_type = :UV, output_type = :UV)
	_add_angular_offset(p₀, offset_angles, PointingType(input_type), PointingType(output_type))
end

# Version with separate angles
add_angular_offset(p₀,θ::ValidAngle, φ::ValidAngle = 0; kwargs...) = add_angular_offset(p₀, (θ, φ); kwargs...)
end

# ╔═╡ 3bc4d363-3a6d-4cd1-b0c2-8a513b53ca55
#=╠═╡
@benchmark add_angular_offset((0.5,0), (10°, 0))
  ╠═╡ =#

# ╔═╡ 7bf3168c-7be0-4bb8-98cc-7a4cab893463
#=╠═╡
let
	u = rand(1000).* .5
	f(x) = _add_angular_offset((x,0.0), (10°, 0°), UV(), ThetaPhi())
	@benchmark map($f, $u)
end
  ╠═╡ =#

# ╔═╡ 80395e68-85be-42f5-ac6a-b7719b957da3
md"""
## Get Angular Offset
"""

# ╔═╡ 75d0a1c9-2535-4dcf-ae78-311aa2d76a80
function _get_angular_distance(p₁, p₂, input_type::PointingType, output_type::PointingType)
	if input_type isa UV && output_type isa UV
		Δu = p₂[1] - p₁[1]
		Δv = p₂[2] - p₁[2]
		return sqrt(Δu^2 + Δv^2)
	end
	uv₂, uv₁ = if input_type isa ThetaPhi
		tp2uv = UVfromThetaPhi()
		@. tp2uv(to_radians((p₁, p₂)))
	else
		SVector{2,Float64}.((p₁, p₂))
	end
	if output_type isa UV
		return norm(uv₂ - uv₁)
	end
	uv2xyz = XYZfromUV()
	p₁_xyz = uv2xyz(uv₁,1)
	p₂_xyz = uv2xyz(uv₂,1)
	return acos(p₁_xyz'p₂_xyz)
end

# ╔═╡ 61738518-f090-4551-88b4-0e36298d93f9
"""
	get_angular_distance(p₁, p₂; input_type=:UV, output_type=:thetaphi)

Compute the angular distance between the target pointing direction `p₂` and the starting pointing direction `p₁`. The two target pointings can be intepreted as either UV or ThetaPhi coordinates depending on the value of the `input_type` kwarg.

The output is a scalar that depending on the value of the `output_type` kwarg represents either the distance in UV, or the angular distance (Δθ) between `p₂` and `p₁`.

For both `input_type` and `output_type`, the following Symbols are supported:
- `:thetaphi`, `:ThetaPhi` and `:θφ` can be used to represent pointing in ThetaPhi
- `:UV` and `:uv` can be used to represent pointing in UV

## Note
All input angles are intepreted as radians unless specified as a Quantity in degrees using `°` from Unitful.

The output angles (when `output_type` is ThetaPhi) are also expressed in radians

This function is similar to [`get_angular_offset`](@ref) but has a slightly faster implementation in case only the distance is required, rather than the full 2D offset.
The following code should evaluate to true
```julia
uv1 = (.3,.4)
uv2 = (-.2,.5)
offset = get_angular_offset(uv1, uv2; input_type = :uv, output_type = :thetaphi)
Δθ = get_angular_distance(uv1, uv2; input_type = :uv, output_type = :thetaphi)
offset[1] ≈ Δθ
```

See also: [`add_angular_offset`](@ref), [`get_angular_offset`](@ref)
"""
function get_angular_distance(p₁, p₂; input_type = :UV, output_type = :thetaphi)
	_get_angular_distance(p₁, p₂, PointingType(input_type), PointingType(output_type))
end

# ╔═╡ 85d67f04-169d-43de-ade1-b163169aff74
#=╠═╡
let
	x = rand(1000)
	f(x) = get_angular_distance((x,0),(.4,0); input_type = :thetaphi, output_type = :thetaphi)
	@benchmark map($f, $x)
end
  ╠═╡ =#

# ╔═╡ ee3aa19f-317e-46f6-8da2-4792a84b7839
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
# OriginTransformations
"""
  ╠═╡ =#

# ╔═╡ 1c9a8798-0b03-4e50-952e-e615192dbd45
"""
	OriginTransformation <: CoordinateTransformations.Transformation

All `OriginTransformations` are used to transform points in the vicinity of the Earth between Coordinate Reference Systems (CRSs) that do not share the same origin.

These are subtyped into `UserCentricTransformation` and `SatCentricTransformation` depending on whether the reference origin of the transformation is located on the user or on the satellite.

Since the points are assumed to be around Earth, all `OriginTransformations` will have their CRS expressed in ECEF coordinates.

All `OriginTransformation` must have the following 3 fields:
- `origin::Svector{3,Float64}`: The SVector containing the ECEF coordinates of the CRS Origin
- `R::RotMatrix3{Float64}`: The rotation matrix that is needed to rotate between the starting CRS to the target CRS
- `ellipsoid::Ellipsoid{Float64}`: The ellipsoid that is used for computing geodetic points from the transformation
"""
abstract type OriginTransformation <: CoordinateTransformations.Transformation end

# ╔═╡ 2e07bdfa-7393-4864-be2f-35b7843f6cc8
abstract type UserCentricTransformation <: OriginTransformation end

# ╔═╡ e7a73ba7-731d-4a58-ac39-6fdebff78d7f
abstract type SatCentricTransformation <: OriginTransformation end

# ╔═╡ de4f6c24-49ca-429d-aa96-11d055027bbb
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## UserCentric Transformations
"""
  ╠═╡ =#

# ╔═╡ 1200a81e-c1fe-4be5-a514-a87fafc9e5fb
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### ECEF <-> ENU
"""
  ╠═╡ =#

# ╔═╡ 51a0cfb6-64a8-450a-b6e5-79f3de6c5725
begin
	# Define the transformations structs and constructors
	@user_transformation ECEF ENU
	
	function (trans::ECEFfromENU)(enu::StaticVector{3})
		ecef = trans.R * enu + trans.origin
	end
	function (trans::ENUfromECEF)(ecef::StaticVector{3})
		enu = trans.R * (ecef - trans.origin)
	end
end

# ╔═╡ 2b32d7e7-1519-4fe9-bfa8-d3c6a57b237f
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### ERA <-> ECEF
"""
  ╠═╡ =#

# ╔═╡ de6f9a8b-efc1-4666-88fe-31005efcd06e
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
The transformations defined here allow going from the ECEF coordinates of a satellite to the elevation range and azimuth as seen from a point on ground (which is the tropocentric origin used for the transformation).

The satellite position is expected in ECEF because the altitude of a satellite in orbit above the reference ellipsoid changes with latitude (if the ellipsoid is not a sphere), so by forcing the user to provide ECEF coordinates one has to think about the transformation and there is less risk of putting the same reference orbit altitude regardless of the latitude
"""
  ╠═╡ =#

# ╔═╡ 0bff095f-534e-4342-82c2-931f75e16c18
begin
	@user_transformation ECEF ERA
	function (trans::ECEFfromERA)(era::ERA)
		ecef = trans.R * ENUfromERA()(era) + trans.origin
	end
	function (trans::ERAfromECEF)(ecef::StaticVector{3})
		era = ERAfromENU()(trans.R * (ecef - trans.origin))
	end
end

# ╔═╡ 9bcf7751-ad05-404c-8432-990b436a7634
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### compute\_sat\_position
"""
  ╠═╡ =#

# ╔═╡ f34b95e5-d906-400f-a201-b9d2bf5a1b12
compute_sat_position(llaorecef::Union{LLA, StaticVector{3}}, args...;ellipsoid = wgs84_ellipsoid, kwargs...) = compute_sat_position(ECEFfromERA(llaorecef; ellipsoid), args...;kwargs...)

# ╔═╡ 97e3be69-b480-482b-a1aa-5bf2ede10cbe
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Satellite-Centric transformations
"""
  ╠═╡ =#

# ╔═╡ 95704330-4d7b-44fd-b8c0-d1570812f619
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
The computation of the lat/long position of a point on earth given the view angle (θ,φ or u,v) from the satellite can easily be performed exploiting spherical trigonometry when assuming the earth to be a sphere.

When considering the more appropriate ellipsoid of revolution model, computations become a bit more complex but the formulation can be found in [this recent paper](https://arc.aiaa.org/doi/10.2514/1.G004156)

The paper exploits the cosine directions of the pointing from the spacecraft, expressed in the earth reference frame. These can be obtained directly from the pointing U,V coordinates from the satellite point of view by performing a rotation.

We will identify the SatView CRS with axis names ``U``, ``V`` and ``W``; with axis ``W`` pointing towards the nadir direction, axis ``V`` pointing towards North and ``U`` pointing towards West (so as to have ``UVW`` following the right-hand rule).

Similarly, we will identify the earth reference frame with axes named ``X``, ``Y`` and ``Z`` and following the standard ECEF orientation with ``Z`` exiting the north pole, ``X`` exiting the equator at the longitude of the greenwhich meridian and ``Y`` pointed according to the right-hand rule.

The rotation needed to go from SatView to ECEF coordinates is the one needed to have bring the ``UVW`` axes to coincide with the ``XYZ`` ones (except the translation to align the origins).
It easy to prove that this can be achieved by rotating ``UVW`` first around ``U`` counter-clockwise by ``α = 270° - lat_s`` (obtaining the rotated CRS ``U'V'W'``), and then around ``W'`` counter-clockwise by ``γ = 90° - lon_s`` 



Exploiting the definitions of [rotation matrices](https://en.wikipedia.org/wiki/Rotation_matrix#Basic_rotations) and remembering that for change of CRS we are dealing with [*passive*](https://en.wikipedia.org/wiki/Active_and_passive_transformation) transformations, we can define the rotation matrix to translate points defined in the SatView CRS to the earth reference frame CRS as:
$(texeq("
\\mathbf{R}_{S→E} = 
\\begin{bmatrix}
	{\\rm sin}(lon_s) & -{\\rm sin}(lat_s){\\rm cos}(lon_s)& -{\\rm cos}(lat_s){\\rm cos}(lon_s) \\
	-{\\rm cos}(lon_s) & -{\\rm sin}(lat_s){\\rm sin}(lon_s) & - {\\rm cos}(lat_s){\\rm sin}(lon_s) \\
	0 & {\\rm cos}(lat_s) & -{\\rm sin}(lat_s)
\\end{bmatrix}
"))

"""
  ╠═╡ =#

# ╔═╡ 2e788b78-e5e0-4f60-aa8c-ad4f203c982e
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### Earth Intersection
"""
  ╠═╡ =#

# ╔═╡ 2af585a1-05d0-4b5a-9ee6-15eabb40a27c
function _intersection_solutions(pointing_ecef,sat_ecef,a,b)
	# Create the vector containing the parameters of the ellipse
	ellps_coeffs = SA_F64[b,b,a]
	
	# Create the vectors used to conveniently represent the 2nd degree equation of t in the paper (equation above (38))
	v1 = pointing_ecef .* ellps_coeffs
	v2 = sat_ecef .* ellps_coeffs
	
	# Find the variables to solve the quadratic equation
	α = v1'v1
	β = 2v1'v2
	γ = v2'v2 - (a*b)^2
	
	# Compute the discriminant
	Δ = β^2 - 4*α*γ
	
	# If the discriminant is negative, no intersection exists
	Δ < 0 && return NaN,NaN
	

	# Compute the two possible values of t
	t₁ = (-β - √Δ)/2α	
	t₂ = (-β + √Δ)/2α

	# t₁ < t₂ is always true
	return t₁,t₂
end

# ╔═╡ 29d1bcf5-8d65-445e-9fbb-67e742f55acb
"""
	compute_sat_position(era2ecef::ECEFfromERA, el, az; h)
	compute_sat_position(lla::LLA, el, az; h)
	compute_sat_position(ecef::StaticVector{3}, el, az; h)
Given a [`ECEFfromERA`](@ref) transformation, a target elevation `el` and azimuth `az`, and a reference height `h` [m] of the satellite from the orbit center, provides the ECEF coordinates of the satellite that is seen with angles `el` and `az` from the origin of the `era2ecef` instance at an orbit altitude of `h` **(computed from the earth/orbit center, not from the earth surface)**.

If called with an `LLA` or ECEF coordinate as the first argument, it automatically constructs the ECEFfromERA instance `era2ecef` with the provided coordinate as origin.

# Note
The values `el` and `az` are used internally to construct [`ERA`](@ref) objects, so are to be given in [rad] or specified in degrees using `°` from `Unitful`, which is also exported by TelecomUtils.

See also: [`ERA`](@ref), [`ECEF`](@ref), [`LLA`](@ref), [`ERAfromECEF`](@ref), [`ECEFfromERA`](@ref)
"""
function compute_sat_position(era2ecef::ECEFfromERA, el, az; h)
	# Find the ecef pointing, we simply compute the ecef coordinate of an ERA with unitary radius and subtract the origin coordinate of the era2ecef transformation
	pointing_ecef = era2ecef(ERA(el, 1, az)) - era2ecef.origin
	# Find the intersections with the sphere at distance h from the earth center. Since we are looking outside of earth, there are always two solution the first of which is always negative
	_, t = _intersection_solutions(pointing_ecef, era2ecef.origin, h, h)
	sat_ecef = era2ecef(ERA(el, t, az))
end

# ╔═╡ 1df46c22-c2ab-4384-9436-4b45e5603ed2
# ╠═╡ skip_as_script = true
#=╠═╡
compute_sat_position(LLA(0°,0°,0), 90°, 0°; h = 7e6) 
  ╠═╡ =#

# ╔═╡ 4cea8d15-9bb9-455c-b8bf-10b8d9a2d4af
begin
# Get the ECEF coordinates of the point where the direction of view from the satellite intercept the earth. The optional kwarg h can be provided to find the intersection at an altitude h above the ellipsoid
function earth_intersection(pointing_ecef,sat_ecef,a,b, ::ExtraOutput; h = 0.0)
	
	t₁,t₂ = _intersection_solutions(pointing_ecef,sat_ecef,a + h,b + h)

	# By definition we have t₂ > t₁ so we can already find the lower positive value of the two
	t = t₁ > 0 ? t₁ : t₂
	
	# If no solution exists, t₁ is NaN while if both solution are negative we
	# assume there is no solution as both intersection are not in direction of
	# the provided pointing. In both cases we return a 3d NaN vector
	if isnan(t) || t < 0 
		return SA_F64[NaN,NaN,NaN], NaN
	end

	# When we reach this point, if a non-zero h was provided, and a solution was found, we also have to make sure that the solution at higher altitude is not actually blocked by the original ellipsoid (earth)
	if h ≠ 0
		t̃₁, t̃₂ =  _intersection_solutions(pointing_ecef,sat_ecef,a,b)
		t̃ = t̃₁ > 0 ? t̃₁ : t̃₂
	
		# If we found two positive solution and the smaller is lower than t, it means that the earth is actually blocking the view to the candidate point, so again we return NaN
		!isnan(t̃) && t̃ > 0 && t̃ < t && t̃₁ ≠ t̃₂ && return SA_F64[NaN,NaN,NaN], NaN
	end
	
	# Compute the ecef coordinates of the intersectinon on earth
	ecef = sat_ecef + t*pointing_ecef, t
end
earth_intersection(pointing_ecef,sat_ecef,a,b; kwargs...) = earth_intersection(pointing_ecef,sat_ecef,a,b, ExtraOutput(); kwargs...)[1]
end

# ╔═╡ ea3e2a47-de2f-4383-8f01-e8fbebbdd605
#=╠═╡
let
	sat_ecef = SA_F64[1e7, 0, 0]
	pointing_ecef = SA_F64[-1, 0, 0]
	sp = SphericalEllipsoid()
	@benchmark earth_intersection($pointing_ecef, $sat_ecef, $(sp.a), $(sp.b))
end	
  ╠═╡ =#

# ╔═╡ d788faa8-df04-4a14-bef0-d76f85a9175e
#=╠═╡
let
	sat_ecef = SA_F64[1e7, 0, 0]
	pointing_ecef = SA_F64[-1, 0, 0]
	sp = SphericalEllipsoid()
	@benchmark earth_intersection($pointing_ecef, $sat_ecef, $(sp.a), $(sp.b); h = 100)
end	
  ╠═╡ =#

# ╔═╡ 5cea5fed-1cee-41f3-bcdf-2d81e96c72d4
#=╠═╡
@benchmark _intersection_solutions(SA_F64[100,0,0], SA_F64[-1,0,0], 10,10)
  ╠═╡ =#

# ╔═╡ 7ab00d88-9f0c-4ad9-a735-6ef845055823
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### ECEF <-> WND
"""
  ╠═╡ =#

# ╔═╡ f634d5d0-bb61-4bd6-9b1c-df75399de739
begin
	# Define the transformations structs and constructors
	@sat_transformation ECEF WND
	
	function (trans::ECEFfromWND)(wnd::StaticVector{3})
		ecef = trans.R * wnd + trans.origin
	end
	# Tuple overload
	(trans::ECEFfromWND)(tup::Tuple{<:Number, <:Number, <:Number}) = trans(SVector(tup))
	
	
	function (trans::WNDfromECEF)(ecef::StaticVector{3})
		wnd = trans.R * (ecef - trans.origin)
	end
	# Tuple overload
	(trans::WNDfromECEF)(tup::Tuple{<:Number, <:Number, <:Number}) = trans(SVector(tup))
end

# ╔═╡ 8b3f7041-ce2f-4d64-a135-9403eacd6385
# ╠═╡ skip_as_script = true
#=╠═╡
WNDfromECEF(LLA(0°,180°,600km))((1e7,0,0))
  ╠═╡ =#

# ╔═╡ ea890a0d-b696-4131-87ea-202ee8199358
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### ECEF <-> UV
"""
  ╠═╡ =#

# ╔═╡ a8761c85-73d2-457c-816f-7db2c83d01e9
begin
	# Define the transformations structs and constructors
	@sat_transformation ECEF UV
	
	function (trans::ECEFfromUV)(uv::Point2D, eo::ExtraOutput;h::Real=0.0)
		# Check that the uv coordinates are valid
		uv² = sum(uv .^ 2)
		@assert uv² <= 1 "u² + v² > 1, the given uv coordinate vector is not valid"
		# Compute the 3d versor identifying the pointing direction from the satellite in the local CRS coordinates
		p̂ = SA_F64[uv..., sqrt(1 - uv²)]
		# Translate the versor in ECEF coordinates
		n̂ = trans.R * p̂
		sat_ecef = trans.origin
		a,b = trans.ellipsoid.a, trans.ellipsoid.b
		ecef, r = earth_intersection(n̂,sat_ecef,a,b,eo;h)
	end
	# SingleOutput
	(trans::ECEFfromUV)(uv::Point2D; kwargs...) = trans(uv, ExtraOutput(); kwargs...)[1]
	
	function (trans::UVfromECEF)(ecef::StaticVector{3}, ::ExtraOutput)
		# Check if the given ecef coordinate is visible from the satellite position or is obstructed from earth
		pdiff = (ecef - trans.origin)
		
		# Find the magnitude of the difference to compare with the intersection solutions
		t = norm(pdiff)
		
		# Find the intersection points with the ellipsoid
		t₁,t₂ = _intersection_solutions(pdiff./t,trans.origin,trans.ellipsoid.a,trans.ellipsoid.b)
		
		# If both t₁ and t₂ are NaN, it means that no intersection with the ellipsoid is found and so there is no earth blockage
		# If t <= t₁ also no blockage is present
		# If t > t₁ then the earth is blocking the view point so we return NaN
		
		# The 1e-3 is there because the computed distance might have some error that is usually way below one mm, and 1mm shouldn't change anything for our required precision
		t₁ > 0 && t > t₁+1e-3 && return SA_F64[NaN,NaN], NaN
		
		# Find the coordinates in the satellite CRS
		xyz = trans.R * pdiff

		# If the target is behind (so D is negative), we assume that it's not visible
		xyz[3] < 0 && return SA_F64[NaN,NaN], NaN
		
		# Find the slant range between the satellite and the point
		r = norm(xyz)
		
		# Normalize the xyz vector
		uv = SVector(xyz[1],xyz[2]) ./  r
		
		# Return both the uv coordinates and the distance to the target point
		return uv, r
	end
	# Default version without range
	(trans::UVfromECEF)(ecef::StaticVector{3}) = trans(ecef,ExtraOutput())[1]
	
	# Tuple overload
	(trans::UVfromECEF)(tup::Tuple{<:Number, <:Number, <:Number}, args...) = trans(SVector(tup),args...)
end

# ╔═╡ f5b5c788-ad21-478d-972f-5bc2d7fd2768
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
### LLA <-> UV
"""
  ╠═╡ =#

# ╔═╡ 4908872b-0894-454e-afed-0efdc0c3a84f
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
We define here the transformations to switch between the satellite point of view in UV and the geodesic coordinates (LLA) of points on or above earth.
The computation is performed accounting for a custom ellipsoid shape of the earth (defaults to the WGS84 one) and an optional target height (above the reference ellipsoid) can be provided when going from UV to LLA.
This target height is used to find the correct geodesic coordinate lat,long when extending the satellite view direction to find the intersection (the same pointing direction results in different lat,long values depending on the target height).
"""
  ╠═╡ =#

# ╔═╡ d569a837-f315-47d4-9624-1bafe9996493
begin
	# Define the transformations structs and constructors
	@sat_transformation UV LLA
	
	function (trans::LLAfromUV)(uv::Point2D,eo::ExtraOutput; h::Real=0.0)
		ecef, r = ECEFfromUV(trans.origin,trans.R,trans.ellipsoid)(uv,eo;h)
		lla = LLAfromECEF(trans.ellipsoid)(ecef)
		lla, r
	end
	# Single Output
	(trans::LLAfromUV)(uv::Point2D;kwargs...) = trans(uv,ExtraOutput();kwargs...)[1]
	
	function (trans::UVfromLLA)(lla::LLA, ex::ExtraOutput)
		ecef = ECEFfromLLA(trans.ellipsoid)(lla)
		uv, xyz = UVfromECEF(trans.origin,trans.R,trans.ellipsoid)(ecef, ex)
	end
	# Single output method
	(trans::UVfromLLA)(lla::LLA) = trans(lla,ExtraOutput())[1]
	
	# Tuple overload
	(trans::UVfromLLA)(tup::Tuple{<:Number, <:Number, <:Number},args...) = trans(LLA(tup...),args...)
end

# ╔═╡ 0194ef46-0946-41a4-8909-e3e08328b4b6
function _get_angular_offset(p₁, p₂, input_type::PointingType, output_type::PointingType)
	if input_type isa UV && output_type isa UV
		Δu = p₂[1] - p₁[1]
		Δv = p₂[2] - p₁[2]
		return SVector{2, Float64}(Δu, Δv)
	end
	uv2tp = ThetaPhifromUV()
	tp2uv = inv(uv2tp)
	uv, θφ = if input_type isa ThetaPhi
		tp = to_radians.((p₁, p₂))
		uv = tp2uv.(tp)
		uv, tp
	else
		uv = (p₁, p₂)
		tp = uv2tp.(uv)
		uv, tp
	end
	if output_type isa UV
		Δu = uv[2][1] - uv[1][1]
		Δv = uv[2][2] - uv[1][2]
		return SVector{2, Float64}(Δu, Δv)
	end
		
	
	R = angle_offset_rotation(θφ[1]) # We take p₁ as reference
	p₂_xyz = XYZfromUV()(uv[2], 1) # We create the 3D vector corresponding to p₂
	perturbation = R' * p₂_xyz
	u, v, w = perturbation
	th = acos(w)
	phi = atan(v,u)
	return SA_F64[th, phi]
end

# ╔═╡ afe59330-5f20-43a7-8ed5-a9113602e3bc
"""
	get_angular_offset(p₁, p₂; input_type=:UV, output_type=:thetaphi)
Compute the angular offset required to reach the target pointing direction `p₂` from starting pointing direction `p₁`. The two target pointings can be intepreted as either UV or ThetaPhi coordinates depending on the value of the `input_type` kwarg.

The output is provided as an SVector{2,Float64} and depending on the value of the `output_type` kwarg it represents either the Δu, Δv offset in UV, or the angular offset (θ and φ) required to reach `p₂` from `p₁`.

For both `input_type` and `output_type`, the following Symbols are supported:
- `:thetaphi`, `:ThetaPhi` and `:θφ` can be used to represent pointing in ThetaPhi
- `:UV` and `:uv` can be used to represent pointing in UV

## Note
All input angles are intepreted as radians unless specified as a Quantity in degrees using `°` from Unitful.

The output angles (when `output_type` is ThetaPhi) are also expressed in radians

This function performs the inverse operation of [`add_angular_offset`](@ref) so the following code should return true
```julia
uv1 = (.3,.4)
uv2 = (-.2,.5)
offset = get_angular_offset(uv1, uv2; input_type = :uv, output_type = :thetaphi)
p = add_angular_offset(uv1, offset; input_type = :uv, output_type = :uv)
all(p .≈ uv2)
```

Check out [`get_angular_distance`](@ref) for a slightly faster implementation in case you only require the angular distance rather than the 2D offset.

See also: [`add_angular_offset`](@ref), [`get_angular_distance`](@ref)
"""
function get_angular_offset(p₁, p₂; input_type = :UV, output_type = :thetaphi)
	_get_angular_offset(p₁, p₂, PointingType(input_type), PointingType(output_type))
end

# ╔═╡ f5577c80-ffdd-44ae-bc05-2baed9de1234
begin
	export LLAfromECEF, ECEFfromLLA, LLAfromUV, UVfromLLA, ECEFfromENU, ENUfromECEF, ERAfromENU, ENUfromERA, ERAfromECEF, ECEFfromERA, ECEFfromUV, UVfromECEF
	export compute_sat_position
	export add_angular_offset, get_angular_offset, get_angular_distance
end

# ╔═╡ 97792ee2-be00-4521-b8bb-574ce03a003a
#=╠═╡
  let
	  for i in 1:100
        tp1 = SVector{2}(rand()*90°, (rand()-.5)*360°) |> Tuple
        tp2 = SVector{2}(rand()*90°, (rand()-.5)*360°) |> Tuple
        offset = get_angular_offset(tp1, tp2; input_type=:thetaphi, output_type=:thetaphi) |> todeg
        tp_target = add_angular_offset(tp1, offset; input_type = :thetaphi, output_type = :thetaphi) |> todeg
        (tp_target[1] ≈ tp2[1] && test_φ(tp_target[2], tp2[2])) || error("The forward-reverse offset test with angles failed with $((;tp1, tp2, tp_target, offset))")
	  end
  end
  ╠═╡ =#

# ╔═╡ f2510c0d-56f3-41be-8617-a4cb81e9aba8
#=╠═╡
let
	uv1 = (.3,.4)
uv2 = (-.2,.5)
offset = get_angular_offset(uv1, uv2; input_type = :uv, output_type = :thetaphi)
D = get_angular_distance(uv1, uv2; input_type = :uv, output_type = :thetaphi)
p = add_angular_offset(uv1, offset; input_type = :uv, output_type = :uv)
@test all(p .≈ uv2) && offset[1] ≈ D
end
  ╠═╡ =#

# ╔═╡ 1b5241eb-c2ed-482f-9a99-f9c951a5e853
#=╠═╡
let
	x = rand(1000)
	f(x) = get_angular_offset((x,0),(.4,0); input_type = :thetaphi, output_type = :thetaphi)
	@benchmark map($f, $x)
end
  ╠═╡ =#

# ╔═╡ 3700617d-e48f-49df-a41a-d49fabf1c52e
#=╠═╡
let
	uv2tp = ThetaPhifromUV()
	uv1 = (.4,0)
	uv2 = (.3,.5)
	tp1 = uv2tp(uv1)
	tp2 = uv2tp(uv2)
	input_type = :uv
	output_type = :thetaphi
	dist = get_angular_distance(uv1, uv2; input_type, output_type)
	offset = get_angular_offset(uv1,uv2; input_type, output_type)
	@test dist ≈ offset[1]
end
  ╠═╡ =#

# ╔═╡ d292f0d3-6a35-4f35-a5f6-e15e1c29f0f1
md"""
# Tests
"""

# ╔═╡ 4e0ecf07-e70f-4c3b-9af1-71c35167e7a8
md"""
## ECEF <-> LLA
"""

# ╔═╡ e7a82a42-9852-4ac3-8612-938004bf24de
# ╠═╡ skip_as_script = true
#=╠═╡
@test ECEFfromLLA()(LLA(30°,45°,100km)) |> LLAfromECEF() ≈ LLA(30°,45°,100km)
  ╠═╡ =#

# ╔═╡ 89eb0e56-e1d5-4497-8de2-3eed528f6358
#=╠═╡
@benchmark $ECEFfromLLA()($LLA(10°,10°,1000km))
  ╠═╡ =#

# ╔═╡ 61ace485-dc58-42dd-a58f-1cd13e1f6444
#=╠═╡
@benchmark $LLAfromECEF()(SA_F64[1e7,1e6,1e6])
  ╠═╡ =#

# ╔═╡ 76da884a-60ff-4b24-bd1f-7d5d8824ab35
md"""
## compute\_sat\_positions
"""

# ╔═╡ b07c6df9-586e-4a4c-be16-cc4ac7b1f704
# ╠═╡ skip_as_script = true
#=╠═╡
let
	h = 6371e3 + 735e3
	lla = LLA(10°, 25°, 0km)
	era2ecef = ECEFfromERA(lla)
	el, az = 35°, 80°
	satecef = compute_sat_position(era2ecef, el, az;h)
	invera = inv(era2ecef)(satecef)
	@test invera.el ≈ el && invera.az ≈ az
end
  ╠═╡ =#

# ╔═╡ 4e2b42d8-cd4f-4e29-b519-b7139a83be02
md"""
## Earth Intersection
"""

# ╔═╡ c9c9402e-c80d-4a31-9a24-f6363be60e7c
# ╠═╡ skip_as_script = true
#=╠═╡
begin
	sp = SphericalEllipsoid()
	l2e = ECEFfromLLA(sp)
end
  ╠═╡ =#

# ╔═╡ c34baa23-8483-4626-a17f-3d46ca162934
#=╠═╡
let
	# We check the erath intersection at the tangent between satellite and earth ellipsoid
	sat_lla = LLA(0,0,600km)
	sat_ecef = l2e(sat_lla)
	eoe_scan = asin(sp.a / (sp.a + sat_lla.alt))
	u = sin(eoe_scan)
	uv2lla = LLAfromUV(sat_lla; ellipsoid = sp)
	uv2lla((u * (1-1e-8),0))
	uv2lla((u,0);h= 100)
	uv2lla((u * (1-1e-6),0);h= 100)
end
  ╠═╡ =#

# ╔═╡ f836fc27-91aa-49f6-a67c-94f8c1f4a607
md"""
## LLA <-> UV
"""

# ╔═╡ 7eff4ca6-5e48-49ac-95cc-5256f8f4e0f7
#=╠═╡
let
	sat_lla = LLA(0°, 0°, 600km)
	target_lla = [
		LLA(0°, 1°, 0km), # Right - U Negative, V 0
		LLA(1°, 0°, 0km), # Top - U 0, V Positive
		LLA(0°, -1°, 0km), # Left - U Positive, V 0
		LLA(-1°, 0°, 0km), # Bottom - U 0, V Negative
	]
	target_uv = map(target_lla) do lla
		UVfromLLA(sat_lla; ellipsoid=sp)(lla) |> normalize
	end
	@test all([
		target_uv[1] == [-1,0],
		target_uv[2] == [0,1],
		target_uv[3] == [1,0],
		target_uv[4] == [0,-1],
	])
end
  ╠═╡ =#

# ╔═╡ cbafedbf-adea-4249-b681-fc2f4816ebb9
#=╠═╡
let
	sat_lla = LLA(0°, 0°, 600km)
	target_lla = LLA(0°, 0°, 610km)
	target_uv =	UVfromLLA(sat_lla; ellipsoid=sp)(target_lla)
	@test all(isnan.(target_uv))
end
  ╠═╡ =#

# ╔═╡ c4a101fc-b7d2-41cb-9252-9fbad7811957
#=╠═╡
let
	lla2ecef = ECEFfromLLA(sp)
	sat_lla = LLA(0°, 0°, 600km)
	sat_ecef = lla2ecef(sat_lla)
	R = _rotation_matrix(Val(:ECEFfromWND), sat_lla.lat, sat_lla.lon) * RotY(180°)
	target_lla = LLA(0°, 0°, 610km)
	target_uv =	UVfromLLA(sat_ecef, R', sp)(target_lla)
	@test all(map(!isnan,target_uv))
end
  ╠═╡ =#

# ╔═╡ 84c178cb-72bb-4aae-8ce0-5284b7b4a58d
#=╠═╡
let         
	# We find the pointing that corresponds to Edge of Earth and check various combination in its vicinity
	sat_lla = LLA(0,0,600km)
	sat_ecef = l2e(sat_lla)
	eoe_scan = asin(sp.a / (sp.a + sat_lla.alt))
	u = sin(eoe_scan)
	uv2lla = LLAfromUV(sat_lla; ellipsoid = sp)
	@test !isnan(uv2lla((u * (1-eps()),0))) # We should find a solution because we are pointing slightly less than EoE
	@test isnan(uv2lla((u * (1+eps()),0))) # We should not find a solution because we are pointing slightly more than EoE
	@test !isnan(uv2lla((u * (1-eps()),0), h = 100e3)) # We should find a solution because we are looking at 100km above earth
	@test !isnan(uv2lla((u * (1+eps()),0), h = 100e3)) # We should find a solution because we are looking at 100km above earth
	@test isnan(uv2lla((u * (1-eps()),0), h = 700e3)) # We should not find a solution because we are looking at 100km above the satellite alitude and with an angle slightly lower than eoe scan, so the corresponding valid point in the pointing direction is located behind earth
	@test !isnan(uv2lla((u * (1+eps()),0), h = 700e3)) # We should find a solution because we are pointing more than eoe_scan so the earth is not blocking the view of the corresponding point
end
  ╠═╡ =#

# ╔═╡ 83a29fe9-ad7a-4e7d-ba05-e6b3ce45c0c3
#=╠═╡
let
	sat_lla = LLA(0,0,600km)
	uv2lla = LLAfromUV(sat_lla; ellipsoid = sp)
	lla2uv = inv(uv2lla)
	target_uv = SA_F64[0.1,0.1]
	target_lla, r = uv2lla(target_uv, ExtraOutput())
	uv2, r2 = lla2uv(target_lla, ExtraOutput())
	@test uv2 ≈ target_uv && r2 ≈ r
end
  ╠═╡ =#

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
PlutoDevMacros = "a0499f29-c39b-4c5c-807c-88074221b949"
PlutoExtras = "ed5d0301-4775-4676-b788-cf71e66ff8ed"
PlutoTest = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"

[compat]
BenchmarkTools = "~1.3.2"
PlutoDevMacros = "~0.5.2"
PlutoExtras = "~0.7.4"
PlutoTest = "~0.2.2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "d9a9701b899b30332bbcb3e1679c41cce81fb0e8"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.3.2"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.2+0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[Parsers]]
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "478ac6c952fddd4399e71d4779797c538d0ff2bf"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.8"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.0"

[[PlutoDevMacros]]
deps = ["HypertextLiteral", "InteractiveUtils", "MacroTools", "Markdown", "Pkg", "Random", "TOML"]
git-tree-sha1 = "2ec9ca2a56ab69334ab54c79c347a9d04afae9f5"
uuid = "a0499f29-c39b-4c5c-807c-88074221b949"
version = "0.5.3"

[[PlutoExtras]]
deps = ["AbstractPlutoDingetjes", "HypertextLiteral", "InteractiveUtils", "Markdown", "OrderedCollections", "PlutoDevMacros", "PlutoUI", "REPL"]
git-tree-sha1 = "15e75e48e51416d33bab70943923a62a0b63f137"
uuid = "ed5d0301-4775-4676-b788-cf71e66ff8ed"
version = "0.7.4"

[[PlutoTest]]
deps = ["HypertextLiteral", "InteractiveUtils", "Markdown", "Test"]
git-tree-sha1 = "17aa9b81106e661cffa1c4c36c17ee1c50a86eda"
uuid = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
version = "0.2.2"

[[PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "b478a748be27bd2f2c73a7690da219d0844db305"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.51"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "7eb1686b4f04b82f96ed7a4ea5890a4f0c7a09f1"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.0"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

[[SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[Tricks]]
git-tree-sha1 = "aadb748be58b492045b4f56166b5188aa63ce549"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.7"

[[URIs]]
git-tree-sha1 = "074f993b0ca030848b897beff716d93aca60f06a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.2"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.7.0+0"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╠═d852d113-2be1-4580-92dd-bf4082d0df11
# ╠═36f00194-59ac-4e1a-a746-f41c9057e972
# ╠═f43c934c-84c8-4c3d-b4d9-2b716753d89c
# ╠═f41cdadb-808d-4714-983a-b871151ff32f
# ╠═f5577c80-ffdd-44ae-bc05-2baed9de1234
# ╠═059edd4a-b3b7-4db2-9ecd-ca8a36021d2e
# ╠═e43a64ba-d776-42dd-97be-2be24a2769a7
# ╠═91045805-53e1-457a-b7d1-db5e6df5af19
# ╠═f41cdadb-808d-4714-983a-b871151ff1c0
# ╟─f5577c80-ffdd-44ae-bc05-2baed9de552d
# ╟─b2c827b1-2177-4b81-bdea-ea89242152ea
# ╟─3cc3b232-01e8-4064-8a2a-abe14aa6e5c0
# ╠═00d31f8c-dd75-4d8f-83b6-d8e976b040d0
# ╠═f91fbe7d-137f-4e05-a7c7-0486db54e39e
# ╠═41896117-5597-40e0-b6a1-27bba86398f5
# ╠═f582bd71-774b-4745-adb0-5c2bbd00d515
# ╠═436a9f46-4a05-41e1-b95a-62deb6337a8d
# ╠═fc816e38-ac19-40d7-a2ab-925b97b48910
# ╠═46730818-1bb8-4c79-8b6f-f8cf0188c918
# ╠═17d1271f-713d-4a85-b6ef-43e2632b74cf
# ╠═965e7534-cc27-4657-b3cf-5a5b36be2a9c
# ╠═5113cbdb-6c07-4258-9d19-2d2a6b596fcd
# ╠═40363971-4729-435a-b3ae-515ac30634b0
# ╠═99a35555-52e7-4e45-b265-d3868da813a8
# ╠═2222735a-0a6b-43cd-81ea-24f9288ffa59
# ╠═9cdffed4-bc8e-4c3f-8d3e-196b687815f6
# ╠═4d406319-9640-47c5-915c-0e291c30bd15
# ╠═2af0a4dd-bc00-4563-91ed-7ba1caf6a0d6
# ╠═3a537233-2f08-451c-9cb4-dcd3723cd6c8
# ╠═a01e22f4-02cd-4427-8762-a3db2dd15112
# ╠═7bcb5a5d-6ee3-4cde-8d3c-97699c765fd3
# ╠═e19edbc3-c268-4294-b551-f6dd6964316a
# ╠═fcb9caa8-2ee3-469a-8bb7-d462ab4162bd
# ╠═0ac44137-9d7f-4746-868e-ae09b628f5e0
# ╠═ef2c3b39-5487-42ec-a006-20e0794ed21e
# ╟─42c1d4da-2ac8-4b44-92af-8c5f0a3958e9
# ╠═31089d3a-e122-4f4a-bf6a-33bd6a7bff3f
# ╠═1b0bb6f6-648d-46c8-b45b-85fbac0b2ed9
# ╠═88ba47dd-5845-4c36-83bc-d02c3cabcd63
# ╟─017d2193-9a59-4186-b897-04232d61e02a
# ╠═28a20e60-c694-408d-8aee-5aa35c498878
# ╠═c122e5eb-bc43-429c-b663-bc3574f2d029
# ╠═50184617-bc1f-49cd-ae9d-712e15250398
# ╠═caad1fed-fea3-4571-a295-75aaa9929862
# ╠═f52f8229-4320-4b17-ab8c-cfc430d4fa1b
# ╠═9ee2004f-f50a-44c1-8bdb-290471a2617e
# ╠═9cb24af1-0023-4a5b-9b01-9e1aa3ba2e6e
# ╠═92d05476-cbb4-44eb-bba6-055aed1a1a17
# ╠═97792ee2-be00-4521-b8bb-574ce03a003a
# ╠═1a723bf6-686c-4c05-a876-463010285757
# ╠═b02ee0d4-3bfa-4131-90cc-2bbfef7ef586
# ╠═b79949d3-de02-4608-8248-89ad72e85fb9
# ╠═3bc4d363-3a6d-4cd1-b0c2-8a513b53ca55
# ╠═7bf3168c-7be0-4bb8-98cc-7a4cab893463
# ╟─80395e68-85be-42f5-ac6a-b7719b957da3
# ╠═afe59330-5f20-43a7-8ed5-a9113602e3bc
# ╠═f2510c0d-56f3-41be-8617-a4cb81e9aba8
# ╠═0194ef46-0946-41a4-8909-e3e08328b4b6
# ╠═1b5241eb-c2ed-482f-9a99-f9c951a5e853
# ╠═85d67f04-169d-43de-ade1-b163169aff74
# ╠═61738518-f090-4551-88b4-0e36298d93f9
# ╠═75d0a1c9-2535-4dcf-ae78-311aa2d76a80
# ╠═3700617d-e48f-49df-a41a-d49fabf1c52e
# ╟─ee3aa19f-317e-46f6-8da2-4792a84b7839
# ╟─1c9a8798-0b03-4e50-952e-e615192dbd45
# ╠═2e07bdfa-7393-4864-be2f-35b7843f6cc8
# ╠═e7a73ba7-731d-4a58-ac39-6fdebff78d7f
# ╠═de4f6c24-49ca-429d-aa96-11d055027bbb
# ╠═1200a81e-c1fe-4be5-a514-a87fafc9e5fb
# ╠═51a0cfb6-64a8-450a-b6e5-79f3de6c5725
# ╠═2b32d7e7-1519-4fe9-bfa8-d3c6a57b237f
# ╠═de6f9a8b-efc1-4666-88fe-31005efcd06e
# ╠═0bff095f-534e-4342-82c2-931f75e16c18
# ╠═9bcf7751-ad05-404c-8432-990b436a7634
# ╠═29d1bcf5-8d65-445e-9fbb-67e742f55acb
# ╠═f34b95e5-d906-400f-a201-b9d2bf5a1b12
# ╠═1df46c22-c2ab-4384-9436-4b45e5603ed2
# ╟─97e3be69-b480-482b-a1aa-5bf2ede10cbe
# ╟─95704330-4d7b-44fd-b8c0-d1570812f619
# ╠═2e788b78-e5e0-4f60-aa8c-ad4f203c982e
# ╠═4cea8d15-9bb9-455c-b8bf-10b8d9a2d4af
# ╠═ea3e2a47-de2f-4383-8f01-e8fbebbdd605
# ╠═d788faa8-df04-4a14-bef0-d76f85a9175e
# ╠═2af585a1-05d0-4b5a-9ee6-15eabb40a27c
# ╠═5cea5fed-1cee-41f3-bcdf-2d81e96c72d4
# ╠═7ab00d88-9f0c-4ad9-a735-6ef845055823
# ╠═f634d5d0-bb61-4bd6-9b1c-df75399de739
# ╠═8b3f7041-ce2f-4d64-a135-9403eacd6385
# ╟─ea890a0d-b696-4131-87ea-202ee8199358
# ╠═a8761c85-73d2-457c-816f-7db2c83d01e9
# ╠═f5b5c788-ad21-478d-972f-5bc2d7fd2768
# ╠═4908872b-0894-454e-afed-0efdc0c3a84f
# ╠═d569a837-f315-47d4-9624-1bafe9996493
# ╠═d292f0d3-6a35-4f35-a5f6-e15e1c29f0f1
# ╠═4e0ecf07-e70f-4c3b-9af1-71c35167e7a8
# ╠═e7a82a42-9852-4ac3-8612-938004bf24de
# ╠═89eb0e56-e1d5-4497-8de2-3eed528f6358
# ╠═61ace485-dc58-42dd-a58f-1cd13e1f6444
# ╠═76da884a-60ff-4b24-bd1f-7d5d8824ab35
# ╠═b07c6df9-586e-4a4c-be16-cc4ac7b1f704
# ╟─4e2b42d8-cd4f-4e29-b519-b7139a83be02
# ╠═c9c9402e-c80d-4a31-9a24-f6363be60e7c
# ╠═c34baa23-8483-4626-a17f-3d46ca162934
# ╟─f836fc27-91aa-49f6-a67c-94f8c1f4a607
# ╠═7eff4ca6-5e48-49ac-95cc-5256f8f4e0f7
# ╠═cbafedbf-adea-4249-b681-fc2f4816ebb9
# ╠═c4a101fc-b7d2-41cb-9252-9fbad7811957
# ╠═84c178cb-72bb-4aae-8ce0-5284b7b4a58d
# ╠═83a29fe9-ad7a-4e7d-ba05-e6b3ce45c0c3
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
