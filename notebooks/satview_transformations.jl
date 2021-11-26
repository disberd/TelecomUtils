### A Pluto.jl notebook ###
# v0.17.2

# using Markdown
# using InteractiveUtils

# ╔═╡ 36f00194-59ac-4e1a-a746-f41c9057e972
begin
	using CoordinateTransformations
	using StaticArrays
	using LinearAlgebra
	using Unitful
	using Unitful.DefaultSymbols
	using Rotations
	using Parameters
	using SatelliteToolbox: geodetic_to_ecef, ecef_to_geodetic, wgs84_ellipsoid
	using DocStringExtensions
	using Parameters
	import Proj4
end

# ╔═╡ f43c934c-84c8-4c3d-b4d9-2b716753d89c
#=╠═╡ notebook_exclusive
begin
	using BenchmarkTools
	using PlutoTest
	using PlutoUtils
	using PlutoDevMacros
end
  ╠═╡ notebook_exclusive =#

# ╔═╡ d852d113-2be1-4580-92dd-bf4082d0df11
#=╠═╡ notebook_exclusive
md"""
# Packages
"""
  ╠═╡ notebook_exclusive =#

# ╔═╡ 059edd4a-b3b7-4db2-9ecd-ca8a36021d2e
#=╠═╡ notebook_exclusive
ToC()
  ╠═╡ notebook_exclusive =#

# ╔═╡ e43a64ba-d776-42dd-97be-2be24a2769a7
#=╠═╡ notebook_exclusive
initialize_eqref()
  ╠═╡ notebook_exclusive =#

# ╔═╡ 91045805-53e1-457a-b7d1-db5e6df5af19
#=╠═╡ notebook_exclusive
md"""
## Load previous notebook
"""
  ╠═╡ notebook_exclusive =#

# ╔═╡ f41cdadb-808d-4714-983a-b871151ff1c0
@plutoinclude "satview_basics.jl" "all"

# ╔═╡ f41cdadb-808d-4714-983a-b871151ff32f
#=╠═╡ notebook_exclusive
md"""
# Exports
"""
  ╠═╡ notebook_exclusive =#

# ╔═╡ f5577c80-ffdd-44ae-bc05-2baed9de1234
export LLAfromECEF, ECEFfromLLA, LLAfromUV, UVfromLLA, ECEFfromENU, ENUfromECEF, ERAfromENU, ENUfromERA, ERAfromECEF, ECEFfromERA, ECEFfromUV, UVfromECEF

# ╔═╡ f5577c80-ffdd-44ae-bc05-2baed9de552d
#=╠═╡ notebook_exclusive
md"""
# Helper Functions
"""
  ╠═╡ notebook_exclusive =#

# ╔═╡ b2c827b1-2177-4b81-bdea-ea89242152ea
#=╠═╡ notebook_exclusive
md"""
## Rotation Matrix
"""
  ╠═╡ notebook_exclusive =#

# ╔═╡ 00d31f8c-dd75-4d8f-83b6-d8e976b040d0
# Generic definition, the @inline here was necessary to avoid allocations, see https://discourse.julialang.org/t/dispatch-on-value-allocating/26337/11
@inline _rotation_matrix(s::Symbol,lat,lon)::RotMatrix3{Float64} = _rotation_matrix(Val(s),lat,lon)

# ╔═╡ 3cc3b232-01e8-4064-8a2a-abe14aa6e5c0
#=╠═╡ notebook_exclusive
md"""
### User-Centric
"""
  ╠═╡ notebook_exclusive =#

# ╔═╡ f91fbe7d-137f-4e05-a7c7-0486db54e39e
begin
	"""
	_rotation_matrix(::Union{Val{:ENUfromECEF},Val{:ERAfromECEF}},lat,lon)
	
	Compute the rotation matrix to compute the tropocentric coordinates with tropocentric origin in the point located at geodetic coordinates `lat` and `lon` expressed in radians or Unitful Angles (both `rad` and `°`)
	"""
@inline	function _rotation_matrix(::Union{Val{:ENUfromECEF},Val{:ERAfromECEF}},lat,lon)::RotMatrix3{Float64}
		# Precompute the sines and cosines
		sλ, cλ = sincos(lat)
		sφ, cφ = sincos(lon)
		
		# Generate the rotation matrix as a StaticArray
		return SA_F64[
			-sλ      cλ      0
			-sφ*cλ  -sφ*sλ   cφ
			 cφ*cλ   cφ*sλ   sφ
			] |> RotMatrix
	end
	_rotation_matrix(::Union{Val{:ECEFfromENU},Val{:ECEFfromERA}},lat,lon)::RotMatrix3{Float64} = inv(_rotation_matrix(Val(:ENUfromECEF),lat,lon))
end

# ╔═╡ 46730818-1bb8-4c79-8b6f-f8cf0188c918
#=╠═╡ notebook_exclusive
md"""
### Satellite-Centric
"""
  ╠═╡ notebook_exclusive =#

# ╔═╡ 17d1271f-713d-4a85-b6ef-43e2632b74cf
begin
	# Define the relevant rotation matrix
		function _rotation_matrix(::Union{Val{:ECEFfromUV},Val{:ECEFfromWND},Val{:LLAfromUV}},lat,lon)::RotMatrix3{Float64}
		# Precompute the sines and cosines
		sλ, cλ = sincos(lat)
		sφ, cφ = sincos(lon)
		
		# Generate the rotation matrix as a StaticArray
		return SA_F64[
			 sφ -sλ*cφ -cλ*cφ
			-cφ -sλ*sφ -cλ*sφ
			 0   cλ    -sλ
			] |> RotMatrix
	end
	_rotation_matrix(::Union{Val{:UVfromECEF},Val{:WNDfromECEF},Val{:UVfromLLA}},lat,lon)::RotMatrix3{Float64} = inv(_rotation_matrix(Val(:ECEFfromUV),lat,lon))
end

# ╔═╡ 965e7534-cc27-4657-b3cf-5a5b36be2a9c
#=╠═╡ notebook_exclusive
md"""
## Code Generation
"""
  ╠═╡ notebook_exclusive =#

# ╔═╡ 5113cbdb-6c07-4258-9d19-2d2a6b596fcd
_origin_transformation_docstring(srcname,dstname) = """
Convert a point from $srcname coordinates to $dstname ones

# Fields
- `origin::SVector{3,Float64}`: ECEF coordinates of the reference CRS origin
- `R::RotMatrix3{Float64}`: Rotation matrix to align the source to the destination CRS axes
- `ellipsoid::Ellipsoid{Float64}`: Reference ellipsoid used for the transformation between ECEF and other coordinates
"""

# ╔═╡ 40363971-4729-435a-b3ae-515ac30634b0
#=╠═╡ notebook_exclusive
md"""
### Basic Transformation
"""
  ╠═╡ notebook_exclusive =#

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
#=╠═╡ notebook_exclusive
md"""
# Generic Transformations
"""
  ╠═╡ notebook_exclusive =#

# ╔═╡ 3a537233-2f08-451c-9cb4-dcd3723cd6c8
#=╠═╡ notebook_exclusive
md"""
The transformations here do not depend on the specific position of the user or the satellite.
"""
  ╠═╡ notebook_exclusive =#

# ╔═╡ a01e22f4-02cd-4427-8762-a3db2dd15112
#=╠═╡ notebook_exclusive
md"""
## ECEF <-> LLA 
"""
  ╠═╡ notebook_exclusive =#

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
#=╠═╡ notebook_exclusive
md"""
## ENU <-> ERA
"""
  ╠═╡ notebook_exclusive =#

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
#=╠═╡ notebook_exclusive
md"""
## UV <-> ThetaPhi
"""
  ╠═╡ notebook_exclusive =#

# ╔═╡ ef2c3b39-5487-42ec-a006-20e0794ed21e
begin

"""
	UVfromThetaPhi <: CoordinateTransformations.Transformation
Convert a 2-D point representing θ,φ pointing angles [in rad] in the uv coordinates. Theta (θ) and Phi (φ) follow the ISO convention for spherical coordinates so represent the polar and azimuth angle respectively.
	
See also [`ThetaPhifromUV](@ref)
"""
struct UVfromThetaPhi <: CoordinateTransformations.Transformation end
"""
	ThetaPhifromUV <: CoordinateTransformations.Transformation
Convert a 2-D point representing uv pointing coordinates into the same pointing represented in polar angle θ and azimuth angle φ [expressed in rad].
	
See also [`UVfromThetaPhi](@ref)
"""
struct ThetaPhifromUV <: CoordinateTransformations.Transformation end
	
Base.inv(::UVfromThetaPhi) = ThetaPhifromUV()
Base.inv(::ThetaPhifromUV) = UVfromThetaPhi()
	
function (::ThetaPhifromUV)(uv::StaticVector{2,T}) where T
	u,v = uv
	# We use the ISO Physics convention, so θ is the polar angle
	θ = asin(norm(uv))
	# φ is the angle from the U (W in WND) axis going towards V (N in WND)
	φ = atan(v,u)
	return SVector(θ,φ)
end
function (::UVfromThetaPhi)(θφ::StaticVector{2,T}) where T
	θ,φ = θφ
	v, u = sin(θ) .* sincos(φ)
	return SVector(u,v)
end
end

# ╔═╡ ee3aa19f-317e-46f6-8da2-4792a84b7839
#=╠═╡ notebook_exclusive
md"""
# OriginTransformations
"""
  ╠═╡ notebook_exclusive =#

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
#=╠═╡ notebook_exclusive
md"""
## UserCentric Transformations
"""
  ╠═╡ notebook_exclusive =#

# ╔═╡ 1200a81e-c1fe-4be5-a514-a87fafc9e5fb
#=╠═╡ notebook_exclusive
md"""
### ECEF <-> ENU
"""
  ╠═╡ notebook_exclusive =#

# ╔═╡ 2b32d7e7-1519-4fe9-bfa8-d3c6a57b237f
#=╠═╡ notebook_exclusive
md"""
### ERA <-> ECEF
"""
  ╠═╡ notebook_exclusive =#

# ╔═╡ de6f9a8b-efc1-4666-88fe-31005efcd06e
#=╠═╡ notebook_exclusive
md"""
The transformations defined here allow going from the ECEF coordinates of a satellite to the elevation range and azimuth as seen from a point on ground (which is the tropocentric origin used for the transformation).

The satellite position is expected in ECEF because the altitude of a satellite in orbit above the reference ellipsoid changes with latitude (if the ellipsoid is not a sphere), so by forcing the user to provide ECEF coordinates one has to think about the transformation and there is less risk of putting the same reference orbit altitude regardless of the latitude
"""
  ╠═╡ notebook_exclusive =#

# ╔═╡ 97e3be69-b480-482b-a1aa-5bf2ede10cbe
#=╠═╡ notebook_exclusive
md"""
## Satellite-Centric transformations
"""
  ╠═╡ notebook_exclusive =#

# ╔═╡ 95704330-4d7b-44fd-b8c0-d1570812f619
#=╠═╡ notebook_exclusive
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
  ╠═╡ notebook_exclusive =#

# ╔═╡ 2e788b78-e5e0-4f60-aa8c-ad4f203c982e
#=╠═╡ notebook_exclusive
md"""
### Earth Intersection
"""
  ╠═╡ notebook_exclusive =#

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
	
	return t₁,t₂
end

# ╔═╡ 4cea8d15-9bb9-455c-b8bf-10b8d9a2d4af
# Get the ECEF coordinates of the point where the direction of view from the satellite intercept the earth 
function earth_intersection(pointing_ecef,sat_ecef,a,b)
	
	t₁,t₂ = _intersection_solutions(pointing_ecef,sat_ecef,a,b)
	
	# If no solution exists, t₁ is NaN, so we return a 3d NaN vector
	isnan(t₁) && return SA_F64[NaN,NaN,NaN]
	
	t = abs(t₁) < abs(t₂) ? t₁ : t₂
	
	# Compute the ecef coordinates of the intersectinon on earth
	ecef = sat_ecef + t*pointing_ecef
end

# ╔═╡ 7ab00d88-9f0c-4ad9-a735-6ef845055823
#=╠═╡ notebook_exclusive
md"""
### ECEF <-> WND
"""
  ╠═╡ notebook_exclusive =#

# ╔═╡ ea890a0d-b696-4131-87ea-202ee8199358
#=╠═╡ notebook_exclusive
md"""
### ECEF <-> UV
"""
  ╠═╡ notebook_exclusive =#

# ╔═╡ f5b5c788-ad21-478d-972f-5bc2d7fd2768
#=╠═╡ notebook_exclusive
md"""
### LLA <-> UV
"""
  ╠═╡ notebook_exclusive =#

# ╔═╡ 4908872b-0894-454e-afed-0efdc0c3a84f
#=╠═╡ notebook_exclusive
md"""
We define here the transformations to switch between the satellite point of view in UV and the geodesic coordinates (LLA) of points on or above earth.
The computation is performed accounting for a custom ellipsoid shape of the earth (defaults to the WGS84 one) and an optional target height (above the reference ellipsoid) can be provided when going from UV to LLA.
This target height is used to find the correct geodesic coordinate lat,long when extending the satellite view direction to find the intersection (the same pointing direction results in different lat,long values depending on the target height).
"""
  ╠═╡ notebook_exclusive =#

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
#=╠═╡ notebook_exclusive
WNDfromECEF(LLA(0°,180°,600km))((1e7,0,0))
  ╠═╡ notebook_exclusive =#

# ╔═╡ d292f0d3-6a35-4f35-a5f6-e15e1c29f0f1
#=╠═╡ notebook_exclusive
md"""
# Tests
"""
  ╠═╡ notebook_exclusive =#

# ╔═╡ 4e0ecf07-e70f-4c3b-9af1-71c35167e7a8
#=╠═╡ notebook_exclusive
md"""
## ECEF <-> LLA
"""
  ╠═╡ notebook_exclusive =#

# ╔═╡ e7a82a42-9852-4ac3-8612-938004bf24de
#=╠═╡ notebook_exclusive
@test ECEFfromLLA()(LLA(30°,45°,100km)) |> LLAfromECEF() ≈ LLA(30°,45°,100km)
  ╠═╡ notebook_exclusive =#

# ╔═╡ 89eb0e56-e1d5-4497-8de2-3eed528f6358
#=╠═╡ notebook_exclusive
@benchmark $ECEFfromLLA()($LLA(10°,10°,1000km))
  ╠═╡ notebook_exclusive =#

# ╔═╡ 61ace485-dc58-42dd-a58f-1cd13e1f6444
#=╠═╡ notebook_exclusive
@benchmark $LLAfromECEF()(SA_F64[1e7,1e6,1e6])
  ╠═╡ notebook_exclusive =#

# ╔═╡ d569a837-f315-47d4-9624-1bafe9996493
begin
	# Define the transformations structs and constructors
	@sat_transformation UV LLA
	
	function (trans::LLAfromUV)(uv::StaticVector{2},h::Real=0.0)
		ecef = ECEFfromUV(trans.origin,trans.R,trans.ellipsoid)(uv,h)
		lla = LLAfromECEF(trans.ellipsoid)(ecef)
	end
	# Tuple overload
	(trans::LLAfromUV)(tup::Tuple{<:Number, <:Number},h::Real=0.0) = trans(SVector(tup),h)
	
	function (trans::UVfromLLA)(lla::LLA, ex::ExtraOutput)
		ecef = ECEFfromLLA(trans.ellipsoid)(lla)
		uv, wnd = UVfromECEF(trans.origin,trans.R,trans.ellipsoid)(ecef, ex)
	end
	# Single output method
	(trans::UVfromLLA)(lla::LLA) = trans(lla,ExtraOutput())[1]
	
	# Tuple overload
	(trans::UVfromLLA)(tup::Tuple{<:Number, <:Number, <:Number},args...) = trans(LLA(tup...),args...)
end

# ╔═╡ a8761c85-73d2-457c-816f-7db2c83d01e9
begin
	# Define the transformations structs and constructors
	@sat_transformation ECEF UV
	
	function (trans::ECEFfromUV)(uv::StaticVector{2},h::Real=0.0)
		# Check that the uv coordinates are valid
		uv² = sum(uv .^ 2)
		@assert uv² <= 1 "u² + v² > 1, the given uv coordinate vector is not valid"
		# Compute the 3d versor identifying the pointing direction from the satellite in WND coordinates
		p̂ = SA_F64[uv..., sqrt(1 - uv²)]
		# Translate the versor in ECEF coordinates
		n̂ = trans.R * p̂
		sat_ecef = trans.origin
		a,b = trans.ellipsoid.a, trans.ellipsoid.b
		ecef = earth_intersection(n̂,sat_ecef,a+h,b+h)
	end
	# Tuple overload
	(trans::ECEFfromUV)(tup::Tuple{<:Number, <:Number}, h=0.0) = trans(SVector(tup), h)
	
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
		!isnan(t₁) && t > t₁+1e-3 && return SA_F64[NaN,NaN], NaN
		
		# Find the coordinates in the West-North-Down CRS
		wnd = trans.R * pdiff
		
		# Find the slant range between the satellite and the point
		r = norm(wnd)
		
		# Normalize the wnd vector
		uv = SVector(wnd[1],wnd[2]) ./  r
		
		# Return both the uv coordinates and the slant range
		return uv, r
	end
	# Default version without range
	(trans::UVfromECEF)(ecef::StaticVector{3}) = trans(ecef,ExtraOutput())[1]
	
	# Tuple overload
	(trans::UVfromECEF)(tup::Tuple{<:Number, <:Number, <:Number}, args...) = trans(SVector(tup),args...)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
CoordinateTransformations = "150eb455-5306-5404-9cee-2592286d6298"
DocStringExtensions = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Parameters = "d96e819e-fc66-5662-9728-84c9c7592b0a"
PlutoDevMacros = "a0499f29-c39b-4c5c-807c-88074221b949"
PlutoTest = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
PlutoUtils = "ed5d0301-4775-4676-b788-cf71e66ff8ed"
Proj4 = "9a7e659c-8ee8-5706-894e-f68f43bc57ea"
Rotations = "6038ab10-8711-5258-84ad-4b1120ba62dc"
SatelliteToolbox = "6ac157d9-b43d-51bb-8fab-48bf53814f4a"
StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[compat]
BenchmarkTools = "~1.2.0"
CoordinateTransformations = "~0.6.2"
DocStringExtensions = "~0.8.6"
Parameters = "~0.12.3"
PlutoDevMacros = "~0.3.7"
PlutoTest = "~0.2.0"
PlutoUtils = "~0.4.13"
Proj4 = "~0.7.6"
Rotations = "~1.1.0"
SatelliteToolbox = "~0.9.3"
StaticArrays = "~1.2.13"
Unitful = "~1.9.2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.0-rc2"
manifest_format = "2.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "0bc60e3006ad95b4bb7497698dd7c6d649b9bc06"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.1"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "84918055d15b3114ede17ac6a7182f68870c16f7"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.1"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "61adeb0823084487000600ef8b1c00cc2474cd47"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.2.0"

[[deps.CEnum]]
git-tree-sha1 = "215a9aa4a1f23fbd05b92769fdd62559488d70e9"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.Chain]]
git-tree-sha1 = "cac464e71767e8a04ceee82a889ca56502795705"
uuid = "8be319e6-bccf-4806-a6f7-6fae938471bc"
version = "0.4.8"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "f885e7e7c124f8c92650d61b9477b9ac2ee607dd"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.11.1"

[[deps.ChangesOfVariables]]
deps = ["LinearAlgebra", "Test"]
git-tree-sha1 = "9a1d594397670492219635b35a3d830b04730d62"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.1"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "dce3e3fea680869eaa0b774b2e8343e9ff442313"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.40.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f74e9d5388b8620b4cee35d4c5a618dd4dc547f4"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.3.0"

[[deps.CoordinateTransformations]]
deps = ["LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "681ea870b918e7cff7111da58791d7f718067a19"
uuid = "150eb455-5306-5404-9cee-2592286d6298"
version = "0.6.2"

[[deps.Crayons]]
git-tree-sha1 = "3f71217b538d7aaee0b69ab47d9b7724ca8afa0d"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.0.4"

[[deps.DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

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

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "fe385ec95ac5533650fb9b1ba7869e9bc28cdd0a"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.5"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "2db648b6712831ecb333eae76dbfd1c156ca13bb"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.11.2"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.Glob]]
git-tree-sha1 = "4df9f7e06108728ebf00a0a11edee4b29a482bb2"
uuid = "c27321d9-0574-5035-807b-f59d2c89b15c"
version = "1.3.0"

[[deps.HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "0fa77022fe4b511826b39c894c90daf5fce3334a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.17"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IniFile]]
deps = ["Test"]
git-tree-sha1 = "098e4d2c533924c921f9f9847274f2ad89e018b8"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "61aa005707ea2cebf47c8d780da8dc9bc4e0c512"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.13.4"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "a7254c0acd8e62f1ac75ad24d5db43f5f19f3c65"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.2"

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
git-tree-sha1 = "642a199af8b68253517b80bd3bfd17eb4e84df6e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.3.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d735490ac75c5cb9f1b00d8b5509c11984dc6943"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.0+0"

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

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "340e257aada13f95f98ee352d316c3bed37c8ab9"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "be9eef9f9d78cecb6f262f3c10da151a6c5ab827"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.5"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NaNMath]]
git-tree-sha1 = "bfe47e760d60b82b66b61d2d44128b62e3a369fb"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.5"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "043017e0bdeff61cfbb7afeb558ab29536bbb5ed"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.10.8"

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

[[deps.OptionalData]]
git-tree-sha1 = "d047cc114023e12292533bb822b45c23cb51d310"
uuid = "fbd9d27c-2d1c-5c1c-99f2-7497d746985d"
version = "1.0.0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PROJ_jll]]
deps = ["Artifacts", "JLLWrappers", "LibCURL_jll", "LibSSH2_jll", "Libdl", "Libtiff_jll", "MbedTLS_jll", "Pkg", "SQLite_jll", "Zlib_jll", "nghttp2_jll"]
git-tree-sha1 = "2435e91710d7f97f53ef7a4872bf1f948dc8e5f8"
uuid = "58948b4f-47e0-5654-a9ad-f609743f8632"
version = "700.202.100+0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "ae4bbcadb2906ccc085cf52ac286dc1377dceccc"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.1.2"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlutoDevMacros]]
deps = ["MacroTools", "PlutoHooks"]
git-tree-sha1 = "7392720177703062cb2e2a0115efb77dc5dc818c"
uuid = "a0499f29-c39b-4c5c-807c-88074221b949"
version = "0.3.7"

[[deps.PlutoHooks]]
deps = ["FileWatching", "InteractiveUtils", "Markdown", "UUIDs"]
git-tree-sha1 = "f297787f7d7507dada25f6769fe3f08f6b9b8b12"
uuid = "0ff47ea0-7a50-410d-8455-4348d5de0774"
version = "0.0.3"

[[deps.PlutoTest]]
deps = ["HypertextLiteral", "InteractiveUtils", "Markdown", "Test"]
git-tree-sha1 = "92b8ae1eee37c1b8f70d3a8fb6c3f2d81809a1c5"
uuid = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
version = "0.2.0"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "e071adf21e165ea0d904b595544a8e514c8bb42c"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.19"

[[deps.PlutoUtils]]
deps = ["Chain", "Glob", "HypertextLiteral", "InteractiveUtils", "Markdown", "PlutoDevMacros", "PlutoHooks", "PlutoTest", "PlutoUI", "PrettyTables", "Reexport", "Requires", "UUIDs"]
git-tree-sha1 = "3d3856ecfea340b4ee0c77e5c3228dd1b4478ae1"
uuid = "ed5d0301-4775-4676-b788-cf71e66ff8ed"
version = "0.4.13"

[[deps.PolynomialRoots]]
git-tree-sha1 = "5f807b5345093487f733e520a1b7395ee9324825"
uuid = "3a141323-8675-5d76-9d11-e1df1406c778"
version = "1.0.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00cfd92944ca9c760982747e9a1d0d5d86ab1e5a"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.2"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "d940010be611ee9d67064fe559edbb305f8cc0eb"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.2.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

[[deps.Proj4]]
deps = ["CEnum", "CoordinateTransformations", "PROJ_jll", "StaticArrays"]
git-tree-sha1 = "5f15f1c647b563e49f655fbbfd4e2ade24bd3c64"
uuid = "9a7e659c-8ee8-5706-894e-f68f43bc57ea"
version = "0.7.6"

[[deps.Quaternions]]
deps = ["DualNumbers", "LinearAlgebra"]
git-tree-sha1 = "adf644ef95a5e26c8774890a509a55b7791a139f"
uuid = "94ee1d12-ae83-5a48-8b1c-48b8ff168ae0"
version = "0.4.2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "01d341f502250e81f6fec0afe662aa861392a3aa"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.2"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.ReferenceFrameRotations]]
deps = ["Crayons", "LinearAlgebra", "Printf", "StaticArrays"]
git-tree-sha1 = "d526371cec370888f485756a4bf8284ab531860b"
uuid = "74f56ac7-18b3-5285-802d-d4bd4f104033"
version = "1.0.1"

[[deps.RemoteFiles]]
deps = ["Dates", "FileIO", "HTTP"]
git-tree-sha1 = "54527375d877a64c55190fb762d584f927d6d7c3"
uuid = "cbe49d4c-5af1-5b60-bb70-0a60aa018e1b"
version = "0.4.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "4036a3bd08ac7e968e27c203d45f5fff15020621"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.1.3"

[[deps.Rotations]]
deps = ["LinearAlgebra", "Quaternions", "Random", "StaticArrays", "Statistics"]
git-tree-sha1 = "dbf5f991130238f10abbf4f2d255fb2837943c43"
uuid = "6038ab10-8711-5258-84ad-4b1120ba62dc"
version = "1.1.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.SQLite_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "e601efdf2ce1df61147a1d1fed83cc428af1b53a"
uuid = "76ed43ae-9a5d-5a62-8c75-30186b810ce8"
version = "3.36.1+0"

[[deps.SatelliteToolbox]]
deps = ["Crayons", "Dates", "DelimitedFiles", "Interpolations", "LinearAlgebra", "OptionalData", "Parameters", "PolynomialRoots", "PrettyTables", "Printf", "Reexport", "ReferenceFrameRotations", "RemoteFiles", "SparseArrays", "StaticArrays", "Statistics"]
git-tree-sha1 = "0a2c0f1565a51487fe58c28f528675dba1008432"
uuid = "6ac157d9-b43d-51bb-8fab-48bf53814f4a"
version = "0.9.3"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "f0bccf98e16759818ffc5d97ac3ebf87eb950150"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "1.8.1"

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
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "fed34d0e71b91734bf0a7e10eb1bb05296ddbcd0"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.6.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Unitful]]
deps = ["ConstructionBase", "Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "0992ed0c3ef66b0390e5752fe60054e5ff93b908"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.9.2"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "cc4bf3fdde8b7e3e9fa0351bdeedba1cf3b7f6e6"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.0+0"

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
# ╟─d852d113-2be1-4580-92dd-bf4082d0df11
# ╠═36f00194-59ac-4e1a-a746-f41c9057e972
# ╠═f43c934c-84c8-4c3d-b4d9-2b716753d89c
# ╠═f41cdadb-808d-4714-983a-b871151ff32f
# ╠═f5577c80-ffdd-44ae-bc05-2baed9de1234
# ╠═059edd4a-b3b7-4db2-9ecd-ca8a36021d2e
# ╠═e43a64ba-d776-42dd-97be-2be24a2769a7
# ╟─91045805-53e1-457a-b7d1-db5e6df5af19
# ╠═f41cdadb-808d-4714-983a-b871151ff1c0
# ╟─f5577c80-ffdd-44ae-bc05-2baed9de552d
# ╟─b2c827b1-2177-4b81-bdea-ea89242152ea
# ╟─3cc3b232-01e8-4064-8a2a-abe14aa6e5c0
# ╠═00d31f8c-dd75-4d8f-83b6-d8e976b040d0
# ╠═f91fbe7d-137f-4e05-a7c7-0486db54e39e
# ╟─46730818-1bb8-4c79-8b6f-f8cf0188c918
# ╠═17d1271f-713d-4a85-b6ef-43e2632b74cf
# ╟─965e7534-cc27-4657-b3cf-5a5b36be2a9c
# ╠═5113cbdb-6c07-4258-9d19-2d2a6b596fcd
# ╟─40363971-4729-435a-b3ae-515ac30634b0
# ╠═99a35555-52e7-4e45-b265-d3868da813a8
# ╠═2222735a-0a6b-43cd-81ea-24f9288ffa59
# ╠═9cdffed4-bc8e-4c3f-8d3e-196b687815f6
# ╠═4d406319-9640-47c5-915c-0e291c30bd15
# ╟─2af0a4dd-bc00-4563-91ed-7ba1caf6a0d6
# ╟─3a537233-2f08-451c-9cb4-dcd3723cd6c8
# ╟─a01e22f4-02cd-4427-8762-a3db2dd15112
# ╠═7bcb5a5d-6ee3-4cde-8d3c-97699c765fd3
# ╟─e19edbc3-c268-4294-b551-f6dd6964316a
# ╠═fcb9caa8-2ee3-469a-8bb7-d462ab4162bd
# ╟─0ac44137-9d7f-4746-868e-ae09b628f5e0
# ╠═ef2c3b39-5487-42ec-a006-20e0794ed21e
# ╟─ee3aa19f-317e-46f6-8da2-4792a84b7839
# ╟─1c9a8798-0b03-4e50-952e-e615192dbd45
# ╠═2e07bdfa-7393-4864-be2f-35b7843f6cc8
# ╠═e7a73ba7-731d-4a58-ac39-6fdebff78d7f
# ╟─de4f6c24-49ca-429d-aa96-11d055027bbb
# ╟─1200a81e-c1fe-4be5-a514-a87fafc9e5fb
# ╠═51a0cfb6-64a8-450a-b6e5-79f3de6c5725
# ╟─2b32d7e7-1519-4fe9-bfa8-d3c6a57b237f
# ╟─de6f9a8b-efc1-4666-88fe-31005efcd06e
# ╠═0bff095f-534e-4342-82c2-931f75e16c18
# ╟─97e3be69-b480-482b-a1aa-5bf2ede10cbe
# ╟─95704330-4d7b-44fd-b8c0-d1570812f619
# ╟─2e788b78-e5e0-4f60-aa8c-ad4f203c982e
# ╠═4cea8d15-9bb9-455c-b8bf-10b8d9a2d4af
# ╠═2af585a1-05d0-4b5a-9ee6-15eabb40a27c
# ╟─7ab00d88-9f0c-4ad9-a735-6ef845055823
# ╠═f634d5d0-bb61-4bd6-9b1c-df75399de739
# ╠═8b3f7041-ce2f-4d64-a135-9403eacd6385
# ╟─ea890a0d-b696-4131-87ea-202ee8199358
# ╠═a8761c85-73d2-457c-816f-7db2c83d01e9
# ╟─f5b5c788-ad21-478d-972f-5bc2d7fd2768
# ╟─4908872b-0894-454e-afed-0efdc0c3a84f
# ╠═d569a837-f315-47d4-9624-1bafe9996493
# ╟─d292f0d3-6a35-4f35-a5f6-e15e1c29f0f1
# ╟─4e0ecf07-e70f-4c3b-9af1-71c35167e7a8
# ╠═e7a82a42-9852-4ac3-8612-938004bf24de
# ╠═89eb0e56-e1d5-4497-8de2-3eed528f6358
# ╠═61ace485-dc58-42dd-a58f-1cd13e1f6444
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
