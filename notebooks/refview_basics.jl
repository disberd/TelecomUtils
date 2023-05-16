### A Pluto.jl notebook ###
# v0.19.25

#> custom_attrs = ["hide-enabled"]

using Markdown
using InteractiveUtils

# ╔═╡ 35da0c45-6d91-4ad9-811f-5d633684208e
begin
	using StaticArrays
	using LinearAlgebra
	using Unitful
	using Unitful.DefaultSymbols
	using SatelliteToolbox: Ellipsoid
	using Parameters
	import Proj4
	using DocStringExtensions
end

# ╔═╡ c3a5e97e-00ed-44af-9c67-fa8198900fbd
# ╠═╡ skip_as_script = true
#=╠═╡
begin
	using PlutoDevMacros
	using PlutoExtras
	using PlutoTest
	using BenchmarkTools
end
  ╠═╡ =#

# ╔═╡ e2f50ace-ad14-4e7c-af74-abf3ca9df9fb
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
# Packages
"""
  ╠═╡ =#

# ╔═╡ 82c0d6f9-373d-4866-81eb-9cd2d0981310
# ╠═╡ skip_as_script = true
#=╠═╡
ExtendedTableOfContents()
  ╠═╡ =#

# ╔═╡ 2e6bcf9b-002f-4fbc-80c1-0cfd2ab1d253
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
# Exports
"""
  ╠═╡ =#

# ╔═╡ 1e3da0c9-2f96-4637-9093-ac7f10c1ad27
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
# Helper Functions
"""
  ╠═╡ =#

# ╔═╡ 48c73104-fe4c-4543-942a-6f23b0fd2547
const Point{N} = Union{Tuple{Vararg{<:Number, N}},StaticVector{N,<:Number}}

# ╔═╡ 48c73104-fe4c-4543-0001-6f23b0fd2547
const Point2D = Point{2}

# ╔═╡ 48c73104-fe4c-4543-0002-6f23b0fd2547
const Point3D = Point{3}

# ╔═╡ c4402f72-67ac-4630-a651-da81c1df71bf
"""
	ExtraOutput
Struct used inside `SatView` to use multiple dispatch to give more than the standard output. 
See for example [`get_range`](@ref)
"""
struct ExtraOutput end

# ╔═╡ e796d5be-e011-45d3-ad90-58769feb5e85
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Show/Print
"""
  ╠═╡ =#

# ╔═╡ d890aff1-dbd0-451c-bf14-bde9758c3be0
function _print_angle(io,val,displayname,last=false)
	print(io,"$displayname=")
	print(io,round(val;digits=2) * rad)
	print(io," (")
	print(io,round(rad2deg(val);digits=2) * °)
	print(io,")")
	last || print(io,", ")
end

# ╔═╡ fdbbc8d9-83d6-413e-aff8-aca18f24dfea
function _print_length(io,val,displayname,last=false)
	print(io,"$displayname=")
	mval = val < 1000 ? round(val;digits=2) * m : round(val/1000;digits=2) * km
	print(io,mval)
	last || print(io,", ")
end

# ╔═╡ cc0cae75-ba10-4a62-b0ef-22259e40a083
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Spherical Ellipsoid
"""
  ╠═╡ =#

# ╔═╡ 855553e3-491b-4e8a-a482-95855697e063
"""
	SphericalEllipsoid(r = 6371e3)
Define a spherical ellipsoid of radius `r` to be used for the various transformation in [`SatView`](@ref). Defaults to 6371km radius
"""
SphericalEllipsoid(r = 6371e3) = Ellipsoid(r,0.0)

# ╔═╡ 1f7a7093-ce8e-461f-8b91-69266de86748
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## geod_geodesic
"""
  ╠═╡ =#

# ╔═╡ f83fc65f-5f7b-498d-9ed3-0762565ad710
Proj4.geod_geodesic(e::Ellipsoid) = Proj4.geod_geodesic(e.a, e.f)

# ╔═╡ 4714c6ae-27d9-47db-b12e-126283b10606
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
# EarthModel
"""
  ╠═╡ =#

# ╔═╡ f9cc2880-bab1-4be3-b67b-d43508df8d3b
begin
"""
$(TYPEDEF)
Geometrical model of the Earth ellipsoid to be used for all the view angles between satellite and points on earth.
# Fields
$(TYPEDFIELDS)

A single instance of this structure should be used for all satellites of a given simulation/constellation.\\
Changes to any of the two fields (via `setproperty!`) will trigger an automatic recomputation of the other field.\\
When called without arguments, it defaults to a spherical earth with a radius of 6371 km.

See also: [`Ellipsoid`](@ref), [`SphericalEllipsoid`](@ref), [`SatView`](@ref)
"""
@with_kw mutable struct EarthModel
	"Ellipsoid structure, used for the various point of view conversion."
	ellipsoid::Ellipsoid{Float64} = SphericalEllipsoid()
	"Extended geod structure, used for the inverse geodesic problem to compute distance and azimuth between points on earth. (Relies on GeographicLib)." 
	geod::Proj4.geod_geodesic = Proj4.geod_geodesic(ellipsoid)
end
EarthModel(e::Ellipsoid) = EarthModel(;ellipsoid = e)
end

# ╔═╡ 1f2b9b74-de46-401d-8a46-0434b9f9aca1
function Base.setproperty!(value::EarthModel, name::Symbol, x)
	if name === :ellipsoid
		setfield!(value, name, x)
		setfield!(value, :geod, Proj4.geod_geodesic(x))
	else
		setfield!(value, name, x)
		setfield!(value, :ellipsoid, Ellipsoid(x.a, x.f))
	end
end

# ╔═╡ 6a5cb372-60cb-4ffc-b4f0-22e4016104e7
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
# Angle Types
"""
  ╠═╡ =#

# ╔═╡ e832c1b7-8c04-4146-90f6-1628e91fea2a
const UnitfulAngleType = Union{typeof(°),typeof(rad)}

# ╔═╡ 64cf1b8b-6686-4eeb-a3cc-786300ea7c7d
const UnitfulAngleQuantity = Quantity{<:Real,<:Any,<:UnitfulAngleType}

# ╔═╡ ab2a30e0-f87b-4e65-8e56-78c283f3eb98
function _check_radians_func(limit = π) 
	f(x::Real) = abs(x) <= limit
	f(x::UnitfulAngleQuantity) = true
end

# ╔═╡ a166e0d3-bb39-4183-8669-5cb7747007d3
_check_radians(x; limit = π) = @assert all(_check_radians_func(limit), x) "Angles directly provided as numbers must be expressed in radians and satisfy -$limit ≤ x ≤ $limit
Consider using `°` from Unitful (Also re-exported by TelecomUtils) if you want to pass numbers in degrees, by doing `x * °`." 

# ╔═╡ 4c3585d9-5671-4658-a8e5-4900daf51aa4
#=╠═╡
let
	x = rand(1000)
	@benchmark _check_radians($x)
end
  ╠═╡ =#

# ╔═╡ 8f31a1f3-fe78-4aad-ae1e-91d08f85960e
#=╠═╡
let
	x = rand(1000) .* 10°
	@benchmark _check_radians(Ref($x)[])
end
  ╠═╡ =#

# ╔═╡ 1eb9b8de-fb9e-4d46-8217-78346eb2f44b
md"""
## Default Units
"""

# ╔═╡ c7ff7ef2-0e7d-4d36-8463-9c046fd36999
const ValidAngle = Union{UnitfulAngleQuantity, Real}

# ╔═╡ 34b02d15-d9d0-4b34-a1ae-bb2183b8ef47
const ValidDistance = Union{Unitful.Length, Real}

# ╔═╡ c592d893-b693-4ba3-83de-49effcf2159f
begin
function to_radians(x::Real)
	_check_radians(x)
	x
end
to_radians(x::UnitfulAngleQuantity) = uconvert(u"rad", x) |> ustrip
function to_radians(x)
	T = eltype(x)
	if T <: ValidAngle
		return map(to_radians, x)
	else
		error("You can only call `to_radians` with scalar angle values or iterables containing angle values")
	end
end
end

# ╔═╡ ed3cfb95-cacb-4cfa-af25-b03c22f1146f
begin
to_degrees(x::Real) = x
to_degrees(x::UnitfulAngleQuantity) = uconvert(u"°", x) |> ustrip
function to_degrees(x)
	T = eltype(x)
	if T <: ValidAngle
		return map(to_degrees, x)
	else
		error("You can only call `to_degrees` with scalar angle values or iterables containing angle values")
	end
end
end

# ╔═╡ 9f437287-27a2-4e7d-b394-1a4eb2ee2825
begin
function to_meters(x::Real)
	x
end
to_meters(x::Unitful.Length) = uconvert(u"m", x) |> ustrip
function to_meters(x)
	T = eltype(x)
	if T <: ValidDistance
		return map(to_meters, x)
	else
		error("You can only call `to_meters` with scalar length or iterables containing angle values")
	end
end
end

# ╔═╡ 1d023a0c-a93a-451c-a894-1d1f6a4b78a9
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
# SatViewCoordinate types
"""
  ╠═╡ =#

# ╔═╡ b8ce87d4-4768-4e8a-a16e-9e68b00b6617
abstract type SatViewCoordinate end

# ╔═╡ e3c221c6-6c4a-4b5f-93a6-30d3508ac9d2
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## LLA
"""
  ╠═╡ =#

# ╔═╡ 8368ae01-ce53-449e-87dd-8dfa3f29f8f4
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
Here we want to define a structure that contains useful informations and functions to perform conversions between the view from the satellite based on it's orbital position and points on ground
"""
  ╠═╡ =#

# ╔═╡ f951805e-515a-475f-893f-bb8b968e425c
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## ERA
"""
  ╠═╡ =#

# ╔═╡ 86ae20a9-e69c-4d63-9119-395449e9ac09
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
ERA stands for Elevation, Range and Azimuth and is used to express the position of a satellite relative to an observer in spherical coordinates.
The elevation is the angle of the pointing with respect to the local horizon of the observer, meaning the plane where the observer is located that is perpendicular to the gravity vector acting on the observe (or in an alternative definition, the plane where the observer is located that is parallel to the tangent plane to the earth ellipsoid at the given lat and lon positions of the observer.
"""
  ╠═╡ =#

# ╔═╡ 9be2fd5c-4b6c-4e13-b2aa-fb7120a504b7
begin
	"""
	Elevation, Range and Azimuth for a target point on space as seen from a source point on or above the earth surface
	
	# Fields
	- `el::Float64`: Elevation view angle (`0 <= el <= π/2`) between source and target point [rad].
	- `r::Float64`: Range (`r >= 0`) between the source and target points [m].
	- `az::Float64`: Azimuth view angle between source and target point [rad], computed from West to North from the source point perspective. Values provided are automatically converted between -π and π
	
	# Constructors
	
		ERA(el::ValidAngle,r::ValidDistance,az::ValidAngle)

	`ValidAngle` is a either a Real number or a `Unitful.Quantity` of unit either `u"rad"` or `u"°"`.

	`ValidDistance` is either a Real Number, or a subtype of `Unitful.Length`.
	"""
	@with_kw_noshow struct ERA <: SatViewCoordinate
		el::Float64 # Elevation in radians
		r::Float64 # Range in meters
		az::Float64 # Azimuth in radians
		
		function ERA(el::ValidAngle,r::ValidDistance,az::ValidAngle)
			(isnan(el) || isnan(r) || isnan(az)) && return new(NaN,NaN,NaN)  
			@assert el >= 0 && el <= π/2 "Elevation (and azimuth) should be provided in radians and between 0 and π/2.  If you wish to use degrees as inputs please multiply the first and last arguments by `°` (from Unitful.jl), which is re-exported by TelecomUtils and can be written with `\\degree<TAB>`"
			r_m = to_meters(r)
			@assert r_m >= 0 "Range must be positive"
			new(to_radians(el),r_m,rem2pi(to_radians(az),RoundNearest))
		end
	end
	
	# Show
	function Base.show(io::IO,era::ERA)
		print(io,"ERA(")
		_print_angle(io,era.el,"el",false)
		_print_length(io,era.r,"r",false)
		_print_angle(io,era.az,"az",true)
		print(io,")")
	end

function Base.isapprox(x1::ERA, x2::ERA; angle_atol = deg2rad(1e-5), range_atol = 1e-3, atol = nothing, kwargs...)
	@assert atol isa Nothing "You can't provide an absolute tolerance directly for comparison between `ERA` types, please use the independent kwargs `angle_atol` [radians] for the elevation/azimuth atol and `range_atol` [m] for the range one"

	# Range, we use 1mm as default tolerance
	Base.isapprox(x1.r,x2.r; atol = range_atol, kwargs...) || return false
	
	# Angles, we default to a default tolerance of 1e-5 degrees for isapprox
	≈(x,y) = isapprox(x,y;atol = angle_atol, kwargs...)
	≉(x,y) = !≈(x,y)
	x1.el ≈ x2.el || return false
	
	# Don't care about different azimuth if elevation is 90°
	x2.el ≉ π/2 && x1.az ≉ x2.az && return false
	return true
end
end

# ╔═╡ 39af2a26-e882-4596-846c-4699c1a0f3a2
#=╠═╡
@test ERA(0°, 0km, 0°) ≉ ERA(1.1e-5°, 0km, 0°)
  ╠═╡ =#

# ╔═╡ 4dc4e619-9d2e-45ea-9072-baef5680ef28
#=╠═╡
@test ERA(0°, 0km, 0°) ≈ ERA(1e-5°, 0km, 0°)
  ╠═╡ =#

# ╔═╡ 948e5814-60ca-4269-ae33-542b835b4116
#=╠═╡
@test ERA(0°, 0km, 0°) ≈ ERA(1e-5°, 0km + 1e-6km, 0°)
  ╠═╡ =#

# ╔═╡ 31d72681-5c68-4708-a493-be1c5f085847
#=╠═╡
@test ERA(0°, 0km, 0°) ≉ ERA(1e-5°, 0km + 1.1e-6km, 0°)
  ╠═╡ =#

# ╔═╡ 7344190c-7989-4b55-b7be-357f7d6b7370
Base.isnan(era::ERA) = isnan(era.el)

# ╔═╡ f207d849-ebff-4e6c-95bb-50693cb7c9b6
begin
	"""
	Identify a point on or above earth using geodetic coordinates
	
	# Fields
	- `lat::Float64`: Latitude (`-π/2 <= lat <= π/2`) of the point [rad].
	- `lon::Float64`: Longitude of the point [rad].
	- `alt::Float64`: Altitude of the point above the reference earth ellipsoid [m].
	
	# Constructors
		LLA(lat::ValidAngle,lon::ValidAngle,alt::ValidDistance)
		LLA(lat::ValidAngle,lon::ValidAngle) # Defaults to 0m altitude
	
	`ValidAngle` is a either a Real number or a `Unitful.Quantity` of unit either `u"rad"` or `u"°"`.

	`ValidDistance` is either a Real Number, or a subtype of `Unitful.Length`.
	"""
	@with_kw_noshow struct LLA <: SatViewCoordinate
		lat::Float64 # Latitude in radians
		lon::Float64 # Longitude in radians
		alt::Float64 # Altitude in meters
		
		function LLA(lat::ValidAngle,lon::ValidAngle,alt::ValidDistance)
			(isnan(lat) || isnan(lon) || isnan(alt)) && return new(NaN,NaN,NaN)  
			lon_rad = rem2pi(to_radians(lon),RoundNearest)
			@assert abs(lat) <= π/2 "Latitude should be given in radians and between -π/2 and π/2. If you wish to use degrees as inputs please multiply the first two arguments by `°` (from Unitful.jl), which is re-exported by TelecomUtils and can be written with `\\degree<TAB>`"
			new(to_radians(lat),lon_rad,to_meters(alt))
		end
	end
	
	# Constructor without altitude, assume it is 0
	LLA(lat,lon) = LLA(lat,lon,0.0)
end

# ╔═╡ 3033d5c1-d8e0-4d46-0001-7dec4ff7afbd
# Show method for LLA
function Base.show(io::IO,lla::LLA)
	print(io,"LLA(")
	_print_angle(io,lla.lat,"lat",false)
	_print_angle(io,lla.lon,"lon",false)
	_print_length(io,lla.alt,"alt",true)
	print(io,")")
end

# ╔═╡ 3033d5c1-d8e0-4d46-a4b9-7dec4ff7afbd
function Base.isapprox(x1::LLA, x2::LLA; angle_atol = deg2rad(1e-5), alt_atol = 1e-3, atol = nothing, kwargs...)
	@assert atol isa Nothing "You can't provide an absolute tolerance directly for comparing `LLA` objects, please use the independent kwargs `angle_atol` [radians] for the longitude and latitude atol and `alt_atol` [m] for the altitude one"
	# Altitude, we default to an absolute tolerance of 1mm for isapprox
	isapprox(x1.alt,x2.alt; atol = alt_atol, kwargs...) || return false
	# Angles, we default to a default tolerance of 1e-5 degrees for isapprox
	≈(x,y) = isapprox(x,y;atol = angle_atol, kwargs...)
	# Don't care about different longitude if latitude is ±90°
	abs(x1.lat) ≈ π/2 && abs(x2.lat) ≈ π/2 && return true
	# Return true if all the lat and lon are matching
	x1.lat ≈ x2.lat && (x1.lon ≈ x2.lon || abs(x1.lon) ≈ abs(x2.lon) ≈ π) && return true
	return false
end

# ╔═╡ cbbff280-8ad7-4589-8dff-1d401f872233
#=╠═╡
@test LLA(0°, 0°, 0km) ≉ LLA(1.1e-5°, 0°, 0km)
  ╠═╡ =#

# ╔═╡ f7ee6d27-3775-4508-8f51-c61199247e0c
#=╠═╡
@test LLA(0°, 0°, 0km) ≈ LLA(1e-5°, 0°, 0km)
  ╠═╡ =#

# ╔═╡ 1a431285-970d-456f-adfb-2e25ad405a5c
#=╠═╡
@test LLA(0°, 0°, 0km) ≈ LLA(1e-5°, 1e-5°, 1e-6km)
  ╠═╡ =#

# ╔═╡ 337271c9-7263-40e7-a444-251c03feb2f2
#=╠═╡
@test LLA(0°, 0°, 0km) ≉ LLA(1e-5°, 1e-5°, 1.1e-6km)
  ╠═╡ =#

# ╔═╡ 528beffe-0707-4661-8970-def1b1e00ea5
Base.isnan(lla::LLA) = isnan(lla.lat)

# ╔═╡ 11938cb6-46b3-0001-96c0-ef6424d1d0db
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
# Inverse Geodesic Problem
"""
  ╠═╡ =#

# ╔═╡ 11938cb6-46b3-0002-96c0-ef6424d1d0db
begin
"""
    geod_inverse(geod::Proj4.geod_geodesic, lonlat1::AbstractVector{Cdouble}, lonlat2::AbstractVector{Cdouble})
	geod_inverse(geod::Proj4.geod_geodesic, lla1::LLA, lla2::LLA)
	geod_inverse(em::EarthModel, lla1::LLA, lla2::LLA)
Solve the inverse geodesic problem.

# Args

- `g`       → the geod_geodesic object specifying the ellipsoid.
- `lonlat1` → point 1 (degrees), where lat ∈ [-90, 90], lon ∈ [-540, 540) 
- `lonlat2` → point 2 (degrees), where lat ∈ [-90, 90], lon ∈ [-540, 540) 

# Outputs

- `dist` → distance between point 1 and point 2 (meters).
- `azi1` → azimuth at point 1 (degrees) ∈ [-180, 180)
- `azi2` → (forward) azimuth at point 2 (degrees) ∈ [-180, 180)

# Remarks

If either point is at a pole, the azimuth is defined by keeping the longitude fixed,
writing lat = 90 +/- eps, and taking the limit as eps -> 0+.
"""
function geod_inverse(geod::Proj4.geod_geodesic, lonlat1::AbstractVector{Cdouble}, lonlat2::AbstractVector{Cdouble})
	dist = Ref{Cdouble}()
	azi1 = Ref{Cdouble}()
	azi2 = Ref{Cdouble}()
	ccall((:geod_inverse, Proj4.libproj), Cvoid, (Ptr{Cvoid},Cdouble,Cdouble,Cdouble,
			Cdouble,Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
			pointer_from_objref(geod), lonlat1[2], lonlat1[1], lonlat2[2], lonlat2[1], dist, azi1, azi2)
	dist[], azi1[], azi2[]
end

# Version with LLA inputs
function geod_inverse(geod::Proj4.geod_geodesic, lla1::LLA, lla2::LLA)
	lonlat1 = rad2deg.(SA_F64[lla1.lon,lla1.lat])
	lonlat2 = rad2deg.(SA_F64[lla2.lon,lla2.lat])
	geod_inverse(geod,lonlat1,lonlat2)
end

# Version with EarthModel as first input
geod_inverse(em::EarthModel, lla1::LLA, lla2::LLA) = geod_inverse(em.geod, lla1, lla2)
end

# ╔═╡ 4c06d21c-ac14-4522-bf25-2e0a1ed2d6b9
begin
	export Ellipsoid, SphericalEllipsoid
	export EarthModel
	export UnitfulAngleQuantity, UnitfulAngleType, °
	export km
	export LLA, ERA
	export geod_inverse
end

# ╔═╡ 52fafdd7-503a-4665-a86f-ddd9fd6552ea
# ╠═╡ skip_as_script = true
#=╠═╡
em = EarthModel()
  ╠═╡ =#

# ╔═╡ 87dec4f1-3842-4291-ad1c-1a384a197508
#=╠═╡
let
	lla1 = LLA(1°, 2°, 0km)
	lla2 = LLA(1°, 1°, 0km)
	@benchmark geod_inverse($em, $lla1, $lla2)
end
  ╠═╡ =#

# ╔═╡ 11938cb6-46b3-499b-96c0-ef6424d1d0db
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
# Tests
"""
  ╠═╡ =#

# ╔═╡ d2c248b1-c48e-437b-a910-edcc59b4424f
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## LLA
"""
  ╠═╡ =#

# ╔═╡ 7b306ed5-4bda-465d-abf2-4d07cb4642c1
# ╠═╡ skip_as_script = true
#=╠═╡
@test LLA(10°,10°,1000) ≈ LLA((10+100*eps())*°,10°,1000)
  ╠═╡ =#

# ╔═╡ c45f6aac-bff3-4c98-8bd5-c91a98c9eef7
# ╠═╡ skip_as_script = true
#=╠═╡
@test LLA(90°,10°,1000) ≈ LLA(90°,130°,1000)
  ╠═╡ =#

# ╔═╡ 1f1505a6-c6af-4fdd-9583-6e783e76de4f
# ╠═╡ skip_as_script = true
#=╠═╡
@test LLA(40°,-180°,1000) ≈ LLA(40°,180°,1000)
  ╠═╡ =#

# ╔═╡ 980e48bd-9dd1-4195-8086-40785a7f43e1
# ╠═╡ skip_as_script = true
#=╠═╡
@test LLA(10°,10°,1000) !== LLA((10+100*eps())*°,10°,1000)
  ╠═╡ =#

# ╔═╡ fc146750-46b2-4084-a6d2-0d91c0e104e6
# ╠═╡ skip_as_script = true
#=╠═╡
LLA(1,1,NaN) |> isnan
  ╠═╡ =#

# ╔═╡ fc3e96a6-7af8-4ef0-a4ca-b33509aa512f
md"""
## ERA
"""

# ╔═╡ ed2f2fa3-be19-4e2d-ad28-36cb053ed6bd
# ╠═╡ skip_as_script = true
#=╠═╡
@test ERA(10°,1000,20°) == ERA(10°,1km,deg2rad(20)*rad)
  ╠═╡ =#

# ╔═╡ 6151e4d1-a7a6-4fa7-bc77-0e7f3b2a3cc0
# ╠═╡ skip_as_script = true
#=╠═╡
@test ERA(90°,1000,20°) ≈ ERA(90°,1km,deg2rad(90)*rad)
  ╠═╡ =#

# ╔═╡ b5449c6a-9c43-4f9c-81f7-f01277acb109
# ╠═╡ skip_as_script = true
#=╠═╡
ERA(1,1,NaN) |> isnan
  ╠═╡ =#

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
DocStringExtensions = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Parameters = "d96e819e-fc66-5662-9728-84c9c7592b0a"
PlutoDevMacros = "a0499f29-c39b-4c5c-807c-88074221b949"
PlutoExtras = "ed5d0301-4775-4676-b788-cf71e66ff8ed"
PlutoTest = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
Proj4 = "9a7e659c-8ee8-5706-894e-f68f43bc57ea"
SatelliteToolbox = "6ac157d9-b43d-51bb-8fab-48bf53814f4a"
StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[compat]
BenchmarkTools = "~1.3.2"
DocStringExtensions = "~0.8.6"
Parameters = "~0.12.3"
PlutoDevMacros = "~0.4.5"
PlutoExtras = "~0.6.0"
PlutoTest = "~0.2.2"
Proj4 = "~0.7.6"
SatelliteToolbox = "~0.9.4"
StaticArrays = "~1.3.3"
Unitful = "~1.10.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.0"
manifest_format = "2.0"
project_hash = "6d114ac95dd2810d72264a44c51ed45a82cf50a8"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "0310e08cb19f5da31d08341c6120c047598f5b9c"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.5.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

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
git-tree-sha1 = "d9a9701b899b30332bbcb3e1679c41cce81fb0e8"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.3.2"

[[deps.CEnum]]
git-tree-sha1 = "eb4cb44a499229b3b8426dcfb5dd85333951ff90"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.2"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c6d890a52d2c4d55d326439580c3b8d0875a77d9"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.7"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "61fdd77467a5c3ad071ef8277ac6bd6af7dd4c04"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.6.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.2+0"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "fb21ddd70a051d882a1686a5a550990bbe371a95"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.4.1"

[[deps.CoordinateTransformations]]
deps = ["LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "681ea870b918e7cff7111da58791d7f718067a19"
uuid = "150eb455-5306-5404-9cee-2592286d6298"
version = "0.6.2"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "e8119c1a33d267e16108be441a287a6981ba1630"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.14.0"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "7be5f99f7d15578798f338f5433b6c432ea8037b"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.0"

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
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "b7bc05649af456efc75d178846f47006c2c4c3c7"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.13.6"

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

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "03a9b9718f5682ecb107ac9f7308991db4ce395b"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.7"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "82d7c9e310fe55aa54996e6f7f94674e2a38fcb4"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.9"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

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
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "946b56b2135c6c10bbb93efad8a78b699b6383ab"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.6"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.0"

[[deps.PlutoDevMacros]]
deps = ["HypertextLiteral", "InteractiveUtils", "MacroTools", "Markdown", "Random", "Requires"]
git-tree-sha1 = "b4b23b981704ac3e2c771a389c2899e69306c091"
uuid = "a0499f29-c39b-4c5c-807c-88074221b949"
version = "0.4.8"

[[deps.PlutoExtras]]
deps = ["AbstractPlutoDingetjes", "HypertextLiteral", "InteractiveUtils", "Markdown", "OrderedCollections", "PlutoDevMacros", "PlutoUI"]
git-tree-sha1 = "14128b3beefcf20a62f6f92be570d770b027f04f"
uuid = "ed5d0301-4775-4676-b788-cf71e66ff8ed"
version = "0.6.0"

[[deps.PlutoTest]]
deps = ["HypertextLiteral", "InteractiveUtils", "Markdown", "Test"]
git-tree-sha1 = "17aa9b81106e661cffa1c4c36c17ee1c50a86eda"
uuid = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
version = "0.2.2"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "eadad7b14cf046de6eb41f13c9275e5aa2711ab6"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.49"

[[deps.PolynomialRoots]]
git-tree-sha1 = "5f807b5345093487f733e520a1b7395ee9324825"
uuid = "3a141323-8675-5d76-9d11-e1df1406c778"
version = "1.0.0"

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

[[deps.Proj4]]
deps = ["CEnum", "CoordinateTransformations", "PROJ_jll", "StaticArrays"]
git-tree-sha1 = "5f15f1c647b563e49f655fbbfd4e2ade24bd3c64"
uuid = "9a7e659c-8ee8-5706-894e-f68f43bc57ea"
version = "0.7.6"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "dc84268fe0e3335a62e315a3a7cf2afa7178a734"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.3"

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
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SQLite_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "2c761a91fb503e94bd0130fcf4352166c3c555bc"
uuid = "76ed43ae-9a5d-5a62-8c75-30186b810ce8"
version = "3.40.0+1"

[[deps.SatelliteToolbox]]
deps = ["Crayons", "Dates", "DelimitedFiles", "Interpolations", "LinearAlgebra", "OptionalData", "Parameters", "PolynomialRoots", "PrettyTables", "Printf", "Reexport", "ReferenceFrameRotations", "RemoteFiles", "SparseArrays", "StaticArrays", "Statistics"]
git-tree-sha1 = "1831cced8785398bf38577e8cf46380d349cf4c9"
uuid = "6ac157d9-b43d-51bb-8fab-48bf53814f4a"
version = "0.9.4"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "74eaf352c0cef1e32ce7394bcc359d9199a28cf7"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.3.6"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "c79322d36826aa2f4fd8ecfa96ddb47b174ac78d"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

[[deps.URIs]]
git-tree-sha1 = "ac00576f90d8a259f2c9d823e91d1de3fd44d348"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.1"

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
git-tree-sha1 = "b95e0b8a8d1b6a6c3e0b3ca393a7a285af47c264"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.10.1"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.7.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╠═e2f50ace-ad14-4e7c-af74-abf3ca9df9fb
# ╠═35da0c45-6d91-4ad9-811f-5d633684208e
# ╠═c3a5e97e-00ed-44af-9c67-fa8198900fbd
# ╠═82c0d6f9-373d-4866-81eb-9cd2d0981310
# ╠═2e6bcf9b-002f-4fbc-80c1-0cfd2ab1d253
# ╠═4c06d21c-ac14-4522-bf25-2e0a1ed2d6b9
# ╠═1e3da0c9-2f96-4637-9093-ac7f10c1ad27
# ╠═48c73104-fe4c-4543-942a-6f23b0fd2547
# ╠═48c73104-fe4c-4543-0001-6f23b0fd2547
# ╠═48c73104-fe4c-4543-0002-6f23b0fd2547
# ╠═c4402f72-67ac-4630-a651-da81c1df71bf
# ╠═e796d5be-e011-45d3-ad90-58769feb5e85
# ╠═d890aff1-dbd0-451c-bf14-bde9758c3be0
# ╠═fdbbc8d9-83d6-413e-aff8-aca18f24dfea
# ╠═cc0cae75-ba10-4a62-b0ef-22259e40a083
# ╠═855553e3-491b-4e8a-a482-95855697e063
# ╠═1f7a7093-ce8e-461f-8b91-69266de86748
# ╠═f83fc65f-5f7b-498d-9ed3-0762565ad710
# ╠═4714c6ae-27d9-47db-b12e-126283b10606
# ╠═f9cc2880-bab1-4be3-b67b-d43508df8d3b
# ╠═1f2b9b74-de46-401d-8a46-0434b9f9aca1
# ╟─6a5cb372-60cb-4ffc-b4f0-22e4016104e7
# ╠═e832c1b7-8c04-4146-90f6-1628e91fea2a
# ╠═64cf1b8b-6686-4eeb-a3cc-786300ea7c7d
# ╠═ab2a30e0-f87b-4e65-8e56-78c283f3eb98
# ╠═a166e0d3-bb39-4183-8669-5cb7747007d3
# ╠═4c3585d9-5671-4658-a8e5-4900daf51aa4
# ╠═8f31a1f3-fe78-4aad-ae1e-91d08f85960e
# ╟─1eb9b8de-fb9e-4d46-8217-78346eb2f44b
# ╠═c7ff7ef2-0e7d-4d36-8463-9c046fd36999
# ╠═34b02d15-d9d0-4b34-a1ae-bb2183b8ef47
# ╠═c592d893-b693-4ba3-83de-49effcf2159f
# ╠═ed3cfb95-cacb-4cfa-af25-b03c22f1146f
# ╠═9f437287-27a2-4e7d-b394-1a4eb2ee2825
# ╟─1d023a0c-a93a-451c-a894-1d1f6a4b78a9
# ╠═b8ce87d4-4768-4e8a-a16e-9e68b00b6617
# ╟─e3c221c6-6c4a-4b5f-93a6-30d3508ac9d2
# ╟─8368ae01-ce53-449e-87dd-8dfa3f29f8f4
# ╠═f207d849-ebff-4e6c-95bb-50693cb7c9b6
# ╠═3033d5c1-d8e0-4d46-0001-7dec4ff7afbd
# ╠═3033d5c1-d8e0-4d46-a4b9-7dec4ff7afbd
# ╠═cbbff280-8ad7-4589-8dff-1d401f872233
# ╠═f7ee6d27-3775-4508-8f51-c61199247e0c
# ╠═1a431285-970d-456f-adfb-2e25ad405a5c
# ╠═337271c9-7263-40e7-a444-251c03feb2f2
# ╠═528beffe-0707-4661-8970-def1b1e00ea5
# ╟─f951805e-515a-475f-893f-bb8b968e425c
# ╟─86ae20a9-e69c-4d63-9119-395449e9ac09
# ╠═9be2fd5c-4b6c-4e13-b2aa-fb7120a504b7
# ╠═39af2a26-e882-4596-846c-4699c1a0f3a2
# ╠═4dc4e619-9d2e-45ea-9072-baef5680ef28
# ╠═948e5814-60ca-4269-ae33-542b835b4116
# ╠═31d72681-5c68-4708-a493-be1c5f085847
# ╠═7344190c-7989-4b55-b7be-357f7d6b7370
# ╠═11938cb6-46b3-0001-96c0-ef6424d1d0db
# ╠═11938cb6-46b3-0002-96c0-ef6424d1d0db
# ╠═52fafdd7-503a-4665-a86f-ddd9fd6552ea
# ╠═87dec4f1-3842-4291-ad1c-1a384a197508
# ╠═11938cb6-46b3-499b-96c0-ef6424d1d0db
# ╠═d2c248b1-c48e-437b-a910-edcc59b4424f
# ╠═7b306ed5-4bda-465d-abf2-4d07cb4642c1
# ╠═c45f6aac-bff3-4c98-8bd5-c91a98c9eef7
# ╠═1f1505a6-c6af-4fdd-9583-6e783e76de4f
# ╠═980e48bd-9dd1-4195-8086-40785a7f43e1
# ╠═fc146750-46b2-4084-a6d2-0d91c0e104e6
# ╠═fc3e96a6-7af8-4ef0-a4ca-b33509aa512f
# ╠═ed2f2fa3-be19-4e2d-ad28-36cb053ed6bd
# ╠═6151e4d1-a7a6-4fa7-bc77-0e7f3b2a3cc0
# ╠═b5449c6a-9c43-4f9c-81f7-f01277acb109
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
