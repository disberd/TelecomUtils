### A Pluto.jl notebook ###
# v0.19.25

#> custom_attrs = ["hide-enabled"]

using Markdown
using InteractiveUtils

# ╔═╡ 7d23f727-1907-4965-a940-cc873f6b2191
using PlutoDevMacros

# ╔═╡ 73e00cef-9734-439a-b89b-7c1d99aab74e
# ╠═╡ skip_as_script = true
#=╠═╡
begin
	using BenchmarkTools
	using PlutoTest
	using PlutoExtras
	using PlutoPlotly
	using PlutoUI
end
  ╠═╡ =#

# ╔═╡ 7729ce27-df74-4393-ab70-c4e2864c85f5
@fromparent begin
	import *
	using >.Rotations
	using >.StaticArrays
	using >.LinearAlgebra
	using >.DocStringExtensions
end

# ╔═╡ 3bda5426-c0de-493f-9514-30b6fe762463
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
# Packages
"""
  ╠═╡ =#

# ╔═╡ dcc81988-903b-4707-a70c-09c38682c80f
# ╠═╡ skip_as_script = true
#=╠═╡
ExtendedTableOfContents()
  ╠═╡ =#

# ╔═╡ 3fd1046c-fabf-4264-9638-ba41301b1804
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## load other notebook
"""
  ╠═╡ =#

# ╔═╡ a57e3983-21de-4a2e-a227-8265fee6b56b
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
# Exports
"""
  ╠═╡ =#

# ╔═╡ a5112e66-c2a2-4ed2-9951-5f97bc1745d5
md"""
# ReferenceView
"""

# ╔═╡ 9e65ad17-95bd-46fa-bc18-9d5e8c501d9a
md"""
## face_rotation
"""

# ╔═╡ f1d1295e-7fa1-44d0-bdac-8b7830da8a61
@enum Faces begin
	PositiveX = 1
	PositiveY = 2
	PositiveZ = 3
	NegativeX = -1
	NegativeY = -2
	NegativeZ = -3
end

# ╔═╡ 9e8f0786-1daa-4b9e-9172-fc0767582c7e
begin
function to_face(s::Symbol)
	s === :PositiveX && return PositiveX
	s === :PositiveY && return PositiveY
	s === :PositiveZ && return PositiveZ
	s === :NegativeX && return NegativeX
	s === :NegativeY && return NegativeY
	s === :NegativeZ && return NegativeZ
	error("Unrecognized")
end
to_face(n::Int) = Faces(n)
to_face(face::Faces) = face
end

# ╔═╡ 93222642-f2a6-4de7-8c92-21c96ef009a4
function face_rotation(face::Faces)::RotMatrix{3, Float64}
	face === PositiveZ && return one(RotMatrix{3, Float64}) # No Rotation
	face === NegativeZ && return RotY(180°)
	face === PositiveX && return RotY(+90°)
	face === NegativeX && return RotY(-90°)
	face === PositiveY && return RotX(-90°)
	face === NegativeY && return RotX(+90°)
end

# ╔═╡ cb0b070b-f70a-458e-bf72-0a4d2e93ec41
face_rotation(face) = face_rotation(to_face(face))

# ╔═╡ 28c868ca-9954-47ff-8949-a41fb7fc6d41
#=╠═╡
@benchmark Ref(face_rotation(:NegativeX))[]
  ╠═╡ =#

# ╔═╡ f4c3876a-a81b-42f3-870a-43526e4c116e
md"""
## boresight versor
"""

# ╔═╡ cc51ab70-0167-477f-a62d-88567e94fed9
begin
	function boresight_versor(face::Faces)
		face === PositiveZ && return SA_F64[0,0,1]
		face === NegativeZ && return SA_F64[0,0,-1]
		face === PositiveX && return SA_F64[1,0,0]
		face === NegativeX && return SA_F64[-1,0,0]
		face === PositiveY && return SA_F64[0,1,0]
		face === NegativeY && return SA_F64[0,-1,0]
	end
	boresight_versor(face) = boresight_versor(to_face(face))

	# convert SVectors and NTuple
	boresight_versor(s::Union{SVector{3}, NTuple{3, <:Number}}) = SVector{3, Float64}(s) |> normalize
	# Handle normal vectors
	function boresight_versor(v::AbstractVector{<:Number})
		@assert length(v) == 3 "Only Vectors with exactly three elements are supported for specifying the boresight versor. Consider using Tuples or SVector from StaticArrays"
		# Revert to SVector method
		boresight_versor(SVector{3, Float64}(v))
	end
	boresight_versor(x,y,z) = boresight_versor(SVector{3, Float64}(x,y,z))
end

# ╔═╡ c54622f5-835d-4957-bdc6-88afed9c9d2a
# ╠═╡ skip_as_script = true
#=╠═╡
boresight_versor(-3)
  ╠═╡ =#

# ╔═╡ 90811965-8bf0-4a4c-86f4-9bc5c86be71c
# ╠═╡ skip_as_script = true
#=╠═╡
boresight_versor((1,2,3))
  ╠═╡ =#

# ╔═╡ 278a536f-82dc-4f4e-8000-cc6a064ab2ee
# ╠═╡ skip_as_script = true
#=╠═╡
boresight_versor(1,2,3)
  ╠═╡ =#

# ╔═╡ 5dc71a85-a20d-4448-98b4-5065a249df1d
md"""
## Struct Definition
"""

# ╔═╡ c6ee08ba-3546-48ea-9801-edc00dfd25f0
begin
	"""
	rv = ReferenceView{T}(lla_or_ecef::Union{LLA, Point3D},earthmodel::EarthModel; kwargs...)
	
A `ReferenceView` instance `rv` tracks the instantaneous position and orientation of an object in space or on earth. It is used to provide convenience methods to assess pointing, distance and visibility between obejcts. It can also be used to extract 3D coordinates (LLA or ECEF) based on pointing angles from the `ReferenceView` instance.
	
When initializing `rv`, two convenience aliases can be used:
- `UserView === ReferenceView{:User}` can be used to represent users
- `SatView === ReferenceView{:Satellite}` can be used to represent satellites

The main difference between the two is that the default local CRS for computing pointing and visibilities is West-North-Down (+Z towards nadir) for `SatView` while it's East-North-Up (+Z towards Zenith) for `UserView`.

# Fields
$TYPEDFIELDS

When multiple `ReferenceView` objects have to be modelled/tracked, the `earthmodel` field of all of them should point to the same `EarthModel` instance, so it should be generated once and then passed to the constructor of all the ReferenceView instances to ensure that all the objects are referring to the same earth model.\\
Doing so will enforce the same Ellipsoid is shared between all satellites even when it is changed from one of the SatView instances.

To understand how the fields `rot_z`, `rot_y` and `rot_x` are used to model the
rotation applied to the default local CRS to the intended one, see the
documentation of [`crs_rotation`](@ref)

The possible values for the reference `face` (defined with respect to the axes
of the local CRS) used for pointing/visibility computations can be found looking
at [`change_reference_face!`](@ref)

See also: [`change_position!`](@ref), [`change_attitude!`](@ref),
[`change_reference_face!`](@ref), [`get_range`](@ref), [`get_pointing`](@ref),
[`get_lla`](@ref), [`get_ecef`](@ref), [`geod_inverse`](@ref),
[`get_distance_on_earth`](@ref), [`get_nadir_beam_diameter`](@ref),
[`ExtraOutput`](@ref), [`get_visibility`](@ref),
[`get_mutual_visibility`](@ref).
	"""
	Base.@kwdef mutable struct ReferenceView{T}
		"ECEF coordinates of the current satellite position"
		ecef::SVector{3,Float64}
		"LLA coordinates of the current satellite position"
		lla::LLA
		"Reference EarthModel used for the projections and for the geodesic computations"
		earthmodel::EarthModel
		"Rotation Matrix to go from the default local CRS to the ECEF CRS"
		R::RotMatrix3{Float64}
		"Rotation angle [rad] about the z axis to bring the default local CRS to the current one"
		rot_z::Float64 = 0
		"Rotation angle [rad] about the y axis to bring the default local CRS to the current one"
		rot_y::Float64 = 0
		"Rotation angle [rad] about the x axis to bring the default local CRS to the current one"
		rot_x::Float64 = 0
		"Reference race of the satellite used for the pointing computations"
		face::Faces = PositiveZ
	end
	# Version allow to specify face as Symbol
	function ReferenceView{T}(ecef, lla, earthmodel, R, rot_z, rot_y, rot_x, face::Union{Symbol, Int}) where T 
		ReferenceView{T}(ecef, lla, earthmodel, R, rot_z, rot_y, rot_x, to_face(face))
	end

	const SatView = ReferenceView{:Satellite}
	const UserView = ReferenceView{:User}

	function _refview_rotation(s::Symbol)
		s === :Satellite && return :ECEFfromWND
		s === :User && return :ECEFfromENU
		error("The provided type of ReferenceView $s is not supported")
	end

	# Custom constructor
	function ReferenceView{T}(lla_or_ecef::Union{LLA, Point3D},em::EarthModel; kwargs...) where T
		lla, ecef = if lla_or_ecef isa LLA
			lla = lla_or_ecef
			ecef = ECEFfromLLA(em.ellipsoid)(lla)
			lla, ecef
		else
			ecef = lla_or_ecef
			lla = LLAfromECEF(em.ellipsoid)(ecef)
			lla, ecef
		end
		R = _rotation_matrix(_refview_rotation(T), lla.lat, lla.lon)
		ReferenceView{T}(;ecef,lla,earthmodel=em,R, kwargs...)
	end
end

# ╔═╡ a34ca715-e7f2-4fa7-ba94-8ad3f9b1f4cd
function Base.getproperty(sv::ReferenceView, name::Symbol)
	if name ∈ (:ellipsoid, :geod)
		return getfield(getfield(sv,:earthmodel),name)
	else
		return getfield(sv, name)
	end
end

# ╔═╡ 227c0f3e-6e45-48ac-92d4-d1c1c730e4e0
function Base.setproperty!(sv::V, name::Symbol, x) where V <: ReferenceView
	error("You can't change fields of $V directly, use change_position! or change_attitude!")
end

# ╔═╡ 6eb9424e-3dd2-46d4-b4d2-81596bb81668
# ╠═╡ skip_as_script = true
#=╠═╡
begin
	em = EarthModel()
	sv = SatView(LLA(0,0,600km), em)
end
  ╠═╡ =#

# ╔═╡ 5372f3fe-699a-4f00-8e8e-36cbea224963
#=╠═╡
let
	lla = LLA(0,0,600km)
	ecef = ECEFfromLLA(em.ellipsoid)(lla)
	@benchmark SatView($lla, $em)
end
  ╠═╡ =#

# ╔═╡ 5f1fd82d-f441-4a1b-9840-773a8635d3db
#=╠═╡
@benchmark Ref(getproperty($sv, :geod))[]
  ╠═╡ =#

# ╔═╡ fb418edb-4937-4668-acf3-e22e3a818a99
md"""
## CRS / Attitude rotation
"""

# ╔═╡ 12f31d2c-1eac-47bc-bda2-a32395ff8265
"""
	R = crs_rotation(;rot_x=0.0, rot_y=0.0, rot_z=0.0)
Create the z-y'-x'' _intrinsic_ [Tait-Bryan](https://en.wikipedia.org/wiki/Euler_angles#Tait%E2%80%93Bryan_angles) rotation matrix used to transform/rotate a starting CRS to another one.
This is used within [`ReferenceView`](@ref) to account for satellite attitude or user pointing when defining the final satellite/user centric CRS.

Taking as example a `SatView`, the default CRS without any attitude rotation is always the `WND` one centered on the satellite, with the `X` axis pointing towards West, `Y` axis towards north and `Z` axis towards nadir.

Let's call this default CRS as `WND` and let's assume to have a satellite in a LEO orbit with an inclination of `70` degrees.
Let's assume that the satellite is at its ascending node, and that we want to rotate the CRS of the satellite so that the `Y` axis points in the direction of propagation. To create this new CRS that we will call `XYZ`, we have to rotate the `WND` one by `20` degrees clock-wise around the `Z` axis of the `WND` CRS.

Assuming to have a point `xyz` in coordinates relative to the `XYZ` CRS, we can express the same point as `wnd` relative to the `WND` CRS with the following equation:
```
wnd = crs_rotation(;rot_z=20°) * xyz
```
"""
crs_rotation(;rot_x = 0.0, rot_y = 0.0, rot_z = 0.0) = RotZYX(rot_z,rot_y,rot_x)

# ╔═╡ 091e4ec2-ea9e-411e-8f39-73aeb73c0214
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Change Position
"""
  ╠═╡ =#

# ╔═╡ 3c2781b8-e8df-446f-95ed-597112edabec
begin
	"""
	change_position!(rv::ReferenceView, ecef, lla::LLA, R, rot_z, rot_x, rot_y)
	change_position!(rv::ReferenceView, ecef::Point3D, lla::LLA; rot_z = 0.0, rot_y = 0.0, rot_x = 0.0)
	change_position!(rv::ReferenceView, lla::LLA; rot_z = 0.0, rot_y = 0.0, rot_x = 0.0)
	change_position!(rv::ReferenceView, ecef::Point3D; rot_z = 0.0, rot_y = 0.0, rot_x = 0.0)
	
Change the position of a [`ReferenceView`](@ref) object `rv`, also returning as output the modified
`rv`.  The function mutates the `ecef`, `lla`, `R`, `rot_z`, `rot_y`, `rot_x` fields of the `rv`
object with the values provided as arguments (when using the 1st method above).\\
For the remaining methods, the missing arguments are computed automatically.

One would normally use one of the last 3 methods so that the rotation matrix (and eventually either ECEF or LLA) are correctly computed by the function.

The first method avoids any computations but does not validate that provided inputs are
correct/consistent and refer to the same position in space. For this reason it should
only be used if those values are correctly pre-computed elsewhere and one wants to avoid the duplicate computations. 

Note: `rot_z`, `rot_y` and `rot_x` are used to obtain the actual reference CRS
to use for the pointing computations. These angles are used to transform the
default CRS (which is WND for SatView and ENU for UserView) into the desired
CRS. This is done performing a z-y-x _intrinsic_ rotation of the default CRS following the
[Tait-Bryan
Definition](https://en.wikipedia.org/wiki/Euler_angles#Tait%E2%80%93Bryan_angles).

See also: [`ReferenceView`](@ref), [`get_range`](@ref), [`get_pointing`](@ref), [`get_lla`](@ref),
[`get_ecef`](@ref), [`get_ecef`](@ref).
	"""
	function change_position!(rv::ReferenceView, ecef, lla::LLA, R, rot_z, rot_y, rot_x)
		setfield!(rv,:ecef,ecef)
		setfield!(rv,:lla,lla)
		setfield!(rv,:R,R)
		setfield!(rv,:rot_z,rot_z)
		setfield!(rv,:rot_y,rot_y)
		setfield!(rv,:rot_x,rot_x)
		return rv
	end
	
	function change_position!(rv::ReferenceView, ecef::StaticVector{3}, lla::LLA; rot_z = 0.0, rot_y = 0.0, rot_x = 0.0)
		R = _rotation_matrix(:ECEFfromUV, lla.lat, lla.lon)
		R *= crs_rotation(;rot_z, rot_y, rot_x)
		change_position!(rv, ecef, lla, R, rot_z, rot_y, rot_x)
	end
	function change_position!(rv::ReferenceView, lla::LLA; kwargs...)
		ecef = ECEFfromLLA(rv.ellipsoid)(lla)
		change_position!(rv, ecef, lla;kwargs...)
	end
	function change_position!(rv::ReferenceView, ecef::StaticVector{3}; kwargs...)
		lla = LLAfromECEF(rv.ellipsoid)(ecef)
		change_position!(rv, ecef, lla;kwargs...)
	end
end

# ╔═╡ b17a4660-3938-4816-8182-6feae471b50d
#=╠═╡
let
	rot_z = rot_x = rot_y = 0.0
	rv = SatView(LLA(0,0,600km), em)
	new_lla = LLA(0.1,0,600km)
	new_ecef = ECEFfromLLA(em.ellipsoid)(new_lla)
	_R = _rotation_matrix(:ECEFfromUV, new_lla.lat, new_lla.lon)
	R = _R * crs_rotation(;rot_z, rot_y, rot_x)
	@benchmark change_position!($rv, $new_ecef, $new_lla, $R, $rot_z, $rot_y, $rot_x)
end
  ╠═╡ =#

# ╔═╡ b4777d24-5ea7-4524-8d2a-ac5d070d58d5
#=╠═╡
let
	lla = LLA(0,0,600km)
	@benchmark _rotation_matrix(:ECEFfromUV, $lla.lat, $lla.lon)
end
  ╠═╡ =#

# ╔═╡ cddcb68d-5400-4cb9-90ef-b3e149900f8a
#=╠═╡
let
	# The benchmark below is wrong as phantom allocations are appearing from a call to _rotation_matrix (see cell above). The same does not happen outside of Pluto where this function call is just around 30ns
	rot_z = rot_x = rot_y = 0.0
	rv = SatView(LLA(0,0,600km), em)
	new_lla = LLA(0.1,0,600km)
	new_ecef = ECEFfromLLA(em.ellipsoid)(new_lla)
	@benchmark change_position!($rv, $new_ecef, $new_lla)
end
  ╠═╡ =#

# ╔═╡ 64919bbb-3b52-48e1-a7ea-c0916c47994d
#=╠═╡
let
	rv = SatView(LLA(0,0,600km), em)
	xyz = SA_F64[0,100,0] # 100m towards the direction of propagation
	wnd = crs_rotation(;rot_z = 20°) * xyz
end
  ╠═╡ =#

# ╔═╡ 5c1d0beb-c00c-40a8-b255-c887be99e706
md"""
## Change Attitude
"""

# ╔═╡ f127481d-c60e-43d0-ab9f-4d4f35984015
"""
	change_attitude!(rv::ReferenceView; rot_z = 0.0, rot_y = 0.0, rot_x = 0.0)
Modify the attitude (local CRS orientation) of `rv` by recomputing its internal rotation matrix as a product of the default location based one and the attitude rotation based on `rot_z`, `rot_y` and `rot_x`.

See [`crs_rotation`](@ref) for details on how the attitude rotation matrix is defined and computed.

See also [`ReferenceView`](ref), [`change_position!`](@ref), [`change_Reference_face!`](@ref)
"""
function change_attitude!(rv::ReferenceView; rot_z = 0.0, rot_y = 0.0, rot_x = 0.0)
	(;ecef, lla) = rv
	change_position!(rv, ecef, lla ;rot_z, rot_y, rot_x)
end

# ╔═╡ 3a745709-8f5b-4f22-848a-2f9754ab27d8
#=╠═╡
let
	sat_lla = LLA(0°, 0°, 600km)
	rv = SatView(sat_lla, em)
	wnd_R = rv.R
	sat_ecef = rv.ecef
	wnd2ecef = ECEFfromWND(rv.ecef, wnd_R, em.ellipsoid)
	# Compute the ECEF position of a point 100km below the satellite
	ref_ecef = wnd2ecef(SA_F64[0,0,100e3])
	## Now we change the attitude of the satellite
	# 90 degrees around the X axis means that the ref point should be towards the positive Y axis 
	change_attitude!(rv; rot_x = deg2rad(90))
	@test WNDfromECEF(rv.ecef, rv.R', em.ellipsoid)(ref_ecef) ≈ SA_F64[0,100e3,0]
	
	# 90 degrees around the Y axis means that the ref point should be towards -X
	change_attitude!(rv; rot_y = deg2rad(90))
	@test WNDfromECEF(rv.ecef, rv.R', em.ellipsoid)(ref_ecef) ≈ SA_F64[-100e3,0,0]

	# 90 degrees around the Z and 90 degrees around the Y axis means that the ref point should be towards -X
	change_attitude!(rv; rot_z = deg2rad(90), rot_y = deg2rad(90))
	@test WNDfromECEF(rv.ecef, rv.R', em.ellipsoid)(ref_ecef) ≈ SA_F64[-100e3,0,0]
	
	# 45 degrees around the Z and 45 degrees around the Y axis means that the ref point should be at 45 degrees between the -X and +Z
	change_attitude!(rv; rot_z = deg2rad(45), rot_y = deg2rad(45))
	@test WNDfromECEF(rv.ecef, rv.R', em.ellipsoid)(ref_ecef) ≈ SA_F64[-√(100e3^2/2),0,√(100e3^2/2)]
end
  ╠═╡ =#

# ╔═╡ 851f1994-22f4-4347-b106-2fa3ea75ebf2
md"""
## Change Reference Face
"""

# ╔═╡ 74d39d15-8d6c-4fb8-8d04-18f32c393ad4
"""
	change_reference_face!(rv::ReferenceView, face)
Modify the reference face that is used to compute the visibilities from object `rv`.
The face can be specified using an instance of `TelecomUtils.Faces` (not exported), a Symbol or an Int as follows:   
- TelecomUtils.PositiveX **or** :PositiveX **or** 1
- TelecomUtils.PositiveY **or** :PositiveY **or** 2
- TelecomUtils.PositiveZ  **or** :PositiveZ **or** 3
- TelecomUtils.NegativeX **or** :NegativeX **or** -1
- TelecomUtils.NegativeY **or** :NegativeY **or** -2
- TelecomUtils.NegativeZ **or** :NegativeZ **or** -3

See also [`ReferenceView`](ref), [`change_position!`](@ref), [`change_attitude!`](@ref)
"""
function change_reference_face!(rv::ReferenceView, face)
	setfield!(rv, :face, to_face(face))
	return rv
end

# ╔═╡ 8b68c7c2-5cbd-4edd-9739-0d2e8c7f5449
#=╠═╡
let
	rv = SatView(LLA(0,0,500km), em)
	change_reference_face!(rv, :NegativeZ)
end
  ╠═╡ =#

# ╔═╡ 84769564-8ba8-46f5-b494-b0689d9abd65
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Get Range
"""
  ╠═╡ =#

# ╔═╡ 0cde0a71-7f27-4290-88cd-2cccf627926b
begin
"""
	get_range(rv::ReferenceView,uv::Point2D[, ::ExtraOutput]; h = 0.0, face = rv.face, R = nothing)
	get_range(rv::ReferenceView,target::Union{LLA, Point3D, ReferenceView}[, ::ExtraOutput]; face = rv.face, R = nothing)

Get the range [in m] between the reference view `rv` (Satellite or User) and a target point. The target point can be identified in two ways:
- Providing the uv pointing `uv` and the reference altitude `h` [m] above the ellipsoid (first method)
- directly passing in a point expressed either as `LLA`, as a `Point3D` representing its ECEF coordinates, or as another instance of `ReferenceView` (second method).

The pointing is assumed to be referred to the specified `face` of the
`ReferenceView`. When the location identified by the provided pointing is not
visible either because it's blocked by earth or because it's in a direction not
visible from the specified `face`, `NaN` is returned.

See [`change_reference_face!`](@ref) for supported identification of the `face` kwarg.

When called with an instance of `TelecomUtils.ExtraOutput` as last argument, the
function also returns the coordinated of the identified point in the local CRS
of `rv`.

The kwarg `R` represents the 3D Rotation Matrix that translates a vector from
ECEF coordinates to the coordinates of the desired local CRS around `rv`. By
default (if `R === nothing`) this rotation matrix is computed based on the
rotation matrix of the `rv` object and on the selected reference face.

See also: [`ReferenceView`](@ref), [`get_era`](@ref), [`get_mutual_pointing`](@ref), [`get_pointing`](@ref), [`get_lla`](@ref), [`get_ecef`](@ref), [`get_distance_on_earth`](@ref).
"""
function get_range(sv::ReferenceView,uv::Point2D, eo::ExtraOutput; h = 0.0, face = sv.face, R = nothing)
	_R = isnothing(R) ? inv(sv.R * face_rotation(face)) : R 
	# Find the ecef coordinate of the target point on earth
	ecef, r = ECEFfromUV(sv.ecef, _R', sv.ellipsoid)(uv, eo;h)
	# Compute the point coordinates in the ReferenceView local CRS
	xyz = r * SVector(uv[1], uv[2], sqrt(1 - (uv[1]^2 + uv[2]^2)))
	# Return the distance betwteen the satellite and the point
	return r, xyz
end

# LLA or ECEF version
function get_range(sv::ReferenceView, lla_or_ecef::Union{LLA, Point3D}, ::ExtraOutput; face = sv.face, R = nothing)
	ecef = if lla_or_ecef isa LLA
		ECEFfromLLA(sv.ellipsoid)(lla_or_ecef)
	else
		lla_or_ecef
	end
	_R = isnothing(R) ? inv(sv.R * face_rotation(face)) : R 
	Δecef = ecef - sv.ecef
	dist = norm(Δecef)
	# Find if the target point is below the satellite, we do this by checking the last coordinate of the WND coordinates of the point
	xyz = _R * Δecef
	xyz[3] < 0 && return NaN, xyz
	# We have to check that the given lla is visible from the satellite, this happens if either there is no intersection with earth in the direction of pointing, or if the first intersection happens for a range greater than th computed one
	t₁, t₂ = _intersection_solutions(Δecef/dist, sv.ecef, sv.ellipsoid.a, sv.ellipsoid.b)
	ref_t = t₁ > 0 ? t₁ : t₂ # t₁ always < t₂, so if positive that is the shortest distance, if negative t₂ might be also negative or positive, so we choose that
	r = (
		ref_t > 0 && # Check that the intersection is in the direction of pointing
		ref_t < dist + 1e-5 # We aim for a precision of 1e-5 m. if ref_t is lower than the distance, it means that the target is behind earth
		) ? NaN : dist
	return r, xyz
end

# Single Output Version
get_range(rv, p; kwargs...) = get_range(rv, p, ExtraOutput(); kwargs...)[1]

# 2 ReferenceView version
function get_range(rv1::ReferenceView, rv2::ReferenceView, args...; kwargs...)
	@assert rv1.earthmodel === rv2.earthmodel "When computing range between ReferenceView objects, the `EarthModel` they use internally must be the same"
	get_range(rv1, rv2.ecef, args...; kwargs...)
end
end

# ╔═╡ 33e4c937-4443-4d97-a3ed-81479ede1e11
#=╠═╡
get_range(sv, SatView(LLA(0,0,700km), em); face = :NegativeZ)
  ╠═╡ =#

# ╔═╡ bd62bdd6-4de4-449c-b5f1-fb1b4f695cda
#=╠═╡
let
	lla_ref = LLA(0,0,10km)
	@benchmark get_range($sv, $lla_ref)
end
  ╠═╡ =#

# ╔═╡ 21bc362a-960c-4c5f-9118-7e1451edb996
#=╠═╡
let
	lla_ref = LLA(0,0,10km)
	ecef_ref = ECEFfromLLA(sv.ellipsoid)(lla_ref)
	@benchmark get_range($sv, $ecef_ref)
end
  ╠═╡ =#

# ╔═╡ fd3a1dcb-2e64-47f8-a932-452742090ac1
#=╠═╡
let
	sv = SatView(LLA(0,0,600km), em)
	get_range(sv, (.1,0), ExtraOutput(); h = 100e3, face = :PositiveZ)
end
  ╠═╡ =#

# ╔═╡ 3381783a-b6d7-43c7-947d-c1de4437e25b
#=╠═╡
let
	uv = UserView(LLA(0,0,0km), em)
	get_range(uv, (0,0); h = 100e3)
end
  ╠═╡ =#

# ╔═╡ 449b49de-2951-41fc-ba46-89eaa6c52e79
# ╠═╡ skip_as_script = true
#=╠═╡
get_range(sv,LLA(0°,5°,10km),ExtraOutput())
  ╠═╡ =#

# ╔═╡ f99b984a-b7cf-4318-8f74-aacb644bc146
#=╠═╡
let
	sv = SatView(LLA(0,0,600km), em)
	get_range(sv, (0,0); face = :PositiveZ)
end
  ╠═╡ =#

# ╔═╡ b389a1db-c444-481c-b168-3069c7367276
#=╠═╡
let
	sv = SatView(LLA(0,0,600km), em)
	@benchmark get_range($sv,LLA(0,0,10km))
end
  ╠═╡ =#

# ╔═╡ 39a1850b-f64a-4157-8f07-d7a78918fea1
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Get Pointing
"""
  ╠═╡ =#

# ╔═╡ 83634223-87d0-4c31-801a-af8a7f9f678a
begin
function _get_pointing(rv::ReferenceView, lla_or_ecef::Union{LLA, Point3D}, eo::ExtraOutput, pointing_type::PointingType; face = rv.face, R = nothing)
	ecef = if lla_or_ecef isa LLA
		ECEFfromLLA(rv.ellipsoid)(lla_or_ecef)
	else
		lla_or_ecef
	end
	# Check if earth is blocking
	blocked, normalized_pdiff, pdiff = block_data = earth_blocking(rv.ecef,ecef, rv.ellipsoid.a,rv.ellipsoid.b, eo)

	# If there is an intersection, we just return false
	blocked && return SA_F64[NaN, NaN], SA_F64[NaN, NaN, NaN], block_data
	_R = isnothing(R) ? inv(rv.R * face_rotation(face)) : R 
	xyz = _R * pdiff
	x,y,z = normalize(xyz)
	# If the target is behind, we return NaN as first output
	z < 0 && return SA_F64[NaN, NaN], xyz, block_data
	if pointing_type isa ThetaPhi
		return SA_F64[acos(z), atan(y,x)], xyz, block_data
	else
		return SA_F64[x, y], xyz, block_data
	end
end

"""
	get_pointing(rv::ReferenceView, target::Union{LLA, Point3D, ReferenceView}[, ::ExtraOutput]; pointing_type::Symbol=:uv, face = rv.face, R=nothing)
Provide the 2-D angular pointing at which the target point (specified as LLA, ECEF or as another `ReferenceView`) is seen from the ReferenceView object `rv`.

`pointing_type` is used to select whether the output is returned in UV or
ThetaPhi [rad] coordinates. The following symbols are supported for this kwarg:
- `:thetaphi`, `:ThetaPhi` and `:θφ` can be used to represent pointing in ThetaPhi [rad]
- `:UV` and `:uv` can be used to represent pointing in UV

When called with an instance of `TelecomUtils.ExtraOutput` as last argument, the
function also returns the coordinated of the identified point in the local CRS
of `rv` as second output, and a NamedTuple containing information on the earth
blockage status as third argument.

For details on how to modify the reference pointing direction using the kwargs `face` and `R` look at the documentation of [`get_range`](@ref) 

See also: [`get_mutual_pointing`](@ref), [`ReferenceView`](@ref),
[`get_range`](@ref), [`get_pointing`](@ref), [`get_lla`](@ref),
[`get_ecef`](@ref), [`get_distance_on_earth`](@ref), [`get_visibility`](@ref),
[`get_mutual_visibility`](@ref).
"""
function get_pointing(sv::ReferenceView, lla_or_ecef::Union{LLA, Point3D}, eo::ExtraOutput; pointing_type::Symbol=:uv, kwargs...)
	_get_pointing(sv, lla_or_ecef, eo, PointingType(pointing_type); kwargs...)
end

# Single output version
get_pointing(sv, p; kwargs...) = get_pointing(sv, p, ExtraOutput(); kwargs...)[1]

# 2 ReferenceView version
function get_pointing(rv1::ReferenceView, rv2::ReferenceView, args...; kwargs...)
	@assert rv1.earthmodel === rv2.earthmodel "When computing visibility angles between ReferenceView objects, the `EarthModel` they use internally must be the same"
	get_pointing(rv1, rv2.ecef, args...; kwargs...)
end
end

# ╔═╡ d4ff395e-55bd-4240-89e3-a4be2dcf9ebe
#=╠═╡
let
	lla_ref = LLA(1°, 0.5°, 0km)
	ecef_ref = ECEFfromLLA(em.ellipsoid)(lla_ref)
	@benchmark get_pointing($sv, $lla_ref)
	# @code_warntype get_pointing(sv, lla_ref)
end
  ╠═╡ =#

# ╔═╡ cb4d8cd7-ff3e-43e6-ad38-8db343087214
#=╠═╡
let
	lla_ref = LLA(1°, 0.5°, 0km)
	ecef_ref = ECEFfromLLA(em.ellipsoid)(lla_ref)
	@benchmark get_pointing($sv, $ecef_ref)
end
  ╠═╡ =#

# ╔═╡ f3b8b9ba-c428-4868-96b9-1983eabb3774
#=╠═╡
let
	lla_ref = LLA(1°, 0.5°, 0km)
	rv2 = SatView(lla_ref, em)
	@benchmark get_pointing($sv, $rv2)
end
  ╠═╡ =#

# ╔═╡ 951fc2e2-ada9-4aad-876a-4cbd46200f6c
#=╠═╡
@bind test_face let
	map(Symbol, collect(instances(Faces))) |> Select
end
  ╠═╡ =#

# ╔═╡ a93e2354-34d1-4d1d-ac3c-f4c99099820a
#=╠═╡
let
	sat_lla = LLA(0°,0°,600km)
	sv = SatView(sat_lla, em)
	directions = (;
		PositiveZ = LLA(0°,0°,500km),
		PositiveX = LLA(0°,-1°,600km),
		PositiveY = LLA(1°,0°,600km),
		NegativeZ = LLA(0°,0°,700km),
		NegativeX = LLA(0°,1°,600km),
		NegativeY = LLA(-1°,0°,600km),
	)
	face = test_face
	ref_uv = get_pointing(sv, getproperty(directions,face); face = face)
end
  ╠═╡ =#

# ╔═╡ 46eb3a4d-c80b-41a0-9333-aaf6411a010c
#=╠═╡
let
	sv = SatView(LLA(0,0,600km),em)
	lla = map(rand(1000), rand(1000)) do lat,lon
		LLA(lat*°,lon*°,0km)
	end
	ecef = map(lla) do lla
		ECEFfromLLA(sv.ellipsoid)(lla)
	end
	faces = rand([:PositiveZ, :PositiveY, :NegativeX], 1000)
	map((x,y) -> get_pointing(sv,x; pointing_type= :thetaphi, face = y), ecef, faces)
end
  ╠═╡ =#

# ╔═╡ 314c4c2a-6f3b-4a7d-aa83-9073eb67897b
#=╠═╡
let
	rv1 = SatView(LLA(0°, 10°, 800km), em)
	lat = range(0°, 90°; step = 0.1°)
	θ = map(lat) do lat
		rv2 = SatView(LLA(lat, 10°, 600km), em; face = -3)
		# v = get_visibility(rv1, rv2; fov = 60°)
		p = get_pointing(rv1, rv2; pointing_type = :thetaphi)[1] |> rad2deg
	end
	x = to_degrees(lat)
	plot(scatter(;x, y = θ))
	# v,p
end
  ╠═╡ =#

# ╔═╡ a61e66cf-b430-4d84-a753-faf42f4b6337
md"""
## Get Mutual Pointing
"""

# ╔═╡ 98f3f83d-fcde-48ba-8855-c30643776d81
begin
"""
	p₁, p₂ = get_mutual_pointing(rv1::ReferenceView, rv2::ReferenceView[, ::ExtraOutput]; pointing_type::Symbol=:uv, faces = (rv1.face, rv2.face), Rs=(nothing, nothing))

Provide the 2-D angular pointing in both directions between `rv1` and `rv2`:
- `p₁` is the pointing of `rv2` with respect to `rv1`
- `p₂` is the pointing of `rv1` with respect to `rv2`

`pointing_type` is used to select whether the outputs are returned in UV or
ThetaPhi [rad] coordinates. The following symbols are supported for this kwarg:
- `:thetaphi`, `:ThetaPhi` and `:θφ` can be used to represent pointing in ThetaPhi [rad]
- `:UV` and `:uv` can be used to represent pointing in UV

When called with an instance of `TelecomUtils.ExtraOutput` as last argument, the
function also returns the coordinated of the identified point in the local CRS
of `rv1` (or `rv2`). In this case:
- `p₁` is a tuple containing the pointing as first argument, the local CRS coordinates of `rv2` with respect to `rv1` as second argument and a NamedTuple with earth blockage data as third argument.
- `p₂` is a tuple containing the pointing as first argument, the local CRS coordinates of `rv1` with respect to `rv2` as second argument and a NamedTuple with earth blockage data as third argument.

`faces` and `Rs` are tuples containing the values of `face` and `R` for the two
ReferenceView objects. `faces[1]` is used as reference face for `rv1` while
`faces[2]` is used for `rv2`. Similarly for `Rs`.

For details on how to modify the reference pointing direction using `face` and
`R` look at the documentation of [`get_range`](@ref) 

See also: [`ReferenceView`](@ref), [`get_range`](@ref), [`get_pointing`](@ref),
[`get_lla`](@ref), [`get_ecef`](@ref), [`get_distance_on_earth`](@ref),
[`get_visibility`](@ref), [`get_mutual_visibility`](@ref).
"""
function get_mutual_pointing(rv1::ReferenceView, rv2::ReferenceView, eo::ExtraOutput; pointing_type = :uv, faces = (rv1.face, rv2.face), Rs = (nothing, nothing))
	# Pointing of rv2 as seen from rv1
	p₁, xyz₁, block_data = data1 = get_pointing(rv1, rv2, eo; pointing_type, face = faces[1], R = Rs[1])

	block_data.blocked && return data1, (SA_F64[NaN, NaN], SA_F64[NaN, NaN, NaN], block_data)

	block_data = (;block_data..., pdiff = -block_data.pdiff, normalized_pdiff = -block_data.normalized_pdiff)

	R = Rs[2]
	_R = isnothing(R) ? inv(rv2.R * face_rotation(faces[2])) : R 
	xyz = _R * block_data.pdiff
	x,y,z = normalize(xyz)
	# If the target is behind, we return NaN as first output
	z < 0 && return data1, (SA_F64[NaN, NaN], xyz, block_data)
	data2 = if pointing_type isa ThetaPhi
		SA_F64[acos(z), atan(y,x)], xyz, block_data
	else
		SA_F64[x, y], xyz, block_data
	end
	
	return data1, data2
end
get_mutual_pointing(rv1, rv2; kwargs...) = first.(get_mutual_pointing(rv1, rv2, ExtraOutput(); kwargs...))
end

# ╔═╡ cda9503a-5262-4318-9f5d-aea7910ab42e
#=╠═╡
let
	sv1 = sv
	sv2 = SatView(LLA(0°,1°, 1000km), em)
	@benchmark get_mutual_pointing($sv1, $sv2; faces = (-3,3))
end
  ╠═╡ =#

# ╔═╡ 868bb972-aa0b-484d-b906-1628412b1646
#=╠═╡
let
	sv1 = sv
	sv2 = SatView(LLA(0°,1°, 1000km), em)
	@benchmark get_mutual_pointing($sv1, $sv2, ExtraOutput();pointing_type = :thetaphi, faces = (-3,3))
end
  ╠═╡ =#

# ╔═╡ 83098fe4-2089-44be-8db9-6e2db74d886d
#=╠═╡
let
	sv1 = sv
	sv2 = SatView(LLA(0°,1°, 1000km), em)
	get_mutual_pointing(sv1, sv2, ExtraOutput();pointing_type = :thetaphi, faces = (-3,3))
end
  ╠═╡ =#

# ╔═╡ 1de548c0-3830-4161-b80c-be72dd894e85
md"""
## Get Visibility
"""

# ╔═╡ 664be7e6-28f5-4ab3-b8f0-b7618bcd5fe5
# ╠═╡ skip_as_script = true
#=╠═╡
begin
@addmethod function earth_blocking(ecef1, ecef2, a, b, ::ExtraOutput)
	@inline
	# Check if the given ecef coordinate is visible from the satellite position or is obstructed from earth
	pdiff = (ecef2 - ecef1)
		
	# Find the magnitude of the difference to compare with the intersection solutions
	t = norm(pdiff)
	normalized_pdiff = pdiff ./ t
	t₁,t₂ = _intersection_solutions(normalized_pdiff,ecef1,a,b)
	# If there is an intersection, we return true and the normalized_pdiff
	blocked = t₁ > 0 && t > t₁+1e-3
	return (;blocked, normalized_pdiff, pdiff)
end

# Single Output Version
@addmethod earth_blocking(ecef1, ecef2, a, b) = earth_blocking(ecef1, ecef2, a, b, ExtraOutput())[1]
end
  ╠═╡ =#

# ╔═╡ 387d2c76-1a08-441c-97fc-7b0a90a95c9a
begin
"""
	get_visibitiliy(rv::ReferenceView, target::Union{LLA, Point3D, ReferenceView}[, ::ExtraOutput]; boresight = rv.face, fov = 90°)
Returns `true` if `target` is visible from `rv` assuming an antenna boresight
direction specified by `boresight` and a maximum Field of View from the
boresight specified by `fov`.

The `boresight` kwarg can be specified either as a `face` (compatible with the
input types specified in [`change_reference_face!`](@ref)) or as a 3D pointing
vector using either 3 elements Tuple, Vector or SVector. If a 3D vector is
provided, the function normalizes to unitary norm automatically.

The `fov` kwarg must be specified as an angle in radians if provided as a simple
Number. To provide a fov in degress, directly use ° from Unitful re-exported by
TelecomUtils.

When called with an instance of `TelecomUtils.ExtraOutput` as last argument, the
function also returns the pointing angle θ between `rv` and `target` as second
argument, and the direction in ECEF coordinates from `rv` to `target`.

See also: [`get_pointing`](@ref), [`get_mutual_pointing`](@ref),
[`ReferenceView`](@ref), [`get_range`](@ref), [`get_pointing`](@ref),
[`get_lla`](@ref), [`get_ecef`](@ref), [`get_distance_on_earth`](@ref),
[`get_mutual_visibility`](@ref).
"""
function get_visibility(rv::ReferenceView, lla_or_ecef::Union{LLA, Point3D}, eo::ExtraOutput; boresight = rv.face, fov = 90°)
	ecef = if lla_or_ecef isa LLA
		ECEFfromLLA(rv.ellipsoid)(lla_or_ecef)
	else
		lla_or_ecef
	end	
	blocked, normalized_pdiff, _ = block_data = earth_blocking(rv.ecef,ecef, rv.ellipsoid.a,rv.ellipsoid.b, eo)

	# # If there is an intersection, we just return false
	blocked && return false, NaN, normalized_pdiff

	xyz = rv.R' * normalized_pdiff
	θ = acos(dot(xyz, boresight_versor(boresight)))
	return θ <= fov, θ, normalized_pdiff
end
	
	# Double RefView method
	get_visibility(rv1::ReferenceView, rv2::ReferenceView, args...; kwargs...) = get_visibility(rv1, rv2.ecef, args...; kwargs...)

	# Single Output
	get_visibility(rv1::ReferenceView, target; kwargs...) = get_visibility(rv1, target, ExtraOutput(); kwargs...)[1]
end

# ╔═╡ abb3a7ee-ff60-421c-8a85-23c3ef1ce6af
#=╠═╡
let
	rv1 = SatView(LLA(0°, 10°, 800km), em)
	rv2 = SatView(LLA(35°, 10°, 600km), em; face = -3)
	v = get_visibility(rv1, rv2; fov = 60°)
	p = get_pointing(rv1, rv2, ExtraOutput(); pointing_type = :uv)
	v,p
end
  ╠═╡ =#

# ╔═╡ b75df549-ed46-4def-9fe0-e63ac496b3f4
#=╠═╡
let
	rv1 = SatView(LLA(0°, 10°, 800km), em)
	rv2 = LLA(0°, 190°, 600km)
	v = get_visibility(rv1, rv2; fov = 60°)
	p = get_pointing(rv1, rv2, ExtraOutput(); pointing_type = :uv)
	v,p
end
  ╠═╡ =#

# ╔═╡ 864cc72a-6252-45a6-8ed4-81d960c0b080
#=╠═╡
let
	rv1 = SatView(LLA(1°, 10°, 800km), em)
	rv2 = SatView(LLA(0°, 10°, 600km), em; face = -3)
	v = get_visibility(rv1, rv2; fov = 60°)
	p = get_pointing(rv1, rv2; pointing_type = :thetaphi) |> x -> rad2deg.(x)
	v,p
end
  ╠═╡ =#

# ╔═╡ ecd8d6ac-27f6-4a50-81ff-9c8f0b7fed47
#=╠═╡
let
	rv1 = SatView(LLA(50°, 10°, 800km), em)
	rv2 = SatView(LLA(0°, 10°, 600km), em; face = -3)
	@benchmark get_visibility($rv1, $rv2; fov = 60°)
end
  ╠═╡ =#

# ╔═╡ fe0680dd-95fb-40a0-8063-97c009f552dd
#=╠═╡
let
	rv1 = SatView(LLA(50°, 10°, 800km), em)
	rv2 = SatView(LLA(0°, 10°, 600km), em; face = -3)
	@benchmark get_pointing($rv1, $rv2; pointing_type = :thetaphi)
end
  ╠═╡ =#

# ╔═╡ f0565008-d355-4ba4-9269-38a8dc4d88bf
md"""
## Get Mutual Visibility
"""

# ╔═╡ abd76e81-2464-4119-a9bf-6046805f8e57
begin
"""
	get_mutual_visibitiliy(rv1::ReferenceView, rv2::ReferenceView[, ::ExtraOutput]; boresights = (rv1.face, rv2.face), fov = (90°, 90°))
Similar to [`get_visibility`](@ref), returns `true` if `rv1` and `rv2` can see
each other, assuming their antenna boresight directions to be specified by
`boresights` and their maximum Field of View from the boresight to be specified
by `fovs`.

`boresights` and `fovs` are tuples containing the values of boresight and fov for the two
ReferenceView objects. `boresights[1]` is used as reference boresight for `rv1` while
`boresights[2]` is used for `rv2`. Similarly for `fovs`.

Each `boresight` value inside the `boresights` kwarg can be specified either as a `face` (compatible with the
input types specified in [`change_reference_face!`](@ref)) or as a 3D pointing
vector using either 3 elements Tuple, Vector or SVector. If a 3D vector is
provided, the function normalizes to unitary norm automatically.

Each `fov` value inside the `fovs` kwarg must be specified as an angle in
radians if provided as a simple Number. To provide a fov in degress, directly
use ° from Unitful re-exported by TelecomUtils.

When called with an instance of `TelecomUtils.ExtraOutput` as last argument, the
function also returns a Tuple containing the full outputs (as if called with
ExtraOutput) of [`get_visibility`](@ref) for each direction (`rv1` to `rv2` and
vice-versa).

See also: [`get_pointing`](@ref), [`get_mutual_pointing`](@ref),
[`ReferenceView`](@ref), [`get_range`](@ref), [`get_pointing`](@ref),
[`get_lla`](@ref), [`get_ecef`](@ref), [`get_distance_on_earth`](@ref),
[`get_visibility`](@ref).
"""
function get_mutual_visibility(rv1::ReferenceView, rv2::ReferenceView, eo::ExtraOutput; boresights = (rv1.face, rv2.face), fovs = (90°, 90°))
	fwd = get_visibility(rv1, rv2, eo; boresight = boresights[1], fov = fovs[1])

	# We compute the forward visibility, which also checks for earth intersection
	fwd[1] || return false, (fwd, fwd)

	normalized_pdiff = -fwd[3]

	xyz = rv2.R' * normalized_pdiff
	θ = acos(dot(xyz, boresight_versor(boresights[2])))
	rtn = θ <= fovs[2], θ, normalized_pdiff
	return fwd[1] && rtn[1], (fwd, rtn)
end

	# Single Output
	get_mutual_visibility(rv1::ReferenceView, target; kwargs...) = get_mutual_visibility(rv1, target, ExtraOutput(); kwargs...)[1]
end

# ╔═╡ de0372e1-da0c-4da9-8236-790ecd667d07
#=╠═╡
let
	rv1 = SatView(LLA(0°, 10°, 800km), em)
	rv2 = SatView(LLA(1°, 10°, 600km), em; face = -3)
	v = get_mutual_visibility(rv1, rv2, ExtraOutput(); fovs = (60°, 60°))
end
  ╠═╡ =#

# ╔═╡ a9051f97-2561-4317-b215-eac2b4f20aa4
#=╠═╡
let
	rv1 = SatView(LLA(0°, 10°, 800km), em)
	rv2 = SatView(LLA(1°, 10°, 600km), em; face = -3)
	@benchmark get_mutual_visibility($rv1, $rv2; fovs = (60°, 60°))
end
  ╠═╡ =#

# ╔═╡ 1f27b72f-9a3b-4732-a98e-d216af067072
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Get ECEF
"""
  ╠═╡ =#

# ╔═╡ 12330b6f-97b0-4efb-9885-49758bc2f127
begin
function _get_ecef(rv::ReferenceView, pointing::Point2D, eo::ExtraOutput, pointing_type::PointingType; h = 0.0, face = rv.face, R = nothing)
	_R = isnothing(R) ? inv(rv.R * face_rotation(face)) : R 
	uv = if pointing_type isa ThetaPhi
		UVfromThetaPhi()(pointing)
	else
		pointing
	end
	ecef, r = ECEFfromUV(rv.ecef,_R',rv.ellipsoid)(uv, eo;h)
	xyz = XYZfromUV()(uv, r)
	return ecef, xyz
end
"""
	get_ecef(rv::ReferenceView, pointing::Point2D[, ::ExtraOutput]; pointing_type::Symbol=:uv, h = 0.0, face = rv.face, R = nothing)

Computes the ECEF coordinates of the point that is seen by `rv` in the direction
specified by `pointing` and is located at a target altitude `h` [m] above the
earth's surface.

If a valid point can not be found because either earth is blocking the view or
no point at altitude `h` can be seen from the provided pointing direction in the
`rv` local CRS (also accounting for desired face), the function returns a
SVector{3, Float64} filled win NaNs.

`pointing_type` is used to select whether the output is returned in UV or
ThetaPhi [rad] coordinates. The following symbols are supported for this kwarg:
- `:thetaphi`, `:ThetaPhi` and `:θφ` can be used to represent pointing in ThetaPhi [rad]
- `:UV` and `:uv` can be used to represent pointing in UV

When called with an instance of `TelecomUtils.ExtraOutput` as last argument, the
function also returns the coordinated of the identified point in the local CRS
of `rv`.

For details on how to modify the reference pointing direction using the kwargs
`face` and `R` look at the documentation of [`get_range`](@ref) 

See also: [`ReferenceView`](@ref), [`get_range`](@ref), [`get_pointing`](@ref),
[`get_lla`](@ref), [`get_era`](@ref), [`get_distance_on_earth`](@ref),
[`get_visibility`](@ref), [`get_mutual_visibility`](@ref).
"""
function get_ecef(rv::ReferenceView, pointing::Point2D, eo::ExtraOutput; pointing_type::Symbol=:uv, kwargs...)
	_get_ecef(rv, pointing, eo, PointingType(pointing_type); kwargs...)
end

# Single Output Version
get_ecef(rv::ReferenceView, pointing::Point2D;kwargs...) = get_ecef(rv, pointing, ExtraOutput(); kwargs...)[1]
end

# ╔═╡ bcf6ae44-fa8c-4d89-9fb9-01019098d981
#=╠═╡
let
	ref_uv = SA_F64[0.1,0.2]
	@benchmark get_ecef($sv, $ref_uv)
end
  ╠═╡ =#

# ╔═╡ 90d02a56-604d-4267-b827-7ca263892ce2
#=╠═╡
let
	ref_uv = SA_F64[0.1,0.2]
	@benchmark get_ecef($sv, $ref_uv, ExtraOutput())
end
  ╠═╡ =#

# ╔═╡ cc1c1137-a253-49de-8293-5819236a00cf
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Get LLA
"""
  ╠═╡ =#

# ╔═╡ 99fa5fab-70d0-4b85-9a0c-ecbac0441399
begin
"""
	get_lla(rv::ReferenceView, pointing::Point2D[, ::ExtraOutput]; pointing_type::Symbol=:uv, h = 0.0, face = rv.face, R = nothing)

Computes the LLA coordinates of the point that is seen by `rv` in the direction
specified by `pointing` and is located at a target altitude `h` [m] above the
earth's surface.

If a valid point can not be found because either earth is blocking the view or
no point at altitude `h` can be seen from the provided pointing direction in the
`rv` local CRS (also accounting for desired face), the function returns a
SVector{3, Float64} filled win NaNs.

`pointing_type` is used to select whether the output is returned in UV or
ThetaPhi [rad] coordinates. The following symbols are supported for this kwarg:
- `:thetaphi`, `:ThetaPhi` and `:θφ` can be used to represent pointing in ThetaPhi [rad]
- `:UV` and `:uv` can be used to represent pointing in UV

When called with an instance of `TelecomUtils.ExtraOutput` as last argument, the
function also returns the coordinated of the identified point in the local CRS
of `rv`.

For details on how to modify the reference pointing direction using the kwargs
`face` and `R` look at the documentation of [`get_range`](@ref) 

See also: [`ReferenceView`](@ref), [`get_range`](@ref), [`get_pointing`](@ref),
[`get_lla`](@ref), [`get_ecef`](@ref), [`get_era`](@ref),
[`get_distance_on_earth`](@ref), [`get_visibility`](@ref), [`get_mutual_visibility`](@ref).
"""
function get_lla(rv::ReferenceView,pointing::Point2D, eo::ExtraOutput; kwargs...)
	ecef, xyz = get_ecef(rv, pointing, eo; kwargs...)
	return LLAfromECEF(rv.ellipsoid)(ecef), xyz
end

# Single Output
get_lla(rv, pointing; kwargs...) = get_lla(rv, pointing, ExtraOutput(); kwargs...)[1]
end

# ╔═╡ 6ea89a6e-5995-4339-af32-cb4bf8604e25
#=╠═╡
let
	ref_uv = SA_F64[0.1,0.2]
	@benchmark get_lla($sv, $ref_uv)
end
  ╠═╡ =#

# ╔═╡ 3f7c8543-7126-4483-9d36-5fd602ae679a
#=╠═╡
let
	ref_uv = SA_F64[0.1,0.2]
	@benchmark get_lla($sv, $ref_uv, ExtraOutput())
end
  ╠═╡ =#

# ╔═╡ d30176d8-4faf-4ea4-bf11-326b42d15eac
#=╠═╡
let
	sat_lla = LLA(0°,0°,600km)
	sv = SatView(sat_lla, em)
	uv = (0,1)
	face = :NegativeY
	ref_uv = get_lla(sv, uv; face = face)
end
  ╠═╡ =#

# ╔═╡ 8ef02fae-4a93-45de-9040-e045b2e4f59f
#=╠═╡
let
	sat_lla = LLA(0°,0°,600km)
	sv = SatView(sat_lla, em)
	face = :NegativeY
	ref_uv = get_pointing(sv, LLA(0,0,0); face = face)
end
  ╠═╡ =#

# ╔═╡ 2f030209-7388-4d42-aaf4-5436b5ca5bfb
#=╠═╡
let
	sat_lla = LLA(0°,0°,600km)
	sv = SatView(sat_lla, em)
	uv = (0,0)
	face = :NegativeY
	ref_uv = get_lla(sv, uv; h = 700e3, face = face)
end
  ╠═╡ =#

# ╔═╡ 3f030ab2-d9cb-40de-8017-bfe734634551
#=╠═╡
let    
	sp_ell = SphericalEllipsoid()
    lla2ecef = ECEFfromLLA(sp_ell)

    em = EarthModel(sp_ell)
    sat_lla = LLA(0,0,700km)

    sv = SatView(sat_lla, em)
	lla_ref = LLA(1°, 1°, 1km)
	ecef_ref = lla2ecef(lla_ref)
	ref_uv = get_pointing(sv, lla_ref)
	@test get_lla(sv, ref_uv; h = lla_ref.alt) ≈ lla_ref
	@test get_ecef(sv, ref_uv; h = lla_ref.alt) ≈ ecef_ref
end
  ╠═╡ =#

# ╔═╡ 4fafa550-580b-4851-855d-0e46bf14a357
#=╠═╡
let    
	sp_ell = SphericalEllipsoid()
    lla2ecef = ECEFfromLLA(sp_ell)

    em = EarthModel(sp_ell)
    sat_lla = LLA(0,0,700km)

    rv = SatView(sat_lla, em)
	uv = (0.1,0)
	@benchmark get_lla($rv, $uv; h = 100e3)
end
  ╠═╡ =#

# ╔═╡ 8bc60d8d-7b54-4dce-a3e4-e336c0b16d4e
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Get ERA
"""
  ╠═╡ =#

# ╔═╡ f59d34bc-3b19-40ab-b0c6-986e5fa62304
begin
"""
	get_era(uv::UserView, target::Union{LLA, Point3D, ReferenceView}[, ::ExtraOutput]; face = rv.face, R = nothing)

Computes the ERA coordinates of the provided `target` as seen from the `UserView` `uv`.
`target` can be given either as `LLA`/`ECEF` coordinates or directly as another `ReferenceView`

For details on how to modify the reference pointing direction using the kwargs `face` and `R` look at the documentation of [`get_range`](@ref) 

See also: [`ReferenceView`](@ref), [`get_range`](@ref), [`get_pointing`](@ref), [`get_lla`](@ref), [`get_ecef`](@ref), [`get_distance_on_earth`](@ref).
"""
function get_era(uv::UserView, target::Union{LLA, Point3D, ReferenceView}; face = uv.face, R = nothing)	
	# Get ECEF coordinates of the target
	ecef = if target isa LLA
		ECEFfromLLA(uv.ellipsoid)(target)
	elseif target isa ReferenceView
		@assert target.earthmodel === uv.earthmodel "When computing visibility angles between ReferenceView objects, the `EarthModel` they use internally must be the same"
		target.ecef
	else
		target
	end
	_R = isnothing(R) ? inv(uv.R * face_rotation(face)) : R 
	ERAfromECEF(uv.ecef, _R, uv.ellipsoid)(ecef)
end
get_era(sv::SatView, args...; kwargs...) = error("The first argument to `get_era` must be of type `UserView`.")
end

# ╔═╡ ee657a11-c976-4128-8bb4-2336a5ecd319
#=╠═╡
# We test that a non-visible point is nan
@test get_era(UserView(LLA(40°,-39°,0),em),LLA(0,0,700km)) |> isnan
  ╠═╡ =#

# ╔═╡ 2ad13505-0c60-4ccb-b536-e865c24a0396
#=╠═╡
# We test that a visible point is not nan
@test get_era(UserView(LLA(0,0,500km),em),LLA(0,0,600km)) ≈ ERA(90°, 100km, 0°)
  ╠═╡ =#

# ╔═╡ 97c3ab73-5d2b-4871-aaa2-f8d7f1a7204d
#=╠═╡
let
	uv = UserView(LLA(0,0,0), em)
	@benchmark get_era($uv,$(sv.lla))
end
  ╠═╡ =#

# ╔═╡ d34f768d-3a7f-4fc8-a311-098c9aea789a
#=╠═╡
let
	uv = UserView(LLA(0,0,0), em)
	@benchmark get_era($uv,$(sv.ecef))
end
  ╠═╡ =#

# ╔═╡ 036defd5-7a49-44ef-a55f-8d9a12d43c0c
#=╠═╡
let
	uv = UserView(LLA(0,0,0), em)
	@benchmark get_era($uv,$sv)
end
  ╠═╡ =#

# ╔═╡ bbf6f990-40b3-471f-a46c-61f5fd6f5824
# ╠═╡ skip_as_script = true
#=╠═╡
# Visual Test ERA
begin
	function meshgrid(xin,yin)
		nx=length(xin)
		ny=length(yin)
		xout=zeros(ny,nx)
		yout=zeros(ny,nx)
		for jx=1:nx
		    for ix=1:ny
		        xout[ix,jx]=xin[jx]
		        yout[ix,jx]=yin[ix]
		    end
		end
		return xout,yout
	end;
	
	svTest = SatView(LLA(0,0,700km),em)
	gridRes = 0.5
	x_plot = -180:gridRes:180
	y_plot = -90:gridRes:90

	lonMesh,latMesh = meshgrid(x_plot, y_plot)
	el_plot = fill(NaN,size(lonMesh,1),size(lonMesh,2))
	az_plot = fill(NaN,size(lonMesh,1),size(lonMesh,2))
	r_plot = fill(NaN,size(lonMesh,1),size(lonMesh,2))
	for row = 1:size(latMesh,1)
		for col = 1:size(latMesh,2)
			era = get_era(UserView(LLA(y_plot[row]*°, x_plot[col] * °, 0), em),svTest)
			el_plot[row,col] = rad2deg(era.el)
			az_plot[row,col] = rad2deg(era.az)
			r_plot[row,col] = era.r
		end
	end
end
  ╠═╡ =#

# ╔═╡ 74b99a07-9e73-4532-a50d-10221c47f324
#=╠═╡
let
	el_plot = heatmap(
		x = x_plot,
		y = y_plot,
		z = el_plot
	)
	
	plot(el_plot, Layout(
					yaxis_title = "LAT",
					xaxis_title = "LON",
					title = "Elevation Test"))
end
  ╠═╡ =#

# ╔═╡ 1c0aa81c-efa2-4ba4-a3b2-70276d76c4f1
#=╠═╡
let
	az_plot = heatmap(
		x = x_plot,
		y = y_plot,
		z = az_plot
	)
	
	plot(az_plot, Layout(
					yaxis_title = "LAT",
					xaxis_title = "LON",
					title = "Azimuth Test"))
end
  ╠═╡ =#

# ╔═╡ 5023b71d-219e-4f2f-b319-e9899e9702ac
#=╠═╡
let
	r_plot = heatmap(
		x = x_plot,
		y = y_plot,
		z = r_plot
	)
	
	plot(r_plot, Layout(
					yaxis_title = "LAT",
					xaxis_title = "LON",
					title = "Range Test"))
end
  ╠═╡ =#

# ╔═╡ 64370881-a469-4748-97c5-ec27199d529b
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Get Distance on Earth
"""
  ╠═╡ =#

# ╔═╡ 2efb01b8-16b1-4186-94f4-cdfbca1310de
@addmethod geod_inverse(sv::ReferenceView, args...) = geod_inverse(sv.geod,args...)

# ╔═╡ 04ceecfa-2174-468c-933c-30b999e5b1be
function gesu(lla1::LLA, lla2::LLA, ::ExtraOutput, em::EarthModel = EarthModel())
	# @assert abs(lla1.alt) < 1 && abs(lla2.alt) < 1 "The altitude of the provided LLA inputs must be 0 (± 1m)"
	geod_inverse(em.geod, lla1, lla2)[1]
end

# ╔═╡ e0915eab-a53d-4fb2-9029-83793073ac3c
begin
"""
	get_distance_on_earth(sv::SatView, p1::Point2D, p2::Point2D[, ::ExtraOutput]; pointing_type::Symbol=:uv, face = sv.face, R)
	get_distance_on_earth(p1::LLA, p2::LLA[, ::ExtraOutput]; em = EarthModel())
	get_distance_on_earth(p1::UserView, p2::UserView)
Computes the distance [m] between the points on the earth surface (`lla1` and `lla2`) using the reference earth model used by `sv`.

If the points are not provided as LLA instances, but as angular directions (`p1` and `p2`), `lla1` and `lla2` as first computed from `p1` and `p2` using the SatView object `sv` as reference.

When called with angular directions, the optional argument `kind` is used to select whether the pointing is expressed in ThetaPhi (`kind ∈ (:ThetaPhi, :thetaphi, :θφ)`) [rad] or UV coordinates.

If an instance of `ExtraOutput` is provided as 4th argument, the function also returns the (forward) azimuth angle between `lla1` and `lla2` (2nd output, [deg]) and the azimuth angle between `lla2` and `lla1` (third output, [deg])

See also: [`SatView`](@ref), [`get_range`](@ref), [`get_pointing`](@ref), [`get_lla`](@ref), [`get_ecef`](@ref), [`geod_inverse`](@ref), [`get_nadir_beam_diameter`](@ref), [`ExtraOutput`](@ref).
"""
function get_distance_on_earth(lla1::LLA, lla2::LLA, ::ExtraOutput; em::EarthModel = EarthModel())
	@assert abs(lla1.alt) < 1 && abs(lla2.alt) < 1 "The altitude of the provided LLA inputs must be 0 (± 1m)"
	geod_inverse(em, lla1, lla2)
end

# Version with SatView as first argument	
function get_distance_on_earth(sv::SatView, p1::Point2D, p2::Point2D, eo::ExtraOutput; kwargs...)
	lla1 = get_lla(sv, p1; kwargs...)
	lla2 = get_lla(sv, p2; kwargs...)
	get_distance_on_earth(lla1, lla2, eo; em = sv.earthmodel)
end

# Versions with At least one of the two inputs as ReferenceView
get_distance_on_earth(rv::ReferenceView, p::LLA, eo::ExtraOutput) = get_distance_on_earth(rv.lla, p, eo; em = rv.earthmodel)
	
get_distance_on_earth(p::LLA, rv::ReferenceView, eo::ExtraOutput) = get_distance_on_earth(p, rv.lla, eo; em = rv.earthmodel)
	
function get_distance_on_earth(rv1::ReferenceView, rv2::ReferenceView, eo::ExtraOutput)
	@assert rv1.earthmodel === rv2.eartmodel "The two provided ReferenceView must refer to the same EarthModel"
	get_distance_on_earth(rv1.lla, rv2.lla, eo; em = rv1.earthmodel)
end

# Versions with single output
get_distance_on_earth(p1, p2; kwargs...) = get_distance_on_earth(p1, p2, ExtraOutput(); kwargs...)[1]

get_distance_on_earth(sv::SatView, p1::Point2D, p2::Point2D; kwargs...) = get_distance_on_earth(sv, p1, p2, ExtraOutput(); kwargs...)[1]
end

# ╔═╡ a4b39d91-7ac9-4253-b66c-bb541647f4e9
#=╠═╡
let
	lla1 = LLA(1°, 2°, 0km)
	lla2 = LLA(1°, 1°, 0km)
	@benchmark get_distance_on_earth($lla1, $lla2)
end
  ╠═╡ =#

# ╔═╡ 2612961b-0e1e-4615-8959-74ab3bc919f9
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## Get Nadir Beam Diameter
"""
  ╠═╡ =#

# ╔═╡ 30959832-9eb2-48c5-83d5-776d336c9aa7
"""
$SIGNATURES
Computes the diameter [m] on earth of a beam pointed at nadir from the satellite position identified by `sv` and assuming a 3db beamwidth identified by the scan angle `scan_3db` [deg] at which the beam pattern is 3dB below the peak.

The computation computes the diameter along the `U` direction and takes into account the reference ellipsoid of `sv`, so the resulting diameter is dependent on the satellite lat/lon position

See also: [`SatView`](@ref), [`get_range`](@ref), [`get_pointing`](@ref), [`get_lla`](@ref), [`get_ecef`](@ref), [`geod_inverse`](@ref), [`get_distance_on_earth`](@ref). 
"""
function get_nadir_beam_diameter(sv, scan_3db)
	uv_radius = sind(scan_3db)
	p1 = (uv_radius, 0)
	p2 = (-uv_radius, 0)
	get_distance_on_earth(sv, p1, p2)
end

# ╔═╡ b9dacaaf-b55c-46c8-8fd0-ad520505ecbb
begin
	export ReferenceView, SatView, UserView, change_position!, change_attitude!, change_reference_face!, get_range, get_era, get_mutual_pointing, get_pointing, get_lla, get_ecef, get_distance_on_earth, get_nadir_beam_diameter, crs_rotation
	export get_visibility, get_mutual_visibility
end

# ╔═╡ 4af7a092-8f42-4aef-9c09-feab8ebc1d87
# ╠═╡ skip_as_script = true
#=╠═╡
get_nadir_beam_diameter(SatView(LLA(50°,0°,735km), EarthModel()), 55)
  ╠═╡ =#

# ╔═╡ c02d0705-6647-4a44-8ae8-fc256f18c4ce
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
# Tests
"""
  ╠═╡ =#

# ╔═╡ d15726ab-5a28-4a24-b5ed-b3c8ecb6c581
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
## nadir beam diameter
"""
  ╠═╡ =#

# ╔═╡ ae686da9-45d5-4fc2-9cbd-2d828d792407
# ╠═╡ skip_as_script = true
#=╠═╡
@test get_nadir_beam_diameter(SatView(LLA(90°,0°,735km), EarthModel()), 55) ≈ get_nadir_beam_diameter(SatView(LLA(0°,0°,735km), EarthModel()), 55)
  ╠═╡ =#

# ╔═╡ 71d3f92e-d143-40dc-8701-37f9053766ef
# ╠═╡ skip_as_script = true
#=╠═╡
@test get_nadir_beam_diameter(SatView(LLA(90°,0°,735km), EarthModel(wgs84_ellipsoid)), 55) ≉ get_nadir_beam_diameter(SatView(LLA(0°,0°,735km), EarthModel(wgs84_ellipsoid)), 55)
  ╠═╡ =#

# ╔═╡ 1772993e-0a22-4e5d-ae31-bd1115a4313c
md"""
## Get Pointing
"""

# ╔═╡ 0af6931b-0738-4b36-af7b-aca124add95c
# ╠═╡ skip_as_script = true
#=╠═╡
let
	rv1 = SatView(LLA(0,0,600km), EarthModel())
	rv2 = SatView(LLA(0.1,0.1,600km), EarthModel())
	get_pointing(rv1, rv2)
end
  ╠═╡ =#

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
PlutoDevMacros = "a0499f29-c39b-4c5c-807c-88074221b949"
PlutoExtras = "ed5d0301-4775-4676-b788-cf71e66ff8ed"
PlutoPlotly = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
PlutoTest = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
BenchmarkTools = "~1.3.2"
PlutoDevMacros = "~0.5.1"
PlutoExtras = "~0.7.4"
PlutoPlotly = "~0.3.6"
PlutoTest = "~0.2.2"
PlutoUI = "~0.7.51"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.0"
manifest_format = "2.0"
project_hash = "beea89e16737ecde1ef0f3f17f87754b793e3f53"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "d9a9701b899b30332bbcb3e1679c41cce81fb0e8"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.3.2"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "be6ab11021cd29f0344d5c4357b163af05a48cba"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.21.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "600cc5508d66b78aae350f7accdb58763ac18589"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.10"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.2+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

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

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

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

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "0a1b7c2863e44523180fdb3146534e265a91870b"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.23"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

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

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

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
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "478ac6c952fddd4399e71d4779797c538d0ff2bf"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.8"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.0"

[[deps.PlotlyBase]]
deps = ["ColorSchemes", "Dates", "DelimitedFiles", "DocStringExtensions", "JSON", "LaTeXStrings", "Logging", "Parameters", "Pkg", "REPL", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "56baf69781fc5e61607c3e46227ab17f7040ffa2"
uuid = "a03496cd-edff-5a9b-9e67-9cda94a718b5"
version = "0.8.19"

[[deps.PlutoDevMacros]]
deps = ["HypertextLiteral", "InteractiveUtils", "MacroTools", "Markdown", "Pkg", "Random", "TOML"]
git-tree-sha1 = "2ec9ca2a56ab69334ab54c79c347a9d04afae9f5"
uuid = "a0499f29-c39b-4c5c-807c-88074221b949"
version = "0.5.3"

[[deps.PlutoExtras]]
deps = ["AbstractPlutoDingetjes", "HypertextLiteral", "InteractiveUtils", "Markdown", "OrderedCollections", "PlutoDevMacros", "PlutoUI", "REPL"]
git-tree-sha1 = "15e75e48e51416d33bab70943923a62a0b63f137"
uuid = "ed5d0301-4775-4676-b788-cf71e66ff8ed"
version = "0.7.4"

[[deps.PlutoPlotly]]
deps = ["AbstractPlutoDingetjes", "Colors", "Dates", "HypertextLiteral", "InteractiveUtils", "LaTeXStrings", "Markdown", "PlotlyBase", "PlutoUI", "Reexport"]
git-tree-sha1 = "dec81dcd52748ffc59ce3582e709414ff78d947f"
uuid = "8e989ff0-3d88-8e9f-f020-2b208a939ff0"
version = "0.3.6"

[[deps.PlutoTest]]
deps = ["HypertextLiteral", "InteractiveUtils", "Markdown", "Test"]
git-tree-sha1 = "17aa9b81106e661cffa1c4c36c17ee1c50a86eda"
uuid = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
version = "0.2.2"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "b478a748be27bd2f2c73a7690da219d0844db305"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.51"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "259e206946c293698122f63e2b513a7c99a244e8"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.1.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "7eb1686b4f04b82f96ed7a4ea5890a4f0c7a09f1"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.0"

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
deps = ["SHA", "Serialization"]
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

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

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

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "ef28127915f4229c971eb43f3fc075dd3fe91880"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.2.0"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

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

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.Tricks]]
git-tree-sha1 = "aadb748be58b492045b4f56166b5188aa63ce549"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.7"

[[deps.URIs]]
git-tree-sha1 = "074f993b0ca030848b897beff716d93aca60f06a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.2"

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
version = "1.2.13+0"

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
# ╠═3bda5426-c0de-493f-9514-30b6fe762463
# ╠═7d23f727-1907-4965-a940-cc873f6b2191
# ╠═73e00cef-9734-439a-b89b-7c1d99aab74e
# ╠═dcc81988-903b-4707-a70c-09c38682c80f
# ╟─3fd1046c-fabf-4264-9638-ba41301b1804
# ╠═7729ce27-df74-4393-ab70-c4e2864c85f5
# ╠═a57e3983-21de-4a2e-a227-8265fee6b56b
# ╠═b9dacaaf-b55c-46c8-8fd0-ad520505ecbb
# ╟─a5112e66-c2a2-4ed2-9951-5f97bc1745d5
# ╟─9e65ad17-95bd-46fa-bc18-9d5e8c501d9a
# ╠═f1d1295e-7fa1-44d0-bdac-8b7830da8a61
# ╠═9e8f0786-1daa-4b9e-9172-fc0767582c7e
# ╠═93222642-f2a6-4de7-8c92-21c96ef009a4
# ╠═cb0b070b-f70a-458e-bf72-0a4d2e93ec41
# ╠═28c868ca-9954-47ff-8949-a41fb7fc6d41
# ╟─f4c3876a-a81b-42f3-870a-43526e4c116e
# ╠═cc51ab70-0167-477f-a62d-88567e94fed9
# ╠═c54622f5-835d-4957-bdc6-88afed9c9d2a
# ╠═90811965-8bf0-4a4c-86f4-9bc5c86be71c
# ╠═278a536f-82dc-4f4e-8000-cc6a064ab2ee
# ╟─5dc71a85-a20d-4448-98b4-5065a249df1d
# ╠═c6ee08ba-3546-48ea-9801-edc00dfd25f0
# ╠═5372f3fe-699a-4f00-8e8e-36cbea224963
# ╠═a34ca715-e7f2-4fa7-ba94-8ad3f9b1f4cd
# ╠═227c0f3e-6e45-48ac-92d4-d1c1c730e4e0
# ╠═6eb9424e-3dd2-46d4-b4d2-81596bb81668
# ╠═5f1fd82d-f441-4a1b-9840-773a8635d3db
# ╟─fb418edb-4937-4668-acf3-e22e3a818a99
# ╠═12f31d2c-1eac-47bc-bda2-a32395ff8265
# ╟─091e4ec2-ea9e-411e-8f39-73aeb73c0214
# ╠═3c2781b8-e8df-446f-95ed-597112edabec
# ╠═b17a4660-3938-4816-8182-6feae471b50d
# ╠═b4777d24-5ea7-4524-8d2a-ac5d070d58d5
# ╠═cddcb68d-5400-4cb9-90ef-b3e149900f8a
# ╠═64919bbb-3b52-48e1-a7ea-c0916c47994d
# ╟─5c1d0beb-c00c-40a8-b255-c887be99e706
# ╠═f127481d-c60e-43d0-ab9f-4d4f35984015
# ╠═3a745709-8f5b-4f22-848a-2f9754ab27d8
# ╟─851f1994-22f4-4347-b106-2fa3ea75ebf2
# ╠═74d39d15-8d6c-4fb8-8d04-18f32c393ad4
# ╠═8b68c7c2-5cbd-4edd-9739-0d2e8c7f5449
# ╠═84769564-8ba8-46f5-b494-b0689d9abd65
# ╠═33e4c937-4443-4d97-a3ed-81479ede1e11
# ╠═0cde0a71-7f27-4290-88cd-2cccf627926b
# ╠═bd62bdd6-4de4-449c-b5f1-fb1b4f695cda
# ╠═21bc362a-960c-4c5f-9118-7e1451edb996
# ╠═fd3a1dcb-2e64-47f8-a932-452742090ac1
# ╠═3381783a-b6d7-43c7-947d-c1de4437e25b
# ╠═449b49de-2951-41fc-ba46-89eaa6c52e79
# ╠═f99b984a-b7cf-4318-8f74-aacb644bc146
# ╠═b389a1db-c444-481c-b168-3069c7367276
# ╟─39a1850b-f64a-4157-8f07-d7a78918fea1
# ╠═83634223-87d0-4c31-801a-af8a7f9f678a
# ╠═d4ff395e-55bd-4240-89e3-a4be2dcf9ebe
# ╠═cb4d8cd7-ff3e-43e6-ad38-8db343087214
# ╠═f3b8b9ba-c428-4868-96b9-1983eabb3774
# ╠═951fc2e2-ada9-4aad-876a-4cbd46200f6c
# ╠═a93e2354-34d1-4d1d-ac3c-f4c99099820a
# ╠═46eb3a4d-c80b-41a0-9333-aaf6411a010c
# ╠═abb3a7ee-ff60-421c-8a85-23c3ef1ce6af
# ╠═b75df549-ed46-4def-9fe0-e63ac496b3f4
# ╠═314c4c2a-6f3b-4a7d-aa83-9073eb67897b
# ╟─a61e66cf-b430-4d84-a753-faf42f4b6337
# ╠═98f3f83d-fcde-48ba-8855-c30643776d81
# ╠═cda9503a-5262-4318-9f5d-aea7910ab42e
# ╠═868bb972-aa0b-484d-b906-1628412b1646
# ╠═83098fe4-2089-44be-8db9-6e2db74d886d
# ╟─1de548c0-3830-4161-b80c-be72dd894e85
# ╠═664be7e6-28f5-4ab3-b8f0-b7618bcd5fe5
# ╠═387d2c76-1a08-441c-97fc-7b0a90a95c9a
# ╠═864cc72a-6252-45a6-8ed4-81d960c0b080
# ╠═ecd8d6ac-27f6-4a50-81ff-9c8f0b7fed47
# ╠═fe0680dd-95fb-40a0-8063-97c009f552dd
# ╟─f0565008-d355-4ba4-9269-38a8dc4d88bf
# ╠═abd76e81-2464-4119-a9bf-6046805f8e57
# ╠═de0372e1-da0c-4da9-8236-790ecd667d07
# ╠═a9051f97-2561-4317-b215-eac2b4f20aa4
# ╟─1f27b72f-9a3b-4732-a98e-d216af067072
# ╠═12330b6f-97b0-4efb-9885-49758bc2f127
# ╠═bcf6ae44-fa8c-4d89-9fb9-01019098d981
# ╠═90d02a56-604d-4267-b827-7ca263892ce2
# ╟─cc1c1137-a253-49de-8293-5819236a00cf
# ╠═99fa5fab-70d0-4b85-9a0c-ecbac0441399
# ╠═6ea89a6e-5995-4339-af32-cb4bf8604e25
# ╠═3f7c8543-7126-4483-9d36-5fd602ae679a
# ╠═d30176d8-4faf-4ea4-bf11-326b42d15eac
# ╠═8ef02fae-4a93-45de-9040-e045b2e4f59f
# ╠═2f030209-7388-4d42-aaf4-5436b5ca5bfb
# ╠═3f030ab2-d9cb-40de-8017-bfe734634551
# ╠═4fafa550-580b-4851-855d-0e46bf14a357
# ╟─8bc60d8d-7b54-4dce-a3e4-e336c0b16d4e
# ╠═f59d34bc-3b19-40ab-b0c6-986e5fa62304
# ╠═ee657a11-c976-4128-8bb4-2336a5ecd319
# ╠═2ad13505-0c60-4ccb-b536-e865c24a0396
# ╠═97c3ab73-5d2b-4871-aaa2-f8d7f1a7204d
# ╠═d34f768d-3a7f-4fc8-a311-098c9aea789a
# ╠═036defd5-7a49-44ef-a55f-8d9a12d43c0c
# ╠═bbf6f990-40b3-471f-a46c-61f5fd6f5824
# ╠═74b99a07-9e73-4532-a50d-10221c47f324
# ╠═1c0aa81c-efa2-4ba4-a3b2-70276d76c4f1
# ╠═5023b71d-219e-4f2f-b319-e9899e9702ac
# ╠═64370881-a469-4748-97c5-ec27199d529b
# ╠═2efb01b8-16b1-4186-94f4-cdfbca1310de
# ╠═04ceecfa-2174-468c-933c-30b999e5b1be
# ╠═e0915eab-a53d-4fb2-9029-83793073ac3c
# ╠═a4b39d91-7ac9-4253-b66c-bb541647f4e9
# ╠═2612961b-0e1e-4615-8959-74ab3bc919f9
# ╠═30959832-9eb2-48c5-83d5-776d336c9aa7
# ╠═4af7a092-8f42-4aef-9c09-feab8ebc1d87
# ╟─c02d0705-6647-4a44-8ae8-fc256f18c4ce
# ╠═d15726ab-5a28-4a24-b5ed-b3c8ecb6c581
# ╠═ae686da9-45d5-4fc2-9cbd-2d828d792407
# ╠═71d3f92e-d143-40dc-8701-37f9053766ef
# ╟─1772993e-0a22-4e5d-ae31-bd1115a4313c
# ╠═0af6931b-0738-4b36-af7b-aca124add95c
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
