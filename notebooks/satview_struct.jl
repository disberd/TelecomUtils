### A Pluto.jl notebook ###
# v0.17.2

# using Markdown
# using InteractiveUtils

# ╔═╡ 13646410-4e96-11ec-3e3d-99763ba1aeea
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

# ╔═╡ 73e00cef-9734-439a-b89b-7c1d99aab74e
#=╠═╡ notebook_exclusive
begin
	using BenchmarkTools
	using PlutoTest
	using PlutoUtils
	using PlutoDevMacros
end
  ╠═╡ notebook_exclusive =#

# ╔═╡ 3bda5426-c0de-493f-9514-30b6fe762463
#=╠═╡ notebook_exclusive
md"""
# Packages
"""
  ╠═╡ notebook_exclusive =#

# ╔═╡ dcc81988-903b-4707-a70c-09c38682c80f
#=╠═╡ notebook_exclusive
ToC()
  ╠═╡ notebook_exclusive =#

# ╔═╡ 3fd1046c-fabf-4264-9638-ba41301b1804
#=╠═╡ notebook_exclusive
md"""
## load other notebook
"""
  ╠═╡ notebook_exclusive =#

# ╔═╡ 7729ce27-df74-4393-ab70-c4e2864c85f5
@plutoinclude "satview_transformations.jl" "all"

# ╔═╡ de735c56-612c-4ffd-8335-95f20a129390
#=╠═╡ notebook_exclusive
@macroexpand @plutoinclude "satview_transformations.jl" "all"
  ╠═╡ notebook_exclusive =#

# ╔═╡ a57e3983-21de-4a2e-a227-8265fee6b56b
#=╠═╡ notebook_exclusive
md"""
# Exports
"""
  ╠═╡ notebook_exclusive =#

# ╔═╡ 030e15c5-83a8-4a24-836a-96b6f4f0bb04
#=╠═╡ notebook_exclusive
md"""
# SatView
"""
  ╠═╡ notebook_exclusive =#

# ╔═╡ b966fa0c-dc52-4821-bc32-e78dd3272ce1
#=╠═╡ notebook_exclusive
md"""
We want to define a SatView mutable struct that represents a satellite and can be used to compute view angles, ranges and ERA for specific points on ground
"""
  ╠═╡ notebook_exclusive =#

# ╔═╡ e469d967-79e4-4ef2-b635-51a183cb12e7
begin
	"""
	SatView(lla::LLA,earthmodel::EarthModel)
	SatView(ecef::StaticVector{3},earthmodel::EarthModel)
Object representing the instantaneous position of a satellite and used to compute various view angles to and from points on ground, as well as inverse geodesic computations.

# Fields
$TYPEDFIELDS

If multiple satellites have to be tracked, the `EarthModel` instance `earthmodel` should be generated once and then passed to all the SatView instances to ensure that all the satellites are referring to the same earth model.\\
Doing so will enforce the same Ellipsoid is shared between all satellites even when it is changed from one of the SatView instances.

See also: [`change_position!`](@ref), [`get_range`](@ref), [`get_pointing`](@ref), [`get_lla`](@ref), [`get_ecef`](@ref), [`geod_inverse`](@ref), [`get_distance_on_earth`](@ref), [`get_nadir_beam_diameter`](@ref), [`ExtraOutput`](@ref).
	"""
	mutable struct SatView
		"ECEF coordinates of the current satellite position"
		ecef::SVector{3,Float64}
		"LLA coordinates of the current satellite position"
		lla::LLA
		"Reference EarthModel used for the projections and for the geodesic computations"
		earthmodel::EarthModel
		"Rotation Matrix to go from the nadir-pointing CRS to the ECEF CRS"
		R::RotMatrix3{Float64}
	end

	# Custom constructor
	function SatView(lla::LLA,em::EarthModel)
		ecef = ECEFfromLLA(em.ellipsoid)(lla)
		R = _rotation_matrix(:ECEFfromUV, lla.lat, lla.lon)
		SatView(ecef,lla,em,R)
	end
	function SatView(ecef::StaticVector{3},em::EarthModel)
		lla = LLAfromECEF(em.ellipsoid)(ecef)
		R = _rotation_matrix(:ECEFfromUV, lla.lat, lla.lon)
		SatView(ecef,lla,em,R)
	end
end

# ╔═╡ 4794b9dc-4297-402b-94a2-ed686584bb09
function Base.getproperty(sv::SatView, name::Symbol)
	if name ∈ (:ellipsoid, :geod)
		return getfield(getfield(sv,:earthmodel),name)
	else
		return getfield(sv, name)
	end
end

# ╔═╡ dd8d119d-550c-4217-923d-43aaf2b8327b
function Base.setproperty!(sv::SatView, name::Symbol, x)
	if name ∈ (:ellipsoid, :geod)
		return setproperty!(getfield(sv,:earthmodel),name, x)
	else
		return setfield!(sv, name, x)
	end
end

# ╔═╡ 6eb9424e-3dd2-46d4-b4d2-81596bb81668
#=╠═╡ notebook_exclusive
em = EarthModel()
  ╠═╡ notebook_exclusive =#

# ╔═╡ 170f2914-fdf9-46c8-a8e0-9130b046bd60
#=╠═╡ notebook_exclusive
sv = SatView(LLA(0,0,100km), em)
  ╠═╡ notebook_exclusive =#

# ╔═╡ e28cbfc7-408e-49b5-9aeb-bd01c32fba46
#=╠═╡ notebook_exclusive
let
	sv = SatView(LLA(0,0,100km), EarthModel())
	sv.ellipsoid = wgs84_ellipsoid
	sv.geod
end
  ╠═╡ notebook_exclusive =#

# ╔═╡ 5f1fd82d-f441-4a1b-9840-773a8635d3db
#=╠═╡ notebook_exclusive
@benchmark getproperty($sv, :geod)
  ╠═╡ notebook_exclusive =#

# ╔═╡ 091e4ec2-ea9e-411e-8f39-73aeb73c0214
#=╠═╡ notebook_exclusive
md"""
## Change Position
"""
  ╠═╡ notebook_exclusive =#

# ╔═╡ 41370c82-d32a-41ea-a21a-614574292c21
begin
	"""
		change_position!(sv::SatView, ecef, lla::LLA, R)
		change_position!(sv::SatView, lla::LLA)
		change_position!(sv::SatView, ecef::StaticVector{3})
Change the position of a [`SatView`](@ref) object `sv`, also returning as output the modified
`sv`.  The function mutates the `ecef`, `lla` and `R` fields of the `sv`
object with the values provided as arguments (when using the 1st method above).\\
If only `ecef` or `lla` coordinates are provided (2nd and 3rd method above), the remaining
two arguments are computed automatically.

One would normally use either the 2nd or 3rd mehtod so that the two missing components are
correctly computed by the function.\\
The first method avoids computations but does not validate that `ecef`, `lla` and `R` are
correct and refer to the same position in space. For this reason the first method should
only be used if those values are correctly pre-computed elsewhere and one wants to avoid the
duplicate computations. 

See also: [`SatView`](@ref), [`get_range`](@ref), [`get_pointing`](@ref), [`get_lla`](@ref),
[`get_ecef`](@ref), [`get_ecef`](@ref).
	"""
	function change_position!(sv::SatView, ecef, lla::LLA, R)
		setfield!(sv,:ecef,ecef)
		setfield!(sv,:lla,lla)
		setfield!(sv,:R,R)
		return sv
	end
	function change_position!(sv::SatView, lla::LLA)
		ecef = ECEFfromLLA(sv.ellipsoid)(lla)
		R = _rotation_matrix(:ECEFfromUV, lla.lat, lla.lon)
		change_position!(sv, ecef, lla, R)
	end
	function change_position!(sv::SatView, ecef::StaticVector{3})
		lla = LLAfromECEF(sv.ellipsoid)(ecef)
		R = _rotation_matrix(:ECEFfromUV, lla.lat, lla.lon)
		change_position!(sv, ecef, lla, R)
	end
end

# ╔═╡ f77dcbcd-f042-4f7c-b97f-de63637229d0
#=╠═╡ notebook_exclusive
@benchmark change_position!($sv,$(SA_F64[1e7,0,0]))
  ╠═╡ notebook_exclusive =#

# ╔═╡ 78a8e7a4-333d-44ca-a438-fd85d7078300
#=╠═╡ notebook_exclusive
@code_warntype change_position!(sv,SA_F64[1e7,0,0])
  ╠═╡ notebook_exclusive =#

# ╔═╡ 709c44c8-c580-4cf7-8376-9c513eb3bd53
#=╠═╡ notebook_exclusive
@benchmark change_position!($sv,LLA(0,0,100km))
  ╠═╡ notebook_exclusive =#

# ╔═╡ 84769564-8ba8-46f5-b494-b0689d9abd65
#=╠═╡ notebook_exclusive
md"""
## Get Range
"""
  ╠═╡ notebook_exclusive =#

# ╔═╡ 642f2ede-b154-4260-a959-0a47ca4793b7
"""
	get_range(sv::SatView,uv::StaticVector{2,<:Number},h::Real = 0.0)
	get_range(sv::SatView,uv::Tuple{<:Number,<:Number},h::Real = 0.0)
Get the range [in m] between the satellite and a given point identified by the uv pointing and the reference altitude above the ellipsoid.
"""
function get_range(sv::SatView,uv::StaticVector{2,<:Number},h::Real = 0.0)
	# Find the ecef coordinate of the target point on earth
	ecef = ECEFfromUV(sv.ecef, sv.R, sv.ellipsoid)(uv,h)
	# Return the distance betwteen the satellite and the point
	return norm(ecef-sv.ecef)
end

# ╔═╡ 7a75d351-8583-455f-89c4-2d50cf79ea96
get_range(sv::SatView,uv::Tuple{<:Number,<:Number},args...) = get_range(sv,SVector(uv),args...)

# ╔═╡ 556934d8-d3ee-4a43-8f74-0939c5431c6f
"""
	get_range(sv::SatView,lla::LLA)
	get_range(sv::SatView,lla::LLA, ::ExtraOutput)
	get_range(sv::SatView,ecef::StaticVector{3})
	get_range(sv::SatView,ecef::StaticVector{3}, ::ExtraOutput)
Get the range [in m] between the satellite and a given point identified either by its LLA or ECEF coordinates. Returns NaN if the point is not visible either because it is covered by the earth surface or it is *above* the satellite.

If called with an instance of the `ExtraOutput` struct as last argument, it also provides the WND coordinates of the target point as seen from the satellite

See also: [`SatView`](@ref), [`change_position!`](@ref), [`get_pointing`](@ref), [`get_lla`](@ref), [`get_ecef`](@ref).
"""
function get_range(sv::SatView, ecef::StaticVector{3}, ::ExtraOutput)
	Δecef = ecef - sv.ecef
	dist = norm(Δecef)
	# Find if the target point is below the satellite, we do this by checking the last coordinate of the WND coordinates of the point
	wnd = sv.R'Δecef
	wnd[3] < 0 && return NaN, wnd
	# We have to check that the given lla is visible from the satellite, this happens if either there is no intersection with earth in the direction of pointing, or if the first intersection happens for a range greater than th computed one
	t₁, t₂ = _intersection_solutions(Δecef/dist, sv.ecef, sv.ellipsoid.a, sv.ellipsoid.b)
	r = (isnan(t₁) || t₁ > dist-1e-5) ? dist : NaN
	return r,wnd
end

# ╔═╡ 68417643-fa77-4780-9890-b0dac95bdb7f
get_range(sv::SatView, ecef::StaticVector{3}) = get_range(sv,ecef,ExtraOutput())[1]

# ╔═╡ da78f52b-30b6-4faf-bcea-b665c10ff4fe
get_range(sv::SatView, lla::LLA, args...) = get_range(sv,ECEFfromLLA(sv.ellipsoid)(lla),args...)

# ╔═╡ 449b49de-2951-41fc-ba46-89eaa6c52e79
#=╠═╡ notebook_exclusive
get_range(sv,LLA(0°,5°,10km),ExtraOutput())
  ╠═╡ notebook_exclusive =#

# ╔═╡ e7443f5b-a1a8-4866-9a64-ce7587465911
#=╠═╡ notebook_exclusive
get_range(SatView(LLA(0,0,600km),em),(0,0))
  ╠═╡ notebook_exclusive =#

# ╔═╡ 7c07a3c1-c1ec-4b83-b7c6-251edf91273c
#=╠═╡ notebook_exclusive
@benchmark $get_range($sv,LLA(0,0,10km))
  ╠═╡ notebook_exclusive =#

# ╔═╡ 39a1850b-f64a-4157-8f07-d7a78918fea1
#=╠═╡ notebook_exclusive
md"""
## Get Pointing
"""
  ╠═╡ notebook_exclusive =#

# ╔═╡ 51987c04-18f5-46bb-a3ba-5f94907a7960
"""
	get_pointing(sv::SatView,lla::LLA,kind::Symbol=:uv)
	get_pointing(sv::SatView,ecef::StaticVector{3},kind::Symbol=:uv)
Provide the 2-D angular pointing at which the target point (specified as LLA or ECEF) is seen from the satellite.

`kind` is used to select whether the output should be given in UV or ThetaPhi coordinates. The result is provided as ThetaPhi [in rad] if `kind ∈ (:ThetaPhi, :thetaphi, :θφ)`

See also: [`SatView`](@ref), [`get_range`](@ref), [`get_pointing`](@ref), [`get_lla`](@ref), [`get_ecef`](@ref), [`get_distance_on_earth`](@ref).
"""
function get_pointing(sv::SatView,ecef::StaticVector{3},kind::Symbol=:uv)
	uv = UVfromECEF(sv.ecef,sv.R',sv.ellipsoid)(ecef)
	if kind ∈ (:ThetaPhi, :thetaphi, :θφ)
		return ThetaPhifromUV()(uv)
	else
		return uv
	end
end

# ╔═╡ a6db34bc-b846-49aa-8d57-fb32cdce1684
get_pointing(sv::SatView,lla::LLA,args...) = get_pointing(sv,ECEFfromLLA(sv.ellipsoid)(lla),args...)

# ╔═╡ 1758748c-fa4b-4414-a05d-a32970c7a94b
#=╠═╡ notebook_exclusive
get_pointing(SatView(LLA(0,0,600km),em),LLA(1°,1°,0),:thetaphi)
  ╠═╡ notebook_exclusive =#

# ╔═╡ cc1c1137-a253-49de-8293-5819236a00cf
#=╠═╡ notebook_exclusive
md"""
## Get LLA
"""
  ╠═╡ notebook_exclusive =#

# ╔═╡ 1f7bf45c-b33b-4bfe-b82d-05b908ce375e
"""
	get_lla(sv::SatView,pointing::Point2D,kind::Symbol=:uv; h = 0.0)
Computes the LLA coordinates of a point on earth seen from the satellite identified by `sv`, the angular pointing identified by `pointing` and located at an altitude `h` [m] above the reference Earth ellipsoid of `sv`.

The optional argument `kind` is used to select whether the pointing is expressed in ThetaPhi (`kind ∈ (:ThetaPhi, :thetaphi, :θφ)`) [rad] or UV coordinates.

See also: [`SatView`](@ref), [`get_range`](@ref), [`get_pointing`](@ref), [`get_lla`](@ref), [`get_ecef`](@ref), [`get_distance_on_earth`](@ref).
"""
function get_lla(sv::SatView,pointing::Point2D,kind::Symbol=:uv; h = 0.0)
	uv = if kind ∈ (:ThetaPhi, :thetaphi, :θφ)
		UVfromThetaPhi()(pointing)
	else
		pointing
	end
	lla = LLAfromUV(sv.ecef,sv.R,sv.ellipsoid)(uv)
end

# ╔═╡ 1f27b72f-9a3b-4732-a98e-d216af067072
#=╠═╡ notebook_exclusive
md"""
## Get ECEF
"""
  ╠═╡ notebook_exclusive =#

# ╔═╡ 948cc7a1-d85e-4cfe-b2e4-e047bcbac305
"""
	get_ecef(sv::SatView,pointing::Point2D,kind::Symbol=:uv; h = 0.0)
Computes the ECEF coordinates of a point on earth seen from the satellite identified by `sv`, the angular pointing identified by `pointing` and located at an altitude `h` [m] above the reference Earth ellipsoid of `sv`.

The optional argument `kind` is used to select whether the pointing is expressed in ThetaPhi (`kind ∈ (:ThetaPhi, :thetaphi, :θφ)`) [rad] or UV coordinates.

See also: [`SatView`](@ref), [`get_range`](@ref), [`get_pointing`](@ref), [`get_lla`](@ref), [`change_position!`](@ref), [`get_distance_on_earth`](@ref).
"""
function get_ecef(sv::SatView,pointing::Point2D,kind::Symbol=:uv; h = 0.0)
	uv = if kind ∈ (:ThetaPhi, :thetaphi, :θφ)
		UVfromThetaPhi()(pointing)
	else
		pointing
	end
	lla = ECEFfromUV(sv.ecef,sv.R,sv.ellipsoid)(uv)
end

# ╔═╡ 8bc60d8d-7b54-4dce-a3e4-e336c0b16d4e
#=╠═╡ notebook_exclusive
md"""
## Get ERA
"""
  ╠═╡ notebook_exclusive =#

# ╔═╡ b2cb0afd-1220-40bd-8e1b-6df35e3db2f1
function get_era(sv::SatView, lla::LLA)
	ERAfromECEF(lla;ellipsoid = sv.ellipsoid)(sv.ecef)
end

# ╔═╡ 8fccb117-2048-4607-8db1-f8df7f5ef156
function get_era(sv::SatView, ecef::StaticVector{3})
	ERAfromECEF(ecef;ellipsoid = sv.ellipsoid)(sv.ecef)
end

# ╔═╡ ee657a11-c976-4128-8bb4-2336a5ecd319
#=╠═╡ notebook_exclusive
# We test that a non-visible point is nan
@test get_era(SatView(LLA(0,0,600km),em),LLA(1°,105°,0)) |> isnan
  ╠═╡ notebook_exclusive =#

# ╔═╡ 2ad13505-0c60-4ccb-b536-e865c24a0396
#=╠═╡ notebook_exclusive
# We test that a non-visible point is nan
@test get_era(SatView(LLA(0,0,600km),em),LLA(0,0,500km)) ≈ ERA(90°, 100km, 0°)
  ╠═╡ notebook_exclusive =#

# ╔═╡ 97c3ab73-5d2b-4871-aaa2-f8d7f1a7204d
#=╠═╡ notebook_exclusive
@benchmark get_era($sv,LLA(0,0,0))
  ╠═╡ notebook_exclusive =#

# ╔═╡ 64370881-a469-4748-97c5-ec27199d529b
#=╠═╡ notebook_exclusive
md"""
## Get Distance on Earth
"""
  ╠═╡ notebook_exclusive =#

# ╔═╡ 2efb01b8-16b1-4186-94f4-cdfbca1310de
geod_inverse(sv::SatView, args...) = geod_inverse(sv.geod,args...)

# ╔═╡ e0915eab-a53d-4fb2-9029-83793073ac3c
"""
	get_distance_on_earth(sv::SatView, p1::Point2D, p2::Point2D, kind::Symbol=:uv)
	get_distance_on_earth(sv::SatView, p1::Point2D, p2::Point2D, ::ExtraOutput, kind::Symbol=:uv)
	get_distance_on_earth(sv::SatView, lla1::LLA, lla2::LLA)
	get_distance_on_earth(sv::SatView, lla1::LLA, lla2::LLA, ::ExtraOutput)
Computes the distance [m] between the points on the earth surface (`lla1` and `lla2`) using the reference earth model used by `sv`.

If the points are not provided as LLA instances, but as angular directions (`p1` and `p2`), `lla1` and `lla2` as first computed from `p1` and `p2` using the SatView object `sv` as reference.

When called with angular directions, the optional argument `kind` is used to select whether the pointing is expressed in ThetaPhi (`kind ∈ (:ThetaPhi, :thetaphi, :θφ)`) [rad] or UV coordinates.

If an instance of `ExtraOutput` is provided as 4th argument, the function also returns the (forward) azimuth angle between `lla1` and `lla2` (2nd output, [deg]) and the azimuth angle between `lla2` and `lla1` (third output, [deg])

See also: [`SatView`](@ref), [`get_range`](@ref), [`get_pointing`](@ref), [`get_lla`](@ref), [`get_ecef`](@ref), [`geod_inverse`](@ref), [`get_nadir_beam_diameter`](@ref), [`ExtraOutput`](@ref).
"""
function get_distance_on_earth(sv::SatView, p1::Point2D, p2::Point2D, eo::ExtraOutput, kind::Symbol=:uv)
	lla1 = get_lla(sv, p1, kind)
	lla2 = get_lla(sv, p2, kind)
	get_distance_on_earth(sv, lla1, lla2, eo)
end

# ╔═╡ 30d32b7a-c95c-4d80-a60a-a87b27b3bf3c
get_distance_on_earth(sv::SatView, p1::Point2D, p2::Point2D, kind::Symbol=:uv) = get_distance_on_earth(sv, p1, p2, ExtraOutput(), kind)[1]

# ╔═╡ 407101b2-c794-49b0-9f8b-07fb45b80ca9
get_distance_on_earth(sv::SatView, lla1::LLA, lla2::LLA, ::ExtraOutput) = geod_inverse(sv.geod, lla1, lla2)

# ╔═╡ af71267d-b5ee-46b7-bf8d-d740033d35e0
get_distance_on_earth(sv::SatView, lla1::LLA, lla2::LLA) = geod_inverse(sv.geod, lla1, lla2, ExtraOutput())[1]

# ╔═╡ 2612961b-0e1e-4615-8959-74ab3bc919f9
#=╠═╡ notebook_exclusive
md"""
## Get Nadir Beam Diameter
"""
  ╠═╡ notebook_exclusive =#

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
export SatView, change_position!, get_range, get_era, get_pointing, get_lla, get_ecef, get_distance_on_earth, get_nadir_beam_diameter

# ╔═╡ 4af7a092-8f42-4aef-9c09-feab8ebc1d87
#=╠═╡ notebook_exclusive
get_nadir_beam_diameter(SatView(LLA(50°,0°,735km), EarthModel()), 55)
  ╠═╡ notebook_exclusive =#

# ╔═╡ c02d0705-6647-4a44-8ae8-fc256f18c4ce
#=╠═╡ notebook_exclusive
md"""
# Tests
"""
  ╠═╡ notebook_exclusive =#

# ╔═╡ d15726ab-5a28-4a24-b5ed-b3c8ecb6c581
#=╠═╡ notebook_exclusive
md"""
## nadir beam diameter
"""
  ╠═╡ notebook_exclusive =#

# ╔═╡ ae686da9-45d5-4fc2-9cbd-2d828d792407
#=╠═╡ notebook_exclusive
@test get_nadir_beam_diameter(SatView(LLA(90°,0°,735km), EarthModel()), 55) ≈ get_nadir_beam_diameter(SatView(LLA(0°,0°,735km), EarthModel()), 55)
  ╠═╡ notebook_exclusive =#

# ╔═╡ 71d3f92e-d143-40dc-8701-37f9053766ef
#=╠═╡ notebook_exclusive
@test get_nadir_beam_diameter(SatView(LLA(90°,0°,735km), EarthModel(wgs84_ellipsoid)), 55) ≉ get_nadir_beam_diameter(SatView(LLA(0°,0°,735km), EarthModel(wgs84_ellipsoid)), 55)
  ╠═╡ notebook_exclusive =#

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
PlutoDevMacros = "~0.3.10"
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
git-tree-sha1 = "7bb2558b40c6176ac5094542f7e01407fc3b38c1"
uuid = "a0499f29-c39b-4c5c-807c-88074221b949"
version = "0.3.10"

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
git-tree-sha1 = "1e0cb51e0ccef0afc01aab41dc51a3e7f781e8cb"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.20"

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
# ╟─3bda5426-c0de-493f-9514-30b6fe762463
# ╠═13646410-4e96-11ec-3e3d-99763ba1aeea
# ╠═73e00cef-9734-439a-b89b-7c1d99aab74e
# ╠═dcc81988-903b-4707-a70c-09c38682c80f
# ╟─3fd1046c-fabf-4264-9638-ba41301b1804
# ╠═7729ce27-df74-4393-ab70-c4e2864c85f5
# ╠═de735c56-612c-4ffd-8335-95f20a129390
# ╟─a57e3983-21de-4a2e-a227-8265fee6b56b
# ╠═b9dacaaf-b55c-46c8-8fd0-ad520505ecbb
# ╟─030e15c5-83a8-4a24-836a-96b6f4f0bb04
# ╟─b966fa0c-dc52-4821-bc32-e78dd3272ce1
# ╠═e469d967-79e4-4ef2-b635-51a183cb12e7
# ╠═4794b9dc-4297-402b-94a2-ed686584bb09
# ╠═dd8d119d-550c-4217-923d-43aaf2b8327b
# ╠═6eb9424e-3dd2-46d4-b4d2-81596bb81668
# ╠═170f2914-fdf9-46c8-a8e0-9130b046bd60
# ╠═e28cbfc7-408e-49b5-9aeb-bd01c32fba46
# ╠═5f1fd82d-f441-4a1b-9840-773a8635d3db
# ╟─091e4ec2-ea9e-411e-8f39-73aeb73c0214
# ╠═41370c82-d32a-41ea-a21a-614574292c21
# ╠═f77dcbcd-f042-4f7c-b97f-de63637229d0
# ╠═78a8e7a4-333d-44ca-a438-fd85d7078300
# ╠═709c44c8-c580-4cf7-8376-9c513eb3bd53
# ╟─84769564-8ba8-46f5-b494-b0689d9abd65
# ╠═642f2ede-b154-4260-a959-0a47ca4793b7
# ╠═7a75d351-8583-455f-89c4-2d50cf79ea96
# ╠═556934d8-d3ee-4a43-8f74-0939c5431c6f
# ╠═68417643-fa77-4780-9890-b0dac95bdb7f
# ╠═da78f52b-30b6-4faf-bcea-b665c10ff4fe
# ╠═449b49de-2951-41fc-ba46-89eaa6c52e79
# ╠═e7443f5b-a1a8-4866-9a64-ce7587465911
# ╠═7c07a3c1-c1ec-4b83-b7c6-251edf91273c
# ╟─39a1850b-f64a-4157-8f07-d7a78918fea1
# ╠═51987c04-18f5-46bb-a3ba-5f94907a7960
# ╠═a6db34bc-b846-49aa-8d57-fb32cdce1684
# ╠═1758748c-fa4b-4414-a05d-a32970c7a94b
# ╟─cc1c1137-a253-49de-8293-5819236a00cf
# ╠═1f7bf45c-b33b-4bfe-b82d-05b908ce375e
# ╟─1f27b72f-9a3b-4732-a98e-d216af067072
# ╠═948cc7a1-d85e-4cfe-b2e4-e047bcbac305
# ╟─8bc60d8d-7b54-4dce-a3e4-e336c0b16d4e
# ╠═b2cb0afd-1220-40bd-8e1b-6df35e3db2f1
# ╠═8fccb117-2048-4607-8db1-f8df7f5ef156
# ╠═ee657a11-c976-4128-8bb4-2336a5ecd319
# ╠═2ad13505-0c60-4ccb-b536-e865c24a0396
# ╠═97c3ab73-5d2b-4871-aaa2-f8d7f1a7204d
# ╟─64370881-a469-4748-97c5-ec27199d529b
# ╠═2efb01b8-16b1-4186-94f4-cdfbca1310de
# ╠═e0915eab-a53d-4fb2-9029-83793073ac3c
# ╠═30d32b7a-c95c-4d80-a60a-a87b27b3bf3c
# ╠═407101b2-c794-49b0-9f8b-07fb45b80ca9
# ╠═af71267d-b5ee-46b7-bf8d-d740033d35e0
# ╟─2612961b-0e1e-4615-8959-74ab3bc919f9
# ╠═30959832-9eb2-48c5-83d5-776d336c9aa7
# ╠═4af7a092-8f42-4aef-9c09-feab8ebc1d87
# ╟─c02d0705-6647-4a44-8ae8-fc256f18c4ce
# ╟─d15726ab-5a28-4a24-b5ed-b3c8ecb6c581
# ╠═ae686da9-45d5-4fc2-9cbd-2d828d792407
# ╠═71d3f92e-d143-40dc-8701-37f9053766ef
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
