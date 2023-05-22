### A Pluto.jl notebook ###
# v0.19.25

#> custom_attrs = ["hide-enabled"]

using Markdown
using InteractiveUtils

# ╔═╡ 0d343b68-f73e-11ed-17ef-7d0ceae106f8
begin
	using PlutoDevMacros
	using Distances
end

# ╔═╡ bca459c7-a29e-4ab7-ab81-ce3005f47a65
# ╠═╡ skip_as_script = true
#=╠═╡
begin
	using PlutoExtras
	using PlutoUI
	using BenchmarkTools
	using PlutoTest
end
  ╠═╡ =#

# ╔═╡ e38b8433-ed1b-4295-86cc-59ddb588867b
@fromparent begin
	import ^: UVfromThetaPhi, PointingType, ThetaPhi, UV, °, XYZfromUV, LLA, EarthModel, geod_inverse, SphericalEllipsoid, ThetaPhifromUV, to_radians, ValidAngle
	using >.StaticArrays
	using >.Proj4
	using >.Rotations
	using >.LinearAlgebra
end

# ╔═╡ 592722de-5f30-416a-9de5-24c8d417c532
md"""
# Packages
"""

# ╔═╡ 2d5099cb-a208-433d-8bf4-66f66a8baa28
md"""
## Load from Package
"""

# ╔═╡ 4edc2a99-3238-4fad-851e-30ca8cba94d1
#=╠═╡
ExtendedTableOfContents()
  ╠═╡ =#

# ╔═╡ 50729705-7e9f-43db-879f-c2f375c21163
md"""
# New Implementation
"""

# ╔═╡ acca8d87-6f17-4b3c-9efb-5bd84dcd7432
md"""
The implementations here are mostly taken from [this webpage](https://www.movable-type.co.uk/scripts/latlong.html), where formulas for various geodesic computations are presented.

The computation are for lat/lon and assume spherical earth. They can easily be adapted to compute angular distance between pointing considering that θ,φ and lat/lon are almost equivalent:
- θ (according to the ISO/physics convention) is equivalent to the co-latitude (π/2 - latitude)
- φ, when considering the ECI reference frame, is equivalent to the longitude.
"""

# ╔═╡ 9a03ffa9-af46-4e9b-b80d-69b4712490a1
# ╠═╡ disabled = true
#=╠═╡
let 
	a = 3
	@benchmark sincos($a)
end
  ╠═╡ =#

# ╔═╡ b4f84a0b-bcd0-4ccb-b664-f257e60b4c63
function haversine(ll1, ll2)
	φ1, λ1 = ll1
	φ2, λ2 = ll2
	Δφ = (φ2 - φ1)
	Δλ = (λ2 - λ1)
	a = (1 - cos(Δφ))/2 + cos(φ1) * cos(φ2) * (1 - cos(Δλ))/2
	c = 2 * asin(min(sqrt(a),one(a)))
end

# ╔═╡ f230e0c9-5430-459d-8158-748f55cefb2e
# ╠═╡ disabled = true
#=╠═╡
let
	ll1 = deg2rad.((30, 40))
	ll2 = deg2rad.((-30, 120))
	a = @benchmark haversine($ll1, $ll2)
	b = @benchmark SphericalAngle()(reverse($ll1), reverse($ll2))
	a,b
end
  ╠═╡ =#

# ╔═╡ 193516b5-6e61-41c9-9a58-a62d7f0791ae
AngleAxis(π/2, 1,0,0)

# ╔═╡ cd27bd1c-101f-49cc-aad4-2b568fbf1cbc
# ╠═╡ disabled = true
#=╠═╡
let
	θ = rand(100)
	y = @SVector rand(3)
	a = @benchmark map(x -> AngleAxis(x,$y..., false), $θ)
	b = @benchmark map(x -> RotXYZ(x, x, x), $θ)
	a,b
end
  ╠═╡ =#

# ╔═╡ 314d95df-8e5a-4533-9cb1-578cac4e8e4e
let
	x2u = XYZfromUV()
	u2x = inv(x2u)
	t2u = UVfromThetaPhi()
	u2t = inv(t2u)
	tp2ll(tp) = (π/2 - tp[1], tp[2])
	tp = (deg2rad(10), deg2rad(10))
	ll = rad2deg.(tp2ll(tp))
	uv = t2u(tp)
	xyz = x2u(uv,1)
	R = RotZ(deg2rad(5.7251))
	uv2 = u2x(R*xyz)
	tp2 = u2t(uv2)
	ll2 = rad2deg.(tp2ll(tp2))
	ll, ll2
end

# ╔═╡ 81981b91-2c82-478e-bdc9-109c9bbcfb6c
md"""
## dist and bearing
"""

# ╔═╡ 57c80a92-72b2-45ab-bd8d-e481cd0c4d6c
begin
function dist_and_bearing(ll1, ll2)
	φ1, λ1 = ll1
	φ2, λ2 = ll2
	Δλ = (λ2 - λ1)
	sφ₁, cφ₁ = sincos(φ1)
	sφ₂, cφ₂ = sincos(φ2)
	sΔλ, cΔλ = sincos(Δλ)

	dist = acos(sφ₁ * sφ₂ + cφ₁ * cφ₂ * cΔλ)
	bearing = atan(sΔλ * cφ₂, cφ₁ * sφ₂ - sφ₁ * cφ₂* cΔλ)

	(dist, bearing)
end
dist_and_bearing(ll1::LLA, ll2::LLA) = dist_and_bearing((ll1.lat, ll1.lon), (ll2.lat, ll2.lon))
end

# ╔═╡ 46f68c32-f646-4eca-ba9a-d21faa1e797c
function dist_and_bearing2(ll1, ll2)
	φ1, λ1 = ll1
	φ2, λ2 = ll2
	Δλ = (λ2 - λ1)
	Δφ = (φ2 - φ1)
	sφ₁, cφ₁ = sincos(φ1)
	sφ₂, cφ₂ = sincos(φ2)
	sΔλ, cΔλ = sincos(Δλ)

	
	a = (1 - cos(Δφ))/2 + cφ₁ * cφ₂ * (1 - cΔλ)/2
	dist = 2 * asin(min(sqrt(a),one(a)))
	bearing = atan(sΔλ * cφ₂, cφ₁ * sφ₂ - sφ₁ * cφ₂* cΔλ)

	(dist, bearing)
end

# ╔═╡ 61c00eee-7d16-464e-932a-b58b7c72e4c9
em = EarthModel(SphericalEllipsoid(1))

# ╔═╡ 724afbbc-90ab-40aa-a3b8-93ec6acac37d
let
	lla1 = LLA(80°, 10°, 0)
	lla2 = LLA(78.49161°, 40°, 0)
	a = geod_inverse(em, lla1, lla2)
	b = dist_and_bearing(lla1, lla2)
	a,rad2deg.(b)
end

# ╔═╡ 28124248-18e2-4542-a893-389dbf416129
md"""
## Add Offset
"""

# ╔═╡ 97c5cf9f-be20-49bf-a33e-0a0a2b03ba22
begin
function add_offset(latlon, dist, bearing)
	φ₁, λ₁ = latlon
	sφ, cφ = sincos(φ₁)
	sθ, cθ = sincos(bearing)
	sδ, cδ = sincos(dist)
	φ₂ = asin(sφ * cδ + cφ * sδ * cθ)
	λ₂ = λ₁ + atan(sθ * sδ * cφ, cδ - sφ * sin(φ₂))
	ll2 = φ₂, rem2pi(λ₂, RoundNearest)
end
add_offset(lla::LLA, dist, bearing) = LLA(add_offset((lla.lat, lla.lon), dist, bearing)..., lla.alt)
end

# ╔═╡ 9722d87d-e617-47c9-9fb6-4a7a5170a50f
let
	ll1 = deg2rad.((30, 40))
	ll2 = deg2rad.((30, 140))
	offset = dist_and_bearing(ll1, ll2)
	t = add_offset(ll1, offset...)
	rad2deg.(t), rad2deg.(offset)
end

# ╔═╡ d5743ccc-ca6b-4d7c-a6c8-7da0ff52131e
let
	lla1 = LLA(80°, 10°, 0)
	lla2 = LLA(78.49161°, 40°, 0)
	dist, azi1, azi2 = geod_inverse(em, lla1, lla2)
	add_offset(lla1, dist, deg2rad(azi1))	
end

# ╔═╡ e1cf709a-e62a-446d-a661-94d4f808220e
# ╠═╡ disabled = true
#=╠═╡
let
	lla1 = LLA(80°, 10°, 0)
	lla2 = LLA(78.49161°, 40°, 0)
	dist, azi1, azi2 = geod_inverse(em, lla1, lla2)
	@benchmark lla2tp(add_offset($lla1, $dist, deg2rad($azi1)))
end
  ╠═╡ =#

# ╔═╡ 3b895f54-5a15-4265-95b1-153af12b91d2
# ╠═╡ disabled = true
#=╠═╡
let
	ll1 = deg2rad.((30, 40))
	ll2 = deg2rad.((-30, 120))
	d = haversine(ll1, ll2)
	b = bearing(ll1, ll2)
	db = dist_and_bearing(ll1, ll2)
	db2 = dist_and_bearing2(ll1, ll2)
	(
		@benchmark((haversine($ll1, $ll2), bearing($ll1, $ll2))),
		@benchmark(dist_and_bearing($ll1, $ll2)),
		@benchmark(dist_and_bearing2($ll1, $ll2)),
		db, db2
	)
end
  ╠═╡ =#

# ╔═╡ 52cddfb6-754d-4109-95f7-52795a94297e


# ╔═╡ ba74055a-b73e-4c66-8c3f-3bb3f8e75a38
begin
function bearing(ll1, ll2)
	φ1, λ1 = ll1
	φ2, λ2 = ll2
	Δλ = (λ2 - λ1)
	sφ₁, cφ₁ = sincos(φ1)
	sφ₂, cφ₂ = sincos(φ2)
	sΔλ, cΔλ = sincos(Δλ)
	y = sΔλ * cφ₂
	x = cφ₁ * sφ₂ - sφ₁ * cφ₂* cΔλ
	b = atan(y,x)
end
bearing(ll1::LLA, ll2::LLA) = bearing((ll1.lat, ll1.lon), (ll2.lat, ll2.lon))
end

# ╔═╡ a1f2dcce-3180-41b6-b5b5-329e7f4fa1d5
#=╠═╡
let
	ll1 = deg2rad.((30, 40))
	ll2 = deg2rad.((30.001, 40))
	d = haversine(ll1, ll2)
	b = bearing(ll1, ll2)
	db = dist_and_bearing(ll1, ll2)
	(d,b), db
	T = SVector{2, Float64}
	@test T(d,b) ≈ T(db)
end
  ╠═╡ =#

# ╔═╡ a9a81d0c-e5fa-47e1-95ec-2050f7506df9
# ╠═╡ disabled = true
#=╠═╡
let
	ll1 = deg2rad.((30, 40))
	ll2 = deg2rad.((-30, 120))
	@benchmark bearing($ll1, $ll2)
	# a,b
end
  ╠═╡ =#

# ╔═╡ f4cd828b-92c0-4675-ba44-73a34a674052
# ╠═╡ disabled = true
#=╠═╡
let
	y = 1.2
	x = 0.4
	@benchmark atan($y,$x)
end
  ╠═╡ =#

# ╔═╡ 04d347ab-a061-414f-9fac-7371d4589971
let
	todeg(x) = @. rad2deg(x)*°
	ll1 = LLA(10°, 0°)
	ll2 = LLA(10°, 140°)
	dist_and_bearing(ll1, ll2) |> todeg
end

# ╔═╡ 57e6eb7f-adce-4e37-b7ef-040254783fea
# ╠═╡ disabled = true
#=╠═╡
let
	a = SA_F64[.2,.4]
	@benchmark XYZfromUV()($a, 1)
end
  ╠═╡ =#

# ╔═╡ 931261d6-7fd6-4639-a0bc-c931542d7c11
em

# ╔═╡ 761027d3-a447-48e9-bb24-85b8d57d6ce8
begin
function geod_directs(geod::Proj4.geod_geodesic, lonlat::AbstractVector{Cdouble}, azimuth::Cdouble, distance::Cdouble)
	lonlat = copy(lonlat)
	p = pointer(lonlat)
	azi = Ref{Cdouble}()
    ccall((:geod_direct, Proj4.libproj),Cvoid,(Ptr{Cvoid},Cdouble,Cdouble,Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble},
          Ptr{Cdouble}), pointer_from_objref(geod), lonlat[2], lonlat[1], azimuth, distance, p+sizeof(Cdouble), p, azi)
	lonlat, azi[]
end

# Version with LLA inputs
function geod_directs(geod::Proj4.geod_geodesic, lla1::LLA, azi, dist)
	lonlat1 = rad2deg.([lla1.lon,lla1.lat])
	geod_directs(geod,lonlat1,azi, dist)
end

# Version with EarthModel as first input
geod_directs(em::EarthModel, lla1::LLA, azi, dist) = geod_directs(em.geod, lla1, azi, dist)
end

# ╔═╡ e11d2f06-419b-437a-9446-3ca5ccb2b3f1
function lla2tp(lla::LLA)
	(π/2 - lla.lat, lla.lon)
end

# ╔═╡ c2389f7d-072c-4c5c-bd08-35d5b9018dd7
function ll2tp(ll)
	SA_F64[π/2 - ll[1], ll[2]]
end

# ╔═╡ 3d71be7d-ba32-4be0-a1e2-a15a4a5e3a42
function uvdist(ll1, ll2)
	tp1 = ll2tp(ll1)
	tp2 = ll2tp(ll2)
	tp2uv = UVfromThetaPhi()
	uv1 = tp2uv(tp1)
	uv2 = tp2uv(tp2)
	Δuv = norm(uv2 - uv1)
	acos(1 - Δuv)
end

# ╔═╡ 3b710626-7a35-4300-bb34-8b492f89f996
let
	ll1 = deg2rad.((30, 40))
	ll2 = deg2rad.((-30, 120))
	b = uvdist(ll1, ll2)
	a = haversine(ll1, ll2)
	# b = @benchmark SphericalAngle()(reverse($ll1), reverse($ll2))
	a,b
end

# ╔═╡ 86bbf735-d5cb-43c1-8169-fd42bae98809
md"""
# Current Implementation
"""

# ╔═╡ 762fc76f-a3cd-4e99-ac0f-4d61604450f1
md"""
## Add Angular Offset
"""

# ╔═╡ b6d7ac5b-09dd-4a1a-8ff3-d7dfdc290552
todeg(x) = @. $Tuple(rad2deg(x) * °)

# ╔═╡ 21ace099-dbda-482d-bc23-ec253f0a44d2
function test_φ(x, y)
	x̂ = rem2pi(to_radians(x), RoundNearest)
	ŷ = rem2pi(to_radians(y), RoundNearest)
	result = abs(x̂) ≈ abs(ŷ) ≈ π || Base.isapprox(x̂, ŷ; atol=1e-10, rtol=1e-5)
	result || @info "Phi" x y x̂ ŷ Base.isapprox(x̂, ŷ; atol=1e-10, rtol=1e-5)
	result
end

# ╔═╡ d32e2acb-1d68-4d61-a879-55ae60a5b0aa
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

# ╔═╡ 439f3e71-c5f2-4138-81b5-639cd9cb795d
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

# ╔═╡ 6d9a68a0-7b22-4da3-8b03-86e2fc77039f
begin
	test_rot(θ, φ) = RotZYZ(π + φ, -θ, π - φ)
	test_rot(tp) = test_rot(tp...)
end

# ╔═╡ e3bac846-0345-45b5-9ca4-764c54cceb97
#=╠═╡
let
	tp = (10°, 35°)
	a = @benchmark angle_offset_rotation($tp)
	b = @benchmark test_rot($tp)
	a,b
end
  ╠═╡ =#

# ╔═╡ 51419886-08a2-43f4-a08a-06ecbd98ad8f
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

# ╔═╡ 60aefed4-2329-4c04-9475-c3741df59dfb
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

# ╔═╡ 154ae2ee-6e76-4aa5-9113-beb08db4eb03
#=╠═╡
let
	x = (30°, 90°)
	@benchmark angle_offset_rotation($x)
end
  ╠═╡ =#

# ╔═╡ 72388055-2889-4ec8-9748-1ca73314917e
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

# ╔═╡ d47ecc76-3848-46b2-a41b-4db58252c58f
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

# ╔═╡ fe5df978-1321-4abe-b957-b5ac8947a144
begin
function _add_angular_offset2(p₀, offset_angles, input_type::PointingType, output_type::PointingType)
	θ, φ = to_radians(offset_angles)
	sθ, cθ = sincos(θ)
	sφ, cφ = sincos(φ)
	
	# θφ_in = if input_type isa ThetaPhi
	# 	to_radians(p₀)
	# else
	# 	ThetaPhifromUV()(p₀)
	# end
	R = test_rot(p₀)
	perturbation = SA_F64[sθ*cφ, sθ*sφ, cθ]
	p3_out = R * perturbation
	# u,v,w = p3_out
	# @assert w >= 0 "The resulting point has a θ > π/2, so it is located behind the viewer"
	# out = if output_type isa ThetaPhi
	# 	ThetaPhifromUV()(u,v)
	# else
	# 	SVector{2,Float64}(u,v)
	# end
end
end

# ╔═╡ 6158c9ec-f643-4e31-b1ab-17a5c160f09d
begin
function _add_angular_offset3(p₀, offset_angles, input_type::PointingType, output_type::PointingType)
	θ, φ = to_radians(offset_angles)
	sθ, cθ = sincos(θ)
	sφ, cφ = sincos(φ)
	
	# θφ_in = if input_type isa ThetaPhi
	# 	to_radians(p₀)
	# else
	# 	ThetaPhifromUV()(p₀)
	# end
	R = angle_offset_rotation(p₀)
	perturbation = SA_F64[sθ*cφ, sθ*sφ, cθ]
	p3_out = R * perturbation
	# u,v,w = p3_out
	# @assert w >= 0 "The resulting point has a θ > π/2, so it is located behind the viewer"
	# out = if output_type isa ThetaPhi
	# 	ThetaPhifromUV()(u,v)
	# else
	# 	SVector{2,Float64}(u,v)
	# end
end
end

# ╔═╡ 2859fae7-9c53-4ddc-8c56-1343acde5244
#=╠═╡
@benchmark add_angular_offset((0.5,0), (10°, 0))
  ╠═╡ =#

# ╔═╡ 45e55b44-89b8-487e-8596-15684a3dd7f0
#=╠═╡
let
	u = rand(1000).* .5
	f(x) = _add_angular_offset((x,0.0), (10°, 0°), UV(), ThetaPhi())
	@benchmark map($f, $u)
end
  ╠═╡ =#

# ╔═╡ fbba03f0-3699-4e85-bfd0-7faf43f1f426
#=╠═╡
let
	u = rand(1000).* .5
	f(x) = _add_angular_offset2((x,0.0), (10°, 0°), ThetaPhi(), ThetaPhi())
	@benchmark map($f, $u)
end
  ╠═╡ =#

# ╔═╡ ff0f331c-ea50-4041-addb-9da441bb4d36
#=╠═╡
let
	u = [ (rand()* .5, 0) for _ in 1:1000]
	f(x) = _add_angular_offset3(x, (10°, 0°), ThetaPhi(), ThetaPhi())
	f2(x) = _add_angular_offset2(x, (10°, 0°), ThetaPhi(), ThetaPhi())
	@benchmark(map($f, $u)), @benchmark(map($f2, $u))
end
  ╠═╡ =#

# ╔═╡ d56f1852-fa71-451e-ae55-fca58bda9573
md"""
## Get Angular Offset
"""

# ╔═╡ 985726cc-a480-4876-9928-c9c6722854aa
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

# ╔═╡ 485342ef-9172-4a12-899c-79dd3bdfb3e7
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

# ╔═╡ 8a32bc7c-d406-43a3-89a6-60f9c6c1443a
let
	  for i in 1:100
        tp1 = SVector{2}(rand()*90°, (rand()-.5)*360°) |> Tuple
        tp2 = SVector{2}(rand()*90°, (rand()-.5)*360°) |> Tuple
        offset = get_angular_offset(tp1, tp2; input_type=:thetaphi, output_type=:thetaphi) |> todeg
        tp_target = add_angular_offset(tp1, offset; input_type = :thetaphi, output_type = :thetaphi) |> todeg
        (tp_target[1] ≈ tp2[1] && test_φ(tp_target[2], tp2[2])) || error("The forward-reverse offset test with angles failed with $((;tp1, tp2, tp_target, offset))")
	  end
  end

# ╔═╡ 3033a860-0799-4a9d-913a-75bc9eaa01f6
#=╠═╡
let
	x = rand(1000)
	f(x) = get_angular_offset((x,0),(.4,0); input_type = :thetaphi, output_type = :thetaphi)
	@benchmark map($f, $x)
end
  ╠═╡ =#

# ╔═╡ 2d9d6505-d353-417b-9bf6-62d618ac6909
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

# ╔═╡ ea91f87d-45ad-496d-bc8a-cc951f24cf70
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

# ╔═╡ 4076978f-05fb-455f-b9f4-ae965e408c15
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

# ╔═╡ 91c5be8d-6452-4087-8df2-3288bb41cba0
#=╠═╡
let
	x = rand(1000)
	f(x) = get_angular_distance((x,0),(.4,0); input_type = :thetaphi, output_type = :thetaphi)
	@benchmark map($f, $x)
end
  ╠═╡ =#

# ╔═╡ eb2fe2e9-08c6-48bc-a7b9-4cdf148cb796
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

# ╔═╡ 5d1729a5-7f90-4cea-a981-81e74a243c62
md"""
# Tests
"""

# ╔═╡ 55d8b526-7f23-48f2-bdcf-2dea869b7be0
let
		# tp1 = SVector{2}(rand()*50°, (rand()-.5)*360°)
  #       tp2 = SVector{2}(rand()*50°, (rand()-.5)*360°)
		tp1 = SVector{2}(80°, 1°)
        tp2 = SVector{2}(80°, 179°)
        offset = get_angular_offset(tp1, tp2; input_type=:thetaphi, output_type=:thetaphi)
        dist = get_angular_distance(tp1, tp2; input_type=:thetaphi, output_type=:thetaphi)
        offset[1] ≈ dist || error("The distance-offset test in ThetaPhi failed with $((;tp1, tp2, dist, offset))")
end

# ╔═╡ c1d951c0-dd0c-4cbb-8c1f-6206839a47b6
dist_and_bearing(LLA(10°, 1°, 0), LLA(10°, 179°, 0))

# ╔═╡ 203ddd72-48cc-4ec7-b81d-c579ac6a6c37
md"""
# Random
"""

# ╔═╡ 38d173bc-98c5-4c21-be5b-d6c9f3c51f8d
Base.@kwdef struct Constellation
	a::Float64 = 0
	b::Float64 = 1.3
	c::Float64 = 3.2
end

# ╔═╡ 5b0cb831-fc61-43c7-85dc-c1078522064a
Constellation(; b = 5)

# ╔═╡ 6520aa18-c236-446c-b5bd-74abf2811387
Base.@kwdef struct Parameters
	constellation::Constellation = Constellation()
	theta::Float64 = 10
	phi::Float64 = 2
end

# ╔═╡ 6aaeb1ea-94f0-4150-933d-f688fd950d14
function main_function(p::Parameters = Parameters())
	(;phi) = p
	phi + 1
end

# ╔═╡ 26b4c5a6-9b57-4bf0-951d-9f8159694191
main_function()

# ╔═╡ be065ba1-13c6-4893-9b2b-f109d794914e
main_function(Parameters(;phi = 12))

# ╔═╡ b300a251-dad4-45e4-a9b0-3bb665cee3bf
let
	constellation = Constellation(c = 13)
	Parameters(;constellation)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
Distances = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
PlutoDevMacros = "a0499f29-c39b-4c5c-807c-88074221b949"
PlutoExtras = "ed5d0301-4775-4676-b788-cf71e66ff8ed"
PlutoTest = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
BenchmarkTools = "~1.3.2"
Distances = "~0.10.8"
PlutoDevMacros = "~0.5.4"
PlutoExtras = "~0.7.4"
PlutoTest = "~0.2.2"
PlutoUI = "~0.7.51"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.0"
manifest_format = "2.0"
project_hash = "828e055dc57006587cdc2ce0adea2ddd476c8a17"

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

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.2+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Distances]]
deps = ["LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "49eba9ad9f7ead780bfb7ee319f962c811c6d3b2"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.8"

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

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

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

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7302075e5e06da7d000d9bfa055013e3e85578ca"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.9"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.0"

[[deps.PlutoDevMacros]]
deps = ["HypertextLiteral", "InteractiveUtils", "MacroTools", "Markdown", "Pkg", "Random", "TOML"]
git-tree-sha1 = "b4b4a7161e858ad92ffed85753c28284553a54e7"
uuid = "a0499f29-c39b-4c5c-807c-88074221b949"
version = "0.5.4"

[[deps.PlutoExtras]]
deps = ["AbstractPlutoDingetjes", "HypertextLiteral", "InteractiveUtils", "Markdown", "OrderedCollections", "PlutoDevMacros", "PlutoUI", "REPL"]
git-tree-sha1 = "15e75e48e51416d33bab70943923a62a0b63f137"
uuid = "ed5d0301-4775-4676-b788-cf71e66ff8ed"
version = "0.7.4"

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

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "45a7769a04a3cf80da1c1c7c60caf932e6f4c9f7"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.6.0"

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
# ╟─592722de-5f30-416a-9de5-24c8d417c532
# ╠═0d343b68-f73e-11ed-17ef-7d0ceae106f8
# ╠═bca459c7-a29e-4ab7-ab81-ce3005f47a65
# ╟─2d5099cb-a208-433d-8bf4-66f66a8baa28
# ╠═e38b8433-ed1b-4295-86cc-59ddb588867b
# ╠═4edc2a99-3238-4fad-851e-30ca8cba94d1
# ╟─50729705-7e9f-43db-879f-c2f375c21163
# ╟─acca8d87-6f17-4b3c-9efb-5bd84dcd7432
# ╠═9a03ffa9-af46-4e9b-b80d-69b4712490a1
# ╠═b4f84a0b-bcd0-4ccb-b664-f257e60b4c63
# ╠═f230e0c9-5430-459d-8158-748f55cefb2e
# ╠═3d71be7d-ba32-4be0-a1e2-a15a4a5e3a42
# ╠═3b710626-7a35-4300-bb34-8b492f89f996
# ╠═193516b5-6e61-41c9-9a58-a62d7f0791ae
# ╠═cd27bd1c-101f-49cc-aad4-2b568fbf1cbc
# ╠═314d95df-8e5a-4533-9cb1-578cac4e8e4e
# ╟─81981b91-2c82-478e-bdc9-109c9bbcfb6c
# ╠═57c80a92-72b2-45ab-bd8d-e481cd0c4d6c
# ╠═46f68c32-f646-4eca-ba9a-d21faa1e797c
# ╠═61c00eee-7d16-464e-932a-b58b7c72e4c9
# ╠═724afbbc-90ab-40aa-a3b8-93ec6acac37d
# ╠═9722d87d-e617-47c9-9fb6-4a7a5170a50f
# ╟─28124248-18e2-4542-a893-389dbf416129
# ╠═97c5cf9f-be20-49bf-a33e-0a0a2b03ba22
# ╠═d5743ccc-ca6b-4d7c-a6c8-7da0ff52131e
# ╠═e1cf709a-e62a-446d-a661-94d4f808220e
# ╠═3b895f54-5a15-4265-95b1-153af12b91d2
# ╠═52cddfb6-754d-4109-95f7-52795a94297e
# ╠═a1f2dcce-3180-41b6-b5b5-329e7f4fa1d5
# ╠═ba74055a-b73e-4c66-8c3f-3bb3f8e75a38
# ╠═a9a81d0c-e5fa-47e1-95ec-2050f7506df9
# ╠═f4cd828b-92c0-4675-ba44-73a34a674052
# ╠═04d347ab-a061-414f-9fac-7371d4589971
# ╠═57e6eb7f-adce-4e37-b7ef-040254783fea
# ╠═931261d6-7fd6-4639-a0bc-c931542d7c11
# ╠═761027d3-a447-48e9-bb24-85b8d57d6ce8
# ╠═e11d2f06-419b-437a-9446-3ca5ccb2b3f1
# ╠═c2389f7d-072c-4c5c-bd08-35d5b9018dd7
# ╠═86bbf735-d5cb-43c1-8169-fd42bae98809
# ╠═762fc76f-a3cd-4e99-ac0f-4d61604450f1
# ╠═b6d7ac5b-09dd-4a1a-8ff3-d7dfdc290552
# ╠═21ace099-dbda-482d-bc23-ec253f0a44d2
# ╠═d32e2acb-1d68-4d61-a879-55ae60a5b0aa
# ╠═e3bac846-0345-45b5-9ca4-764c54cceb97
# ╠═439f3e71-c5f2-4138-81b5-639cd9cb795d
# ╠═6d9a68a0-7b22-4da3-8b03-86e2fc77039f
# ╠═51419886-08a2-43f4-a08a-06ecbd98ad8f
# ╠═60aefed4-2329-4c04-9475-c3741df59dfb
# ╠═8a32bc7c-d406-43a3-89a6-60f9c6c1443a
# ╠═154ae2ee-6e76-4aa5-9113-beb08db4eb03
# ╠═d47ecc76-3848-46b2-a41b-4db58252c58f
# ╠═72388055-2889-4ec8-9748-1ca73314917e
# ╠═fe5df978-1321-4abe-b957-b5ac8947a144
# ╠═6158c9ec-f643-4e31-b1ab-17a5c160f09d
# ╠═2859fae7-9c53-4ddc-8c56-1343acde5244
# ╠═45e55b44-89b8-487e-8596-15684a3dd7f0
# ╠═fbba03f0-3699-4e85-bfd0-7faf43f1f426
# ╠═ff0f331c-ea50-4041-addb-9da441bb4d36
# ╠═d56f1852-fa71-451e-ae55-fca58bda9573
# ╠═485342ef-9172-4a12-899c-79dd3bdfb3e7
# ╠═4076978f-05fb-455f-b9f4-ae965e408c15
# ╠═985726cc-a480-4876-9928-c9c6722854aa
# ╠═3033a860-0799-4a9d-913a-75bc9eaa01f6
# ╠═91c5be8d-6452-4087-8df2-3288bb41cba0
# ╠═ea91f87d-45ad-496d-bc8a-cc951f24cf70
# ╠═2d9d6505-d353-417b-9bf6-62d618ac6909
# ╠═eb2fe2e9-08c6-48bc-a7b9-4cdf148cb796
# ╟─5d1729a5-7f90-4cea-a981-81e74a243c62
# ╠═55d8b526-7f23-48f2-bdcf-2dea869b7be0
# ╠═c1d951c0-dd0c-4cbb-8c1f-6206839a47b6
# ╟─203ddd72-48cc-4ec7-b81d-c579ac6a6c37
# ╠═38d173bc-98c5-4c21-be5b-d6c9f3c51f8d
# ╠═5b0cb831-fc61-43c7-85dc-c1078522064a
# ╠═6520aa18-c236-446c-b5bd-74abf2811387
# ╠═6aaeb1ea-94f0-4150-933d-f688fd950d14
# ╠═26b4c5a6-9b57-4bf0-951d-9f8159694191
# ╠═be065ba1-13c6-4893-9b2b-f109d794914e
# ╠═b300a251-dad4-45e4-a9b0-3bb665cee3bf
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
