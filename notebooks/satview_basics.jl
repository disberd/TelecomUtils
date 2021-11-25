### A Pluto.jl notebook ###
# v0.17.1

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
end

# ╔═╡ e2f50ace-ad14-4e7c-af74-abf3ca9df9fb
md"""
# Packages
"""

# ╔═╡ c3a5e97e-00ed-44af-9c67-fa8198900fbd
#=╠═╡ notebook_exclusive
begin
	using PlutoDevMacros
	using PlutoUtils
	using PlutoTest
end
  ╠═╡ notebook_exclusive =#

# ╔═╡ 82c0d6f9-373d-4866-81eb-9cd2d0981310
#=╠═╡ notebook_exclusive
ToC()
  ╠═╡ notebook_exclusive =#

# ╔═╡ 2e6bcf9b-002f-4fbc-80c1-0cfd2ab1d253
md"""
# Exports
"""

# ╔═╡ 1e3da0c9-2f96-4637-9093-ac7f10c1ad27
md"""
# Helper Functions
"""

# ╔═╡ c4402f72-67ac-4630-a651-da81c1df71bf
"""
	ExtraOutput
Struct used inside `SatView` to use multiple dispatch to give more than the standard output. 
See for example [`get_range`](@ref)
"""
struct ExtraOutput end

# ╔═╡ e796d5be-e011-45d3-ad90-58769feb5e85
md"""
## Show/Print
"""

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
md"""
## Spherical Ellipsoid
"""

# ╔═╡ 855553e3-491b-4e8a-a482-95855697e063
"""
	SphericalEllipsoid(r = 6371e3)
Define a spherical ellipsoid of radius `r` to be used for the various transformation in [`SatView`](@ref). Defaults to 6371km radius
"""
SphericalEllipsoid(r = 6371e3) = Ellipsoid(r,0.0)

# ╔═╡ 6a5cb372-60cb-4ffc-b4f0-22e4016104e7
md"""
# Angle Types
"""

# ╔═╡ e832c1b7-8c04-4146-90f6-1628e91fea2a
const UnitfulAngleType = Union{typeof(°),typeof(rad)}

# ╔═╡ 64cf1b8b-6686-4eeb-a3cc-786300ea7c7d
const UnitfulAngleQuantity = Quantity{<:Real,<:Any,<:UnitfulAngleType}

# ╔═╡ 1d023a0c-a93a-451c-a894-1d1f6a4b78a9
md"""
# SatViewCoordinate types
"""

# ╔═╡ b8ce87d4-4768-4e8a-a16e-9e68b00b6617
abstract type SatViewCoordinate end

# ╔═╡ e3c221c6-6c4a-4b5f-93a6-30d3508ac9d2
md"""
## LLA
"""

# ╔═╡ 8368ae01-ce53-449e-87dd-8dfa3f29f8f4
md"""
Here we want to define a structure that contains useful informations and functions to perform conversions between the view from the satellite based on it's orbital position and points on ground
"""

# ╔═╡ f951805e-515a-475f-893f-bb8b968e425c
md"""
## ERA
"""

# ╔═╡ 86ae20a9-e69c-4d63-9119-395449e9ac09
md"""
ERA stands for Elevation, Range and Azimuth and is used to express the position of a satellite relative to an observer in spherical coordinates.
The elevation is the angle of the pointing with respect to the local horizon of the observer, meaning the plane where the observer is located that is perpendicular to the gravity vector acting on the observe (or in an alternative definition, the plane where the observer is located that is parallel to the tangent plane to the earth ellipsoid at the given lat and lon positions of the observer.
"""

# ╔═╡ 9be2fd5c-4b6c-4e13-b2aa-fb7120a504b7
begin
	"""
	Elevation, Range and Azimuth for a target point on space as seen from a source point on or above the earth surface
	
	# Fields
	- `el::Float64`: Elevation view angle (`0 <= el <= π/2`) between source and target point [rad].
	- `r::Float64`: Range (`r >= 0`) between the source and target points [m].
	- `az::Float64`: Azimuth view angle between source and target point [rad], computed from West to North from the source point perspective. Values provided are automatically converted between -π and π
	
	# Constructors
	
		ERA(el::Real,r::Real,az::Real)
		ERA(el::UnitfulAngleQuantity,r::Real,az::Real)
		ERA(el::UnitfulAngleQuantity,r::Real,az::UnitfulAngleQuantity)
		ERA(el,r::Unitful.Length,az)
	
	where `UnitfulAngleQuantity` is a `Unitful.Quantity` of unit either `u"rad"` or `u"°"`.
	"""
	@with_kw_noshow struct ERA <: SatViewCoordinate
		el::Float64 # Elevation in radians
		r::Float64 # Range in meters
		az::Float64 # Azimuth in radians
		
		function ERA(el::Real,r::Real,az::Real)
			(isnan(el) || isnan(r) || isnan(az)) && return new(NaN,NaN,NaN)  
			@assert el >= 0 && el <= π/2 "Elevation should be between 0 and π/2"
			@assert r >= 0 "Range must be positive"
			new(el,r,rem2pi(az,RoundNearest))
		end
	end

	# Define a constructor that takes combinations of real numbers and angles/lengths
	ERA(el::UnitfulAngleQuantity,r::Real,az::Real) = ERA(
		uconvert(u"rad",el) |> ustrip,
		r,
		az,
	)
	ERA(el::UnitfulAngleQuantity,r::Real,az::UnitfulAngleQuantity) = ERA(
		el,
		r,
		uconvert(u"rad",az) |> ustrip,
	)
	ERA(el,r::Unitful.Length,az) = ERA(
		el,
		uconvert(u"m",r) |> ustrip,
		az,
	)
	
	# Show
	function Base.show(io::IO,era::ERA)
		print(io,"ERA(")
		_print_angle(io,era.el,"el",false)
		_print_length(io,era.r,"r",false)
		_print_angle(io,era.az,"az",true)
		print(io,")")
	end
end

# ╔═╡ 16782c72-ecb1-48ec-8510-78e2e0689a10
function Base.isapprox(x1::ERA, x2::ERA)
	x1.el ≉ x2.el && return false
	# Don't care about different azimuth if elevation is 90°
	x2.el ≉ π/2 && x1.az ≉ x2.az && return false
	x1.r ≉ x2.r && return false
	return true
end

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
		LLA(lat::Real,lon::Real,alt::Real)
		LLA(lat::UnitfulAngleQuantity,lon::Real,alt::Real)
		LLA(lat::UnitfulAngleQuantity,lon::UnitfulAngleQuantity,alt::Real)
		LLA(lat,lon,alt::Unitful.Length)
		LLA(lat,lon) # Defaults to 0.0 altitude
	
	where `UnitfulAngleQuantity` is a `Unitful.Quantity` of unit either `u"rad"` or `u"°"`.
	"""
	@with_kw_noshow struct LLA <: SatViewCoordinate
		lat::Float64 # Latitude in radians
		lon::Float64 # Longitude in radians
		alt::Float64 # Altitude in meters
		
		function LLA(lat::Real,lon::Real,alt::Real)
			(isnan(lat) || isnan(lon) || isnan(alt)) && return new(NaN,NaN,NaN)  
			l2 = rem2pi(lon,RoundNearest)
			@assert abs(lat) <= π/2 "Latitude should be between -π/2 and π/2"
			new(lat,l2,alt)
		end
	end
	
	# Constructor without altitude, assume it is 0
	LLA(lat,lon) = LLA(lat,lon,0.0)
	
	# Define a constructor that takes combinations of real numbers and angles/lengths
	LLA(lat::UnitfulAngleQuantity,lon::Real,alt::Real) = LLA(
		uconvert(u"rad",lat) |> ustrip,
		lon,
		alt)
	LLA(lat::UnitfulAngleQuantity,lon::UnitfulAngleQuantity,alt::Real) = LLA(
		lat,
		uconvert(u"rad",lon) |> ustrip,
		alt)
	LLA(lat,lon,alt::Unitful.Length) = LLA(
		lat,
		lon,
		uconvert(u"m",alt) |> ustrip)	
end

# ╔═╡ 4c06d21c-ac14-4522-bf25-2e0a1ed2d6b9
begin
	export Ellipsoid, SphericalEllipsoid
	export UnitfulAngleQuantity, UnitfulAngleType
	export LLA, ERA
end

# ╔═╡ 3033d5c1-d8e0-4d46-a4b9-7dec4ff7afbd
function Base.isapprox(x1::LLA, x2::LLA)
	x1.alt ≉ x2.alt && return false
	# Don't care about different longitude if latitude is ±90°
	abs(x1.lat) ≈ π/2 && abs(x2.lat) ≈ π/2 && return true
	# Return true if all the lat and lon are matching
	x1.lat ≈ x2.lat && (x1.lon ≈ x2.lon || abs(x1.lon) ≈ abs(x2.lon) ≈ π) && return true
	return false
end

# ╔═╡ 528beffe-0707-4661-8970-def1b1e00ea5
Base.isnan(lla::LLA) = isnan(lla.lat)

# ╔═╡ 11938cb6-46b3-499b-96c0-ef6424d1d0db
#=╠═╡ notebook_exclusive
md"""
# Tests
"""
  ╠═╡ notebook_exclusive =#

# ╔═╡ d2c248b1-c48e-437b-a910-edcc59b4424f
#=╠═╡ notebook_exclusive
md"""
## LLA
"""
  ╠═╡ notebook_exclusive =#

# ╔═╡ 7b306ed5-4bda-465d-abf2-4d07cb4642c1
#=╠═╡ notebook_exclusive
@test $LLA(10°,10°,1000) ≈ $LLA((10+100*eps())*°,10°,1000)
  ╠═╡ notebook_exclusive =#

# ╔═╡ c45f6aac-bff3-4c98-8bd5-c91a98c9eef7
#=╠═╡ notebook_exclusive
@test $LLA(90°,10°,1000) ≈ $LLA(90°,130°,1000)
  ╠═╡ notebook_exclusive =#

# ╔═╡ 1f1505a6-c6af-4fdd-9583-6e783e76de4f
#=╠═╡ notebook_exclusive
@test $LLA(40°,-180°,1000) ≈ $LLA(40°,180°,1000)
  ╠═╡ notebook_exclusive =#

# ╔═╡ 980e48bd-9dd1-4195-8086-40785a7f43e1
#=╠═╡ notebook_exclusive
@test $LLA(10°,10°,1000) !== $LLA((10+100*eps())*°,10°,1000)
  ╠═╡ notebook_exclusive =#

# ╔═╡ fc146750-46b2-4084-a6d2-0d91c0e104e6
#=╠═╡ notebook_exclusive
LLA(1,1,NaN) |> isnan
  ╠═╡ notebook_exclusive =#

# ╔═╡ fc3e96a6-7af8-4ef0-a4ca-b33509aa512f
#=╠═╡ notebook_exclusive
md"""
## ERA
"""
  ╠═╡ notebook_exclusive =#

# ╔═╡ ed2f2fa3-be19-4e2d-ad28-36cb053ed6bd
#=╠═╡ notebook_exclusive
@test ERA(10°,1000,20°) == ERA(10°,1km,deg2rad(20)*rad)
  ╠═╡ notebook_exclusive =#

# ╔═╡ 6151e4d1-a7a6-4fa7-bc77-0e7f3b2a3cc0
#=╠═╡ notebook_exclusive
@test ERA(90°,1000,20°) ≈ ERA(90°,1km,deg2rad(90)*rad)
  ╠═╡ notebook_exclusive =#

# ╔═╡ b5449c6a-9c43-4f9c-81f7-f01277acb109
#=╠═╡ notebook_exclusive
ERA(1,1,NaN) |> isnan
  ╠═╡ notebook_exclusive =#

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Parameters = "d96e819e-fc66-5662-9728-84c9c7592b0a"
SatelliteToolbox = "6ac157d9-b43d-51bb-8fab-48bf53814f4a"
StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[compat]
Parameters = "~0.12.3"
SatelliteToolbox = "~0.9.3"
StaticArrays = "~1.2.13"
Unitful = "~1.9.2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.0-rc2"
manifest_format = "2.0"

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

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "f885e7e7c124f8c92650d61b9477b9ac2ee607dd"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.11.1"

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

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "2db648b6712831ecb333eae76dbfd1c156ca13bb"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.11.2"

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

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

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

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

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

[[deps.OptionalData]]
git-tree-sha1 = "d047cc114023e12292533bb822b45c23cb51d310"
uuid = "fbd9d27c-2d1c-5c1c-99f2-7497d746985d"
version = "1.0.0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PolynomialRoots]]
git-tree-sha1 = "5f807b5345093487f733e520a1b7395ee9324825"
uuid = "3a141323-8675-5d76-9d11-e1df1406c778"
version = "1.0.0"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "d940010be611ee9d67064fe559edbb305f8cc0eb"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.2.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

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

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

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
# ╟─e2f50ace-ad14-4e7c-af74-abf3ca9df9fb
# ╠═35da0c45-6d91-4ad9-811f-5d633684208e
# ╠═c3a5e97e-00ed-44af-9c67-fa8198900fbd
# ╠═82c0d6f9-373d-4866-81eb-9cd2d0981310
# ╟─2e6bcf9b-002f-4fbc-80c1-0cfd2ab1d253
# ╠═4c06d21c-ac14-4522-bf25-2e0a1ed2d6b9
# ╟─1e3da0c9-2f96-4637-9093-ac7f10c1ad27
# ╠═c4402f72-67ac-4630-a651-da81c1df71bf
# ╟─e796d5be-e011-45d3-ad90-58769feb5e85
# ╠═d890aff1-dbd0-451c-bf14-bde9758c3be0
# ╠═fdbbc8d9-83d6-413e-aff8-aca18f24dfea
# ╟─cc0cae75-ba10-4a62-b0ef-22259e40a083
# ╠═855553e3-491b-4e8a-a482-95855697e063
# ╟─6a5cb372-60cb-4ffc-b4f0-22e4016104e7
# ╠═e832c1b7-8c04-4146-90f6-1628e91fea2a
# ╠═64cf1b8b-6686-4eeb-a3cc-786300ea7c7d
# ╟─1d023a0c-a93a-451c-a894-1d1f6a4b78a9
# ╠═b8ce87d4-4768-4e8a-a16e-9e68b00b6617
# ╟─e3c221c6-6c4a-4b5f-93a6-30d3508ac9d2
# ╟─8368ae01-ce53-449e-87dd-8dfa3f29f8f4
# ╠═f207d849-ebff-4e6c-95bb-50693cb7c9b6
# ╠═3033d5c1-d8e0-4d46-a4b9-7dec4ff7afbd
# ╠═528beffe-0707-4661-8970-def1b1e00ea5
# ╟─f951805e-515a-475f-893f-bb8b968e425c
# ╟─86ae20a9-e69c-4d63-9119-395449e9ac09
# ╠═9be2fd5c-4b6c-4e13-b2aa-fb7120a504b7
# ╠═16782c72-ecb1-48ec-8510-78e2e0689a10
# ╠═7344190c-7989-4b55-b7be-357f7d6b7370
# ╟─11938cb6-46b3-499b-96c0-ef6424d1d0db
# ╟─d2c248b1-c48e-437b-a910-edcc59b4424f
# ╠═7b306ed5-4bda-465d-abf2-4d07cb4642c1
# ╠═c45f6aac-bff3-4c98-8bd5-c91a98c9eef7
# ╠═1f1505a6-c6af-4fdd-9583-6e783e76de4f
# ╠═980e48bd-9dd1-4195-8086-40785a7f43e1
# ╠═fc146750-46b2-4084-a6d2-0d91c0e104e6
# ╟─fc3e96a6-7af8-4ef0-a4ca-b33509aa512f
# ╠═ed2f2fa3-be19-4e2d-ad28-36cb053ed6bd
# ╠═6151e4d1-a7a6-4fa7-bc77-0e7f3b2a3cc0
# ╠═b5449c6a-9c43-4f9c-81f7-f01277acb109
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
