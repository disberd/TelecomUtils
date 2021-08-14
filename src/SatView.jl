### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 590cdbce-fc45-11eb-2fde-1d27628251b7
begin
	using PlutoUtils
	using Proj4
	using DocStringExtensions
	using CoordinateTransformations
	using StaticArrays
	using LinearAlgebra
	using Unitful
	using Unitful.DefaultSymbols
end

# ╔═╡ bcee47c1-81c1-4c77-af09-22c04374fa4f
using Parameters

# ╔═╡ 74422a23-0760-470f-9e1e-43b8c3972f65
hide_cell_shortcut()

# ╔═╡ 77e399b7-0f7e-4ff1-9f8e-fd0f3408e894
ToC()

# ╔═╡ ac5e6327-98a1-478f-be65-05fa1cff717d
md"""
# Functions
"""

# ╔═╡ c023c6f2-ebcc-4bf2-8898-2563bf97ca45
md"""
## Define the angle types 
"""

# ╔═╡ c0a080c2-6e70-45e0-8482-25deb97da3e3
const angle_types = Union{typeof(°),typeof(rad)}

# ╔═╡ 8171a163-006d-4f80-bd2d-4dd05088055c
const angle_quantity_type = Quantity{<:Real,<:Any,<:angle_types}

# ╔═╡ f68e9146-cdac-42dd-afe3-e8bc8c2211d6
md"""
## Transformation betwee degrees and real
"""

# ╔═╡ 189afaa6-1d94-4f34-b44f-66994d728f58
md"""
To be able to interact the transformation coming from Proj4 (which only accept real number) with potentially inputs containing angular unit data, we need to create an artificial transformation that takes whatever vector,tuple or number of angular inputs, and outputs an SVector that contains the same number but converted to degrees and stripped from their unit
"""

# ╔═╡ f3e17b34-142e-41b7-946a-be800283c4e7
begin
	struct Degree2Real <: Transformation end
	struct Real2Degree <: Transformation end
	
	Base.inv(t::Degree2Real) = Real2Degree()
	Base.inv(t::Real2Degree) = Degree2Real()
	
	# Assume that real numbers are already in degrees
	(t::Degree2Real)(d::Real...) = SVector(d...)
	(t::Degree2Real)(d::angle_quantity_type...) = SVector(@. uconvert(u"°",d) |> ustrip)
	(t::Degree2Real)(d::SVector{<:Any,<:angle_quantity_type}) = t(d...)
	(t::Degree2Real)(d::Tuple) = t(d...)
	
	# Do the opposite, accept real numbers and give back Degrees
	(t::Real2Degree)(r::Real...) = SVector(r .* °)
	(t::Real2Degree)(r::SVector{<:Any,<:Real}) = t(r...)
	(t::Real2Degree)(r::Tuple) = t(r...)
end

# ╔═╡ f65769b4-04f5-457b-beb6-5865c97108d4
#=╠═╡ notebook_exclusive
md"""
### Tests
"""
  ╠═╡ notebook_exclusive =#

# ╔═╡ f4c1d490-81f5-4243-8592-276020a80d39
#=╠═╡ notebook_exclusive
Degree2Real()(10,20)
  ╠═╡ notebook_exclusive =#

# ╔═╡ cb73c175-0ae7-488e-9f9c-6460d1b08a11
#=╠═╡ notebook_exclusive
Degree2Real()((10°,20°))
  ╠═╡ notebook_exclusive =#

# ╔═╡ 4eaf7448-3c1f-47b8-a9eb-e8b5aca6693b
#=╠═╡ notebook_exclusive
Degree2Real()(SVector(1rad,.2rad))
  ╠═╡ notebook_exclusive =#

# ╔═╡ fa0b1a6f-44ce-4fb8-9e89-673913b9234d
#=╠═╡ notebook_exclusive
Degree2Real()((10°,20°)) |> inv(Degree2Real())
  ╠═╡ notebook_exclusive =#

# ╔═╡ 8cf22dae-6e53-469f-aa75-f7de1cc79ee4
md"""
## Near-sided perspective to UV
"""

# ╔═╡ 1c9de582-b9ae-4959-ad94-97a9db5802b7
md"""
The outcome of the near-sided perspective projection from the PROJ library provides the coordinates on the projection plane that is tangent to the earth at the sub-satellite point.

To be converted into uv, this coordinates have to be scaled depending on the altitude of the satellite that generated the projection.
"""

# ╔═╡ 6c8dc4fe-0ca2-4038-b1fc-6e00747b3ad8


# ╔═╡ 434ef478-3298-4ea7-b8cb-161181abdb2a
md"""
# Satellite View
"""

# ╔═╡ f0758e99-9f2b-4934-88eb-7e62cdd5c51f
md"""
Here we want to define a structure that contains useful informations and functions to perform conversions between the view from the satellite based on it's orbital position and points on ground
"""

# ╔═╡ 41599efd-45b8-471c-a8ec-2fde36b4f58f
begin
	@with_kw_noshow struct LLA
		lat::typeof(1.0u"°")
		lon::typeof(1.0u"°")
		alt::typeof(1.0u"km")
	end
	
	# Define a constructor that takes combinations of real numbers, with lat/lon defaulting to degrees
	LLA(lat::Real,lon::Real,alt::Real) = LLA(lat*°,lon*°,alt*m)
	LLA(lat::Real,lon::Real,alt::Unitful.Length) = LLA(lat*°,lon*°,alt)
end

# ╔═╡ 0ddb5932-89a2-4c95-ab6b-fa73d866e3af
#=╠═╡ notebook_exclusive
LLA(10,10,100)
  ╠═╡ notebook_exclusive =#

# ╔═╡ e72d6e9e-8234-4090-8c8d-187ff5bce5b8
@with_kw_noshow struct SatView
	lla::LLA
end

# ╔═╡ 3b8ce9f3-137b-46a1-81d5-4334e81df27e
SatView(LLA(1,1,100km))

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CoordinateTransformations = "150eb455-5306-5404-9cee-2592286d6298"
DocStringExtensions = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Parameters = "d96e819e-fc66-5662-9728-84c9c7592b0a"
PlutoUtils = "ed5d0301-4775-4676-b788-cf71e66ff8ed"
Proj4 = "9a7e659c-8ee8-5706-894e-f68f43bc57ea"
StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[compat]
CoordinateTransformations = "~0.6.1"
DocStringExtensions = "~0.8.5"
Parameters = "~0.12.2"
PlutoUtils = "~0.3.4"
Proj4 = "~0.7.5"
StaticArrays = "~1.2.12"
Unitful = "~1.9.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.0-beta2"
manifest_format = "2.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.CEnum]]
git-tree-sha1 = "215a9aa4a1f23fbd05b92769fdd62559488d70e9"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.1"

[[deps.Chain]]
git-tree-sha1 = "cac464e71767e8a04ceee82a889ca56502795705"
uuid = "8be319e6-bccf-4806-a6f7-6fae938471bc"
version = "0.4.8"

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
git-tree-sha1 = "6d1c23e740a586955645500bbec662476204a52c"
uuid = "150eb455-5306-5404-9cee-2592286d6298"
version = "0.6.1"

[[deps.Crayons]]
git-tree-sha1 = "3f71217b538d7aaee0b69ab47d9b7724ca8afa0d"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.0.4"

[[deps.DataAPI]]
git-tree-sha1 = "ee400abb2298bd13bfc3df1c412ed228061a2385"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.7.0"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "a32185f5428d3986f47c2ab78b1f216d5e6cc96f"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.5"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.Glob]]
git-tree-sha1 = "4df9f7e06108728ebf00a0a11edee4b29a482bb2"
uuid = "c27321d9-0574-5035-807b-f59d2c89b15c"
version = "1.3.0"

[[deps.HypertextLiteral]]
git-tree-sha1 = "1e3ccdc7a6f7b577623028e0095479f4727d8ec1"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.8.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

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

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

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
git-tree-sha1 = "2276ac65f1e236e0a6ea70baff3f62ad4c625345"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.2"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "477bf42b4d1496b454c10cce46645bb5b8a0cf2c"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.0.2"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlutoTest]]
deps = ["HypertextLiteral", "InteractiveUtils", "Markdown", "Test"]
git-tree-sha1 = "3479836b31a31c29a7bac1f09d95f9c843ce1ade"
uuid = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
version = "0.1.0"

[[deps.PlutoUI]]
deps = ["Base64", "Dates", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "Suppressor"]
git-tree-sha1 = "44e225d5837e2a2345e69a1d1e01ac2443ff9fcb"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.9"

[[deps.PlutoUtils]]
deps = ["Chain", "Glob", "HypertextLiteral", "InteractiveUtils", "Markdown", "PlutoTest", "PlutoUI", "PrettyTables", "Reexport", "Requires", "UUIDs"]
git-tree-sha1 = "db3eaef2cc68f99bb41a8600f882e016f718f65a"
uuid = "ed5d0301-4775-4676-b788-cf71e66ff8ed"
version = "0.3.4"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00cfd92944ca9c760982747e9a1d0d5d86ab1e5a"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.2"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "0d1245a357cc61c8cd61934c07447aa569ff22e6"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.1.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Proj4]]
deps = ["CEnum", "CoordinateTransformations", "PROJ_jll", "StaticArrays"]
git-tree-sha1 = "a26cc8a99169e41c7d60719d4bddf2ce7adb4069"
uuid = "9a7e659c-8ee8-5706-894e-f68f43bc57ea"
version = "0.7.5"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "5f6c21241f0f655da3952fd60aa18477cf96c220"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.1.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "4036a3bd08ac7e968e27c203d45f5fff15020621"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.1.3"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.SQLite_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "9a0e24b81e3ce02c4b2eb855476467c7b93b8a8f"
uuid = "76ed43ae-9a5d-5a62-8c75-30186b810ce8"
version = "3.36.0+0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "3240808c6d463ac46f1c1cd7638375cd22abbccb"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.2.12"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.Suppressor]]
git-tree-sha1 = "a819d77f31f83e5792a76081eee1ea6342ab8787"
uuid = "fd094767-a336-5f1f-9728-57cf17d0bbfb"
version = "0.2.0"

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
git-tree-sha1 = "d0c690d37c73aeb5ca063056283fde5585a41710"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.5.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

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
git-tree-sha1 = "a981a8ef8714cba2fd9780b22fd7a469e7aaf56d"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.9.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "cc4bf3fdde8b7e3e9fa0351bdeedba1cf3b7f6e6"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.0+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll", "Pkg"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╠═590cdbce-fc45-11eb-2fde-1d27628251b7
# ╠═bcee47c1-81c1-4c77-af09-22c04374fa4f
# ╠═74422a23-0760-470f-9e1e-43b8c3972f65
# ╠═77e399b7-0f7e-4ff1-9f8e-fd0f3408e894
# ╟─ac5e6327-98a1-478f-be65-05fa1cff717d
# ╟─c023c6f2-ebcc-4bf2-8898-2563bf97ca45
# ╠═c0a080c2-6e70-45e0-8482-25deb97da3e3
# ╠═8171a163-006d-4f80-bd2d-4dd05088055c
# ╟─f68e9146-cdac-42dd-afe3-e8bc8c2211d6
# ╟─189afaa6-1d94-4f34-b44f-66994d728f58
# ╠═f3e17b34-142e-41b7-946a-be800283c4e7
# ╟─f65769b4-04f5-457b-beb6-5865c97108d4
# ╠═f4c1d490-81f5-4243-8592-276020a80d39
# ╠═cb73c175-0ae7-488e-9f9c-6460d1b08a11
# ╠═4eaf7448-3c1f-47b8-a9eb-e8b5aca6693b
# ╠═fa0b1a6f-44ce-4fb8-9e89-673913b9234d
# ╟─8cf22dae-6e53-469f-aa75-f7de1cc79ee4
# ╟─1c9de582-b9ae-4959-ad94-97a9db5802b7
# ╠═6c8dc4fe-0ca2-4038-b1fc-6e00747b3ad8
# ╟─434ef478-3298-4ea7-b8cb-161181abdb2a
# ╟─f0758e99-9f2b-4934-88eb-7e62cdd5c51f
# ╠═41599efd-45b8-471c-a8ec-2fde36b4f58f
# ╠═0ddb5932-89a2-4c95-ab6b-fa73d866e3af
# ╠═e72d6e9e-8234-4090-8c8d-187ff5bce5b8
# ╠═3b8ce9f3-137b-46a1-81d5-4334e81df27e
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
