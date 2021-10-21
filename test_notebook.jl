### A Pluto.jl notebook ###
# v0.16.2

using Markdown
using InteractiveUtils

# ╔═╡ 122d494e-326b-11ec-0e29-551537ccbc53
begin
	import Pkg
	Pkg.activate(".")
end

# ╔═╡ 0350a012-0b7d-4f14-89e1-d01b860b2928
begin
	using Revise
	using PlotlyBase
	using PlutoUtils
	using TelecomUtils
end

# ╔═╡ 07fbb59d-cb6c-4066-878d-130ba1951ee1
p = generate_hex_lattice(1;M=10)

# ╔═╡ 3f56bf85-5d42-4d2e-8d68-b2999ac507f9
col = TelecomUtils.generate_colors(p)

# ╔═╡ 4cf7dfc7-b846-48b4-a902-7c7cb21ee6cd
scatter(p;marker = attr(
	color = TelecomUtils.generate_colors(p,4),
	colorscale = "Jet",
),mode = "markers") |> Plot

# ╔═╡ Cell order:
# ╠═122d494e-326b-11ec-0e29-551537ccbc53
# ╠═0350a012-0b7d-4f14-89e1-d01b860b2928
# ╠═07fbb59d-cb6c-4066-878d-130ba1951ee1
# ╠═3f56bf85-5d42-4d2e-8d68-b2999ac507f9
# ╠═4cf7dfc7-b846-48b4-a902-7c7cb21ee6cd
