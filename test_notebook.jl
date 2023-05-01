### A Pluto.jl notebook ###
# v0.19.24

using Markdown
using InteractiveUtils

# ╔═╡ 122d494e-326b-11ec-0e29-551537ccbc53
begin
	import Pkg
	Pkg.activate(Base.current_project(@__FILE__))
	using Revise
end

# ╔═╡ 0350a012-0b7d-4f14-89e1-d01b860b2928
begin
	using TelecomUtils
end

# ╔═╡ 07fbb59d-cb6c-4066-878d-130ba1951ee1
p = generate_hex_lattice(1;M=10)

# ╔═╡ 3f56bf85-5d42-4d2e-8d68-b2999ac507f9
col = TelecomUtils.generate_colors(p)

# ╔═╡ 5e26311c-3ffa-45ae-ba66-5728ddd18aa9
UserView(LLA(0,0,100km), EarthModel())

# ╔═╡ Cell order:
# ╠═122d494e-326b-11ec-0e29-551537ccbc53
# ╠═0350a012-0b7d-4f14-89e1-d01b860b2928
# ╠═07fbb59d-cb6c-4066-878d-130ba1951ee1
# ╠═3f56bf85-5d42-4d2e-8d68-b2999ac507f9
# ╠═5e26311c-3ffa-45ae-ba66-5728ddd18aa9
