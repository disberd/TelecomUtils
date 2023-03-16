### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ d2056fd0-5419-11ec-2226-0b24c663de65
begin
	using Interpolations
end

# ╔═╡ 2c06b5c1-0345-444f-9664-338a60fedc83
# ╠═╡ skip_as_script = true
#=╠═╡
begin
	using BenchmarkTools
	using PlutoUtils
	using PlotlyBase
end
  ╠═╡ =#

# ╔═╡ 1c0cf7ad-9579-4593-ae1f-441788c049ae
# ╠═╡ skip_as_script = true
#=╠═╡
ToC()
  ╠═╡ =#

# ╔═╡ 47a21d25-7710-4a93-8cbd-ad93cee68b74
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
# Packages
"""
  ╠═╡ =#

# ╔═╡ 78ad9ad5-7c06-41e1-940d-1aad1de7d5f2
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
# Types
"""
  ╠═╡ =#

# ╔═╡ 35450d8d-cda5-46ba-884e-636abf023ac1
abstract type CapacityComputationType end

# ╔═╡ a1158d4b-38cc-4c9b-a949-e94e98c2b5f6
"""
	DVBS2x <: CapacityComputationType
Structure used to dispatch on `speff2snr` and `snr2speff` to simulate the MODCODs of the DVB-S2x air interface.
# Note
The threshold values are taken from the error performance chapter (6) of the [DVB-S2x  standard](https://dvb.org/?standard=second-generation-framing-structure-channel-coding-and-modulation-systems-for-broadcasting-interactive-services-news-gathering-and-other-broadband-satellite-applications-part-2-dvb-s2-extensions) assuming AWGN conditions

See also: [`Shannon`](@ref), [`NR5G`](@ref), [`speff2snr`](@ref), [`snr2speff`](@ref)
"""
struct DVBS2x <: CapacityComputationType end

# ╔═╡ 0345b4e1-49e6-4fc2-bc0e-6f179a0b1792
"""
	NR5G <: CapacityComputationType
Structure used to dispatch on `speff2snr` and `snr2speff` to simulate the MODCODs of the 5G NR air interface.
# Note
The thresholds for the 5G approximation are taken from a document by Intel that can be found at [this link](https://www.3gpp.org/ftp/TSG_RAN/WG1_RL1/TSGR1_90/Docs/R1-1712553.zip).
The document present a picture of the spectral efficiency vs SNR for various 5G NR MODCODs.

As the tabular data is not provided in the document, the points are extrapolated from the image using [this tool](https://automeris.io/WebPlotDigitizer/). 

The thresholds extrapolated from the image with the tool are increased by 1dB, since the values from the picture refer to 10% BLER (Block Error Rate), and we want to use something more akin to QEF (Quasi Error Free) as done for the DVB-S2x thresholds.

**For this reason, the results for `snr2speff` and `speff2snr` for the 5G NR waveforms will be a linear interpolation based on the extrapolated values, rather than a step function corresponding to existing MODCODs**

See also: [`DVBS2x`](@ref), [`Shannon`](@ref), [`speff2snr`](@ref), [`snr2speff`](@ref)
"""
struct NR5G <: CapacityComputationType end

# ╔═╡ 94032874-a85b-4632-bc24-c8c0154b518a
"""
	Shannon <: CapacityComputationType
Structure used to dispatch on `speff2snr` and `snr2speff` to decide the method to use for the computation

See also: [`DVBS2x`](@ref), [`NR5G`](@ref), [`speff2snr`](@ref), [`snr2speff`](@ref)
"""
struct Shannon <: CapacityComputationType end

# ╔═╡ 4d74427d-21e3-4736-ab29-ef65c7ea0490
Base.broadcastable(c::CapacityComputationType) = Ref(c)

# ╔═╡ a1c5f1c0-adc2-43f0-bd63-e4c6c7eca52b
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
# DVB-S2x Thresholds
"""
  ╠═╡ =#

# ╔═╡ 5907b0b7-f81b-4cc3-a585-8f52fd191bae
# Thresholds 1:5 are for VLSNR, and 4:5 have a spreading factor of 2
# We put -.1 instead of 0 on the first point to avoid errors in case ones provide speff2snr with exactly 0
const dvbs2x_thresholds_linear = [0.102329299,0.147910839,0.245470892,0.323593657,0.424619564,0.5188000389,0.5821032178,0.6266138647,0.751622894,0.9332543008,1.051961874,1.258925412,1.396368361,1.671090614,2.041737945,2.529297996,2.937649652,2.971666032,3.25836701,3.548133892,3.953666201,4.518559444,4.83058802,5.508076964,6.45654229,6.886522963,6.966265141,7.888601176,8.452788452,9.354056741,10.49542429,11.61448614,12.67651866,12.88249552,14.48771854,14.96235656,16.48162392,18.74994508,20.18366364,23.1206479,25.00345362,30.26913428,35.2370871,38.63669771,45.18559444,49.88844875,52.96634439,64.5654229,72.27698036,76.55966069,90.57326009]

# ╔═╡ 37772abb-c600-4beb-a92e-7da3c3021370
# Spectral Efficiencies related to the DVBS2x MODCODs
const dvbs2x_modcod_speff =[0.095000,0.116111,0.184889,0.225975,0.333333,0.434841,0.490243,0.567805,0.656448,0.789412,0.889135,0.988858,1.088581,1.188304,1.322253,1.487473,1.587196,1.647211,1.713601,1.779991,1.972253,2.10485,2.193247,2.370043,2.458441,2.524739,2.635236,2.637201,2.745734,2.856231,2.966728,3.077225,3.165623,3.289502,3.300184,3.510192,3.620536,3.703295,3.841226,3.951571,4.206428,4.338659,4.603122,4.735354,4.933701,5.06569,5.241514,5.417338,5.593162,5.768987,5.900855]

# ╔═╡ 27cda06f-de16-47f6-be79-324963b155c9
const n_MODCODS_DVBS2x = length(dvbs2x_thresholds_linear)

# ╔═╡ f03f7ebb-dad0-4766-8615-afde6478f59d
# Find the minimum SNR needed for a MODCOD that has at least the specified spectral efficiency
const _dvbs2x_itp_speff = LinearInterpolation((dvbs2x_modcod_speff,),1:n_MODCODS_DVBS2x; extrapolation_bc = 1)

# ╔═╡ 1f389615-64d5-4084-b4a7-e7838a19ff9b
# Find the spectral efficiency that can be achieved with DVB-S2x modcods given a certain SNR 
const _dvbs2x_itp_snr = LinearInterpolation((dvbs2x_thresholds_linear, ), 1:n_MODCODS_DVBS2x; extrapolation_bc = n_MODCODS_DVBS2x)

# ╔═╡ eac20e32-59a3-44c3-9c68-31be56ad8ee9
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
# 5G Thresholds
"""
  ╠═╡ =#

# ╔═╡ ff6dd97b-66db-48e8-bd6d-bfccff28ba39
md"""
The thresholds for the 5G approximation are taken from a document by Intel that can be found at [this link](https://www.3gpp.org/ftp/TSG_RAN/WG1_RL1/TSGR1_90/Docs/R1-1712553.zip).
The document present a picture of the spectral efficiency vs SNR for various 5G NR MODCODs.

As the tabular data is not provided in the document, the points are extrapolated from the image using [this tool](https://automeris.io/WebPlotDigitizer/). 

The thresholds extrapolated from the image with the tool are increased by 1dB, since the values from the picture refer to 10% BLER (Block Error Rate), and we want to use something more akin to QEF (Quasi Error Free) as done for the DVB-S2x thresholds
"""

# ╔═╡ e48e3061-4004-4cee-8a93-7218e3ba5bbc
# We add 1 to the SNR as the thresholds were derived for a 10% BLER, we want to make it QEF
const nr5g_thresholds_linear = 1 .+ [-7, -6.154185022026431, -5.429200755191944, -4.583385777218375, -3.737570799244807, -2.8112020138451843, -2.0862177470106973, -1.3612334801762085, -0.6362492133417224, 0.2498426683448738, 1.0151038388923883, 1.9011957205789844, 2.4650723725613606, 2.9886721208307137, 3.75393329137823, 4.559471365638769, 5.244178728760229, 5.848332284455637, 6.61359345500315, 7.419131529263693, 8.506607929515424, 9.513530522341098, 10.157960981749534, 11.003775959723102, 12.010698552548776, 12.775959723096294, 13.50094398993078, 14.427312775330405, 15.152297042164887, 15.917558212712406, 16.60226557583387, 17.488357457520458, 18.45500314663311, 19.3410950283197, 19.98552548772814, 20.509125235997487, 20.952171176840785, 21.677155443675275, 22.402139710509758, 23.00629326620517, 23.73127753303966, 24.49653870358717, 25.10069225928258] |> x -> 10 .^ (x./10)

# ╔═╡ d2fbd3c2-215e-4859-861e-03a0b2efdc9b
const nr5g_modcods_speff = [0.20496894409937916,  0.22981366459627317,  0.279503105590063,  0.341614906832298,  0.44099378881987583,  0.5031055900621126,  0.6024844720496905,  0.6956521739130421,  0.7888198757763982,  0.9254658385093171,  1.0496894409937898,  1.1739130434782608,  1.2608695652173916,  1.3354037267080754,  1.4596273291925472,  1.633540372670808,  1.7950310559006217,  1.9689440993788825,  2.167701863354038,  2.366459627329193,  2.6149068322981375,  2.863354037267081,  3.024844720496895,  3.2360248447204976,  3.4844720496894412,  3.695652173913044,  3.9440993788819885,  4.229813664596274,  4.465838509316771,  4.6645962732919255,  4.850931677018634,  5.099378881987578,  5.360248447204969,  5.658385093167702,  5.919254658385094,  6.093167701863354,  6.242236024844721,  6.490683229813666,  6.714285714285715,  6.91304347826087,  7.099378881987578,  7.2981366459627335,  7.409937888198758]

# ╔═╡ b85c0bff-162d-402e-8d85-44d46b4af8a8
const n_MODCODS_5GNR = length(nr5g_thresholds_linear)

# ╔═╡ c2062a21-8027-4d49-85d6-6666f761cabb
# Find the spectral efficiency that can be achieved with 5G NR given a certain SNR 
const _nr5g_itp_snr = LinearInterpolation((nr5g_thresholds_linear, ), nr5g_modcods_speff; extrapolation_bc = Flat())

# ╔═╡ 29d90f34-0aeb-48d3-80b4-759bb0b77c2f
# Find the minimum SNR needed for a MODCOD that has at least the specified spectral efficiency
const _nr5g_itp_speff = LinearInterpolation((nr5g_modcods_speff,),nr5g_thresholds_linear; extrapolation_bc = Flat())

# ╔═╡ 712f157b-8cc9-4034-8f54-8b57a7953683
# ╠═╡ skip_as_script = true
#=╠═╡
md"""
# Exported Functions
"""
  ╠═╡ =#

# ╔═╡ 4b14629a-b25f-4f3c-94bd-14e77c044695
"""
	speff2snr(DVBS2x(), η)
	speff2snr(Shannon(), η)
Computes the minimum SNR [linear] that is required to achieve the target spectral efficiency `η` [b/s/Hz].

The first argument specifies whether to use the DVBS2x thresholds or the Shannon formula for the computation

See also: [`DVBS2x`](@ref), [`Shannon`](@ref), [`snr2speff`](@ref)
"""
speff2snr(::DVBS2x, η) = η > dvbs2x_modcod_speff[end] ? NaN : dvbs2x_thresholds_linear[ceil(Int,_dvbs2x_itp_speff(η))]

# ╔═╡ 215ca230-77fd-4bbb-a4f7-8af1a0c9bcb3
"""
	snr2speff(DVBS2x(), linear_snr)
	snr2speff(Shannon(), linear_snr)
Computes the minimum the spectral efficiency `η` [b/s/Hz] that can be achieved given the input SNR [linear] `linear_snr`.

The first argument specifies whether to use the DVBS2x thresholds or the Shannon formula for the computation.

See also: [`DVBS2x`](@ref), [`Shannon`](@ref), [`speff2snr`](@ref)
"""
snr2speff(::DVBS2x, snr) = snr < dvbs2x_thresholds_linear[1] ? 0.0 : dvbs2x_modcod_speff[floor(Int,_dvbs2x_itp_snr(snr))]

# ╔═╡ 57500cea-75e5-473a-83c6-635f3ec3fa95
snr2speff(::NR5G, snr) = snr < nr5g_thresholds_linear[1] ? 0.0 : _nr5g_itp_snr(snr)

# ╔═╡ 774fc0f5-a11f-42eb-9bc0-3faa13ed24f4
speff2snr(::NR5G, η) = η > nr5g_modcods_speff[end] ? NaN : _nr5g_itp_speff(η)

# ╔═╡ 5203bbc5-4746-48cf-b209-ec618b223af4
snr2speff(::Shannon, snr_lin) = log2(1 + snr_lin)

# ╔═╡ 154d61e7-9014-4cc0-b817-3257fcfe1826
speff2snr(::Shannon, speff) = 2^(speff) - 1

# ╔═╡ 6c71d117-c2e5-4eb4-8abb-39109ffa30db
export DVBS2x, NR5G, Shannon, snr2speff, speff2snr

# ╔═╡ 602c9baf-96e5-4e2d-9e84-3bf177d137e9
# ╠═╡ skip_as_script = true
#=╠═╡
a = rand(1000,1000)
  ╠═╡ =#

# ╔═╡ 039bfd43-c1e1-4086-b108-2079ce5e0453
# ╠═╡ skip_as_script = true
#=╠═╡
b = sort(vec(a))
  ╠═╡ =#

# ╔═╡ 1696f2e6-1753-41b4-b6c1-ab8038680c91
# ╠═╡ skip_as_script = true
#=╠═╡
@benchmark sort(vec($a))
  ╠═╡ =#

# ╔═╡ 03b7a36a-8a82-4834-baf1-a4114f6afca3
# ╠═╡ skip_as_script = true
#=╠═╡
@benchmark speff2snr.(Shannon(), $a)
  ╠═╡ =#

# ╔═╡ b959608b-2a49-4397-b974-9094480596ab
# ╠═╡ skip_as_script = true
#=╠═╡
@benchmark snr2speff.(Shannon(), $a)
  ╠═╡ =#

# ╔═╡ d0259918-920a-4b8e-b33d-a93da47425de
# ╠═╡ skip_as_script = true
#=╠═╡
@benchmark speff2snr.(DVBS2x(), $a)
  ╠═╡ =#

# ╔═╡ b6909885-61c0-4403-9639-c00ca7cd55e1
# ╠═╡ skip_as_script = true
#=╠═╡
@benchmark snr2speff.(DVBS2x(), $a)
  ╠═╡ =#

# ╔═╡ 353923b3-ab53-4fdf-8eb3-f7a4d126bdf3
# ╠═╡ skip_as_script = true
#=╠═╡
@benchmark speff2snr.(NR5G(), $a)
  ╠═╡ =#

# ╔═╡ 5953946c-52fc-4297-bc3a-9aa2759c8be4
# ╠═╡ skip_as_script = true
#=╠═╡
@benchmark snr2speff.(NR5G(), $a)
  ╠═╡ =#

# ╔═╡ 080fdaa3-b9d7-4030-8bbc-88abcc5d37a3
# ╠═╡ skip_as_script = true
#=╠═╡
@benchmark snr2speff.(DVBS2x(), $b)
  ╠═╡ =#

# ╔═╡ 2f49954a-8e37-4d84-accc-f5537e734331
# ╠═╡ skip_as_script = true
#=╠═╡
snr = range(-6,20,1000)
  ╠═╡ =#

# ╔═╡ 1a6110f9-4086-43ee-8a12-26755ac0ee8a
# ╠═╡ skip_as_script = true
#=╠═╡
snr_lin = 10 .^ (snr./10)
  ╠═╡ =#

# ╔═╡ 87ca96a1-8f5a-4534-9e55-8f4bc4f0e437
# ╠═╡ skip_as_script = true
#=╠═╡
Plot(
	[
		scatter(x = snr, y = snr2speff.(DVBS2x(), snr_lin), name = "DVB-S2X"),
		scatter(x = snr, y = snr2speff.(NR5G(), snr_lin), name = "5G-NR"),
	], Layout(
		xaxis_title = "SNR (dB)", 
		yaxis_title = "Spectral Efficiency",
		legend = attr(
			x = 0,
			y = 1,
			xanchor = "left",
			bgcolor = "rgba(0,0,0,0)",
		),
		width = 500,
		height = 400,
		paper_bgcolor = "rgba(0,0,0,0)",
	), config = PlotConfig(
		toImageButtonOptions = Dict(
			:filename => "dvb_5g_modcods",
			:scale => 2,
		)
	)
)
  ╠═╡ =#

# ╔═╡ 6c13719f-3b65-4550-88e9-6f1c5cd84394
# ╠═╡ skip_as_script = true
#=╠═╡
speff = range(0,7,1000)
  ╠═╡ =#

# ╔═╡ c68bc281-eb52-43db-b216-fa01350b932a
# ╠═╡ skip_as_script = true
#=╠═╡
snr_lin_min = speff2snr.(DVBS2x(), speff)
  ╠═╡ =#

# ╔═╡ 0dcb4a08-b8ed-4a4b-8be2-7ccda056846a
# ╠═╡ skip_as_script = true
#=╠═╡
Plot(
	[
		scatter(;x=speff, y=10log10.(speff2snr.(DVBS2x(), speff)), name = "DVB-S2x"),
		scatter(;x=speff, y=10log10.(speff2snr.(NR5G(), speff)), name = "5G-NR"),
	], Layout(xaxis_title = "Target Spectral Efficiency", yaxis_title = "Minimum SNR (dB)");
)
  ╠═╡ =#

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
Interpolations = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
PlotlyBase = "a03496cd-edff-5a9b-9e67-9cda94a718b5"
PlutoUtils = "ed5d0301-4775-4676-b788-cf71e66ff8ed"

[compat]
BenchmarkTools = "~1.3.1"
Interpolations = "~0.13.4"
PlotlyBase = "~0.8.18"
PlutoUtils = "~0.5.9"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.0-rc1"
manifest_format = "2.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

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
git-tree-sha1 = "4c10eee4af024676200bc7752e536f858c6b8f93"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.3.1"

[[deps.Chain]]
git-tree-sha1 = "8c4920235f6c561e401dfe569beb8b924adad003"
uuid = "8be319e6-bccf-4806-a6f7-6fae938471bc"
version = "0.5.0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "2dd813e5f2f7eec2d1268c57cf2373d3ee91fcea"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.1"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "1e315e3f4b0b7ce40feded39c73049692126cf53"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.3"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "1fd869cc3875b57347f7027521f561cf46d1fcd8"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.19.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "d08c20eef1f2cbc6e60fd3612ac4340b89fea322"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.9"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "924cdca592bc16f14d2f7006754a621735280b74"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.1.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "fb5f5316dd3fd4c5e7c30a24d50643b73e37cd40"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.10.0"

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

[[deps.Glob]]
git-tree-sha1 = "4df9f7e06108728ebf00a0a11edee4b29a482bb2"
uuid = "c27321d9-0574-5035-807b-f59d2c89b15c"
version = "1.3.0"

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

[[deps.Interpolations]]
deps = ["AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "b7bc05649af456efc75d178846f47006c2c4c3c7"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.13.6"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "b3364212fb5d870f724876ffcd34dd8ec6d98918"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.7"

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
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

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

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "09e4b894ce6a976c354a69041a04748180d43637"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.15"

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
git-tree-sha1 = "1ea784113a6aa054c5ebd95945fa5e52c2f378e7"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.7"

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
deps = ["Dates"]
git-tree-sha1 = "0044b23da09b5608b4ecacb4e5e6c6332f833a7e"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.3.2"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlotlyBase]]
deps = ["ColorSchemes", "Dates", "DelimitedFiles", "DocStringExtensions", "JSON", "LaTeXStrings", "Logging", "Parameters", "Pkg", "REPL", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "180d744848ba316a3d0fdf4dbd34b77c7242963a"
uuid = "a03496cd-edff-5a9b-9e67-9cda94a718b5"
version = "0.8.18"

[[deps.PlutoDevMacros]]
deps = ["MacroTools", "Requires"]
git-tree-sha1 = "994167def8f46d3be21783a76705228430e29632"
uuid = "a0499f29-c39b-4c5c-807c-88074221b949"
version = "0.4.5"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "8d1f54886b9037091edf146b517989fc4a09efec"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.39"

[[deps.PlutoUtils]]
deps = ["AbstractPlutoDingetjes", "Chain", "Colors", "DocStringExtensions", "Glob", "HypertextLiteral", "OrderedCollections", "PlutoDevMacros", "PlutoUI", "PrettyTables", "Reexport", "Requires", "StaticArrays", "UUIDs"]
git-tree-sha1 = "3f8dfe27dbb980ad5e83ecd641ded8eed91f3265"
uuid = "ed5d0301-4775-4676-b788-cf71e66ff8ed"
version = "0.5.9"

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

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["Serialization"]
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

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

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
git-tree-sha1 = "a9e798cae4867e3a41cae2dd9eb60c047f1212db"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.6"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "9f8a5dc5944dc7fbbe6eb4180660935653b0a9d9"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.0"

[[deps.StaticArraysCore]]
git-tree-sha1 = "66fe9eb253f910fe8cf161953880cfdaef01cdf0"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.0.1"

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
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "5ce79ce186cc678bbb5c5681ca3379d1ddae11a1"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.7.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

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
# ╠═1c0cf7ad-9579-4593-ae1f-441788c049ae
# ╠═47a21d25-7710-4a93-8cbd-ad93cee68b74
# ╠═d2056fd0-5419-11ec-2226-0b24c663de65
# ╠═2c06b5c1-0345-444f-9664-338a60fedc83
# ╠═6c71d117-c2e5-4eb4-8abb-39109ffa30db
# ╠═78ad9ad5-7c06-41e1-940d-1aad1de7d5f2
# ╠═35450d8d-cda5-46ba-884e-636abf023ac1
# ╠═a1158d4b-38cc-4c9b-a949-e94e98c2b5f6
# ╟─0345b4e1-49e6-4fc2-bc0e-6f179a0b1792
# ╠═94032874-a85b-4632-bc24-c8c0154b518a
# ╠═4d74427d-21e3-4736-ab29-ef65c7ea0490
# ╠═a1c5f1c0-adc2-43f0-bd63-e4c6c7eca52b
# ╠═5907b0b7-f81b-4cc3-a585-8f52fd191bae
# ╠═37772abb-c600-4beb-a92e-7da3c3021370
# ╠═27cda06f-de16-47f6-be79-324963b155c9
# ╠═f03f7ebb-dad0-4766-8615-afde6478f59d
# ╠═1f389615-64d5-4084-b4a7-e7838a19ff9b
# ╠═eac20e32-59a3-44c3-9c68-31be56ad8ee9
# ╠═ff6dd97b-66db-48e8-bd6d-bfccff28ba39
# ╠═e48e3061-4004-4cee-8a93-7218e3ba5bbc
# ╠═d2fbd3c2-215e-4859-861e-03a0b2efdc9b
# ╠═b85c0bff-162d-402e-8d85-44d46b4af8a8
# ╠═c2062a21-8027-4d49-85d6-6666f761cabb
# ╠═29d90f34-0aeb-48d3-80b4-759bb0b77c2f
# ╠═712f157b-8cc9-4034-8f54-8b57a7953683
# ╠═4b14629a-b25f-4f3c-94bd-14e77c044695
# ╠═215ca230-77fd-4bbb-a4f7-8af1a0c9bcb3
# ╠═57500cea-75e5-473a-83c6-635f3ec3fa95
# ╠═774fc0f5-a11f-42eb-9bc0-3faa13ed24f4
# ╠═5203bbc5-4746-48cf-b209-ec618b223af4
# ╠═154d61e7-9014-4cc0-b817-3257fcfe1826
# ╠═602c9baf-96e5-4e2d-9e84-3bf177d137e9
# ╠═039bfd43-c1e1-4086-b108-2079ce5e0453
# ╠═1696f2e6-1753-41b4-b6c1-ab8038680c91
# ╠═03b7a36a-8a82-4834-baf1-a4114f6afca3
# ╠═b959608b-2a49-4397-b974-9094480596ab
# ╠═d0259918-920a-4b8e-b33d-a93da47425de
# ╠═b6909885-61c0-4403-9639-c00ca7cd55e1
# ╠═353923b3-ab53-4fdf-8eb3-f7a4d126bdf3
# ╠═5953946c-52fc-4297-bc3a-9aa2759c8be4
# ╠═080fdaa3-b9d7-4030-8bbc-88abcc5d37a3
# ╠═2f49954a-8e37-4d84-accc-f5537e734331
# ╠═1a6110f9-4086-43ee-8a12-26755ac0ee8a
# ╠═87ca96a1-8f5a-4534-9e55-8f4bc4f0e437
# ╠═6c13719f-3b65-4550-88e9-6f1c5cd84394
# ╠═c68bc281-eb52-43db-b216-fa01350b932a
# ╠═0dcb4a08-b8ed-4a4b-8be2-7ccda056846a
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
