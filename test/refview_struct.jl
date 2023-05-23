import Unitful: °, km
import TelecomUtils: wgs84_ellipsoid, ExtraOutput
using TelecomUtils
using Test
using StaticArrays

Base.isnan(v::StaticArray) = any(isnan, v)

sp_ell = SphericalEllipsoid()
lla2ecef = ECEFfromLLA(sp_ell)

em = EarthModel(sp_ell)
sat_lla = LLA(0, 0, 700km)

sv = SatView(sat_lla, em)

@testset "SatView Creation" begin
    @test SatView(sat_lla, em; face=:PositiveY).face === TelecomUtils.PositiveY
    @test SatView(sat_lla, em; face=-1).face === TelecomUtils.NegativeX
    @test SatView(sat_lla, em; face=TelecomUtils.PositiveX).face === TelecomUtils.PositiveX
end

@testset "Get Range" begin
    # 2D Pointing based
    @test get_range(sv, (0, 0)) ≈ sat_lla.alt
    @test get_range(sv, (0, 0); h=1e3) ≈ sat_lla.alt - 1e3
    @test get_range(sv, (0, 0); face=:NegativeZ) |> isnan # Default altitude is 0m so there is no intersection from face -Z
    @test get_range(sv, (0, 0); h=sv.lla.alt + 100e3, face=:NegativeZ) ≈ 100e3 # When providing an altitude above the satellite, the intersection from -Z is found

    # ECEF/LLA/ReferenceView based
    @test get_range(sv, SatView(LLA(0, 0, 600km), em)) ≈ 100e3 # 2 ReferenceViews
    @test get_range(sv, SatView(LLA(0, 0, 800km), em)) |> isnan # The target is above the satellite, so not visible from the reference face
    @test get_range(sv, SatView(LLA(0, 0, 800km), em); face=:NegativeZ) ≈ 100e3 # The target is above the satellite, so it's visible from -Z

    # Throws with different earth models
    @test_throws "EarthModel" get_range(sv, SatView(LLA(0, 0, 500km), EarthModel())) # Different EarthModels
end

@testset "Get Pointing/LLA/ECEF" begin
    lla_ref = LLA(1°, 1°, 1km)
    ecef_ref = lla2ecef(lla_ref)
    ref_uv = get_pointing(sv, lla_ref)
    @test ref_uv ≈ get_pointing(sv, ecef_ref)
    @test get_lla(sv, ref_uv; h=lla_ref.alt) ≈ lla_ref
    @test get_ecef(sv, ref_uv; h=lla_ref.alt) ≈ ecef_ref

    above = LLA(sv.lla.lat, sv.lla.lon, sv.lla.alt + 100e3)
    below = LLA(sv.lla.lat, sv.lla.lon, sv.lla.alt - 100e3)
    @test isapprox(get_pointing(sv, below), SA_F64[0, 0]; atol=1e-10)
    @test get_pointing(sv, above) |> isnan
    @test isapprox(get_pointing(sv, above; face=:NegativeZ), SA_F64[0, 0]; atol=1e-10)
    # We test the X and Y faces.
    @test isapprox(get_pointing(sv, below; face=:PositiveX), SA_F64[-1, 0]; atol=1e-10)
    @test isapprox(get_pointing(sv, below; face=:NegativeX), SA_F64[1, 0]; atol=1e-10)
    @test isapprox(get_pointing(sv, below; face=:PositiveY), SA_F64[0, -1]; atol=1e-10)
    @test isapprox(get_pointing(sv, below; face=:NegativeY), SA_F64[0, 1]; atol=1e-10)

    # Test a point with manually computed theta angle
    R_e = em.ellipsoid.a
    sin_ρ = R_e / (R_e + sv.lla.alt)
    θ = asin(sin_ρ)
    lat = acos(sin_ρ)
    thetaphi = get_pointing(sv, LLA(lat - 1e-5, 0, 0); pointing_type=:thetaphi) # We need 1e-5 to allow for tolerance with earth intersection
    @test isapprox(thetaphi, SA_F64[θ, π/2]; atol=1e-5)

    # Test extra output
    _, xyz, block_data = get_pointing(sv, LLA(0, 0, 0), ExtraOutput())
    @test isapprox(xyz, SA_F64[0, 0, sv.lla.alt]; atol=1e-10) && !block_data.blocked
    uv, xyz, block_data = get_pointing(sv, LLA(0, 180°, 0), ExtraOutput())
    @test isnan(uv) && isnan(xyz) && block_data.blocked # The target is blocked by earth
    uv, xyz, block_data = get_pointing(sv, LLA(0, 0, 0), ExtraOutput(); face=:NegativeZ)
    @test isnan(uv) && isapprox(xyz, SA_F64[0, 0, -sv.lla.alt]; atol=1e-10) && !block_data.blocked # The point is behind the satellite reference face
    _, xyz, block_data = get_pointing(sv, LLA(0, 0, 0), ExtraOutput(); face=:PositiveX)
    @test isapprox(xyz, SA_F64[-sv.lla.alt, 0, 0]; atol=1e-10)

    # Test error with different earthmodel
    @test_throws "EarthModel" get_pointing(sv, SatView(LLA(0, 0, 500km), EarthModel())) # Different EarthModels

    # Get ECEF
    @test get_ecef(sv, (0, 0)) ≈ lla2ecef(LLA(0, 0, 0))
    @test get_ecef(sv, (0, 0); face=:NegativeZ) |> isnan
    @test get_ecef(sv, (-1, 0); face=:PositiveX) ≈ lla2ecef(LLA(0, 0, 0))
    @test get_ecef(sv, (1, 0); face=:NegativeX) ≈ lla2ecef(LLA(0, 0, 0))
    @test get_ecef(sv, (0, -1); face=:PositiveY) ≈ lla2ecef(LLA(0, 0, 0))
    @test get_ecef(sv, (0, 1); face=:NegativeY) ≈ lla2ecef(LLA(0, 0, 0))
end

@testset "Get Mutual Pointing" begin
    sv1 = SatView(LLA(0, 0, 800km), em)
    sv2 = SatView(LLA(0, 0, 1000km), em)
    p1, p2 = get_mutual_pointing(sv1, sv2)
    @test isnan(p1) && isapprox(p2, SA_F64[0, 0]; atol=1e-10) # The reference face for sv1 is still nadir, so it doesn't see sv2

    p1, p2 = get_mutual_pointing(sv1, sv2; faces=(:NegativeZ, :PositiveZ))
    @test isapprox(p1, SA_F64[0, 0]; atol=1e-10) && isapprox(p2, SA_F64[0, 0]; atol=1e-10)
end

@testset "Get ERA" begin
    # We test that a non-visible point is NaN
    @test get_era(UserView(LLA(40°, -39°, 0), em), sv) |> isnan
    # We test that a visible point is not NaN, and with the expected value
    @test get_era(UserView(LLA(0, 0, 500km), em), sv) ≈ ERA(90°, 200km, 0°)
    @test_throws "EarthModel" get_era(UserView(LLA(0, 0, 0), EarthModel()), sv)
    @test_throws "UserView" get_era(sv, LLA(0, 0, 0))
end

@testset "Get Distance on Earth" begin
    r = 6371e3
    sp = SphericalEllipsoid(r)
    em = EarthModel(sp)
    sv = SatView(LLA(0, 0, 700km), em)
    target_dist = 2π * r / 360
    @test get_distance_on_earth(LLA(0°, 0°, 0), LLA(1°, 0°, 0); em) ≈ target_dist
end

@testset "Get Nadir Beam Diameter" begin
    # Check that the beam diameter does not depend on latitude on a spherical ellipsoid
    @test get_nadir_beam_diameter(SatView(LLA(90°, 0°, 735km), EarthModel()), 55) ≈ get_nadir_beam_diameter(SatView(LLA(0°, 0°, 735km), EarthModel()), 55)
    # Test that the wgs84 ellipsoid makes a difference in the beam diameter computation
    @test get_nadir_beam_diameter(SatView(LLA(90°, 0°, 735km), EarthModel(wgs84_ellipsoid)), 55) ≉ get_nadir_beam_diameter(SatView(LLA(0°, 0°, 735km), EarthModel(wgs84_ellipsoid)), 55)
end

@testset "Get Visibility" begin
    above = SatView(LLA(sv.lla.lat, sv.lla.lon, sv.lla.alt + 100e3), em)
    below = SatView(LLA(sv.lla.lat, sv.lla.lon, sv.lla.alt - 100e3), em)
    user_30_deg = get_ecef(sv, (30°, 0°); pointing_type = :thetaphi) |> x -> UserView(x, em)
    sat_60_deg = get_ecef(sv, (60°, 0°); h = sv.lla.alt - 200e3, pointing_type = :thetaphi) |> x -> SatView(x, em; face = -3)
    @test get_visibility(sv, below)
    @test !get_visibility(sv, above)
    @test get_visibility(sv, above; boresight = :NegativeZ)
    @test get_visibility(sv, user_30_deg; fov = 30°)
    @test !get_visibility(sv, user_30_deg; fov = 29.99°)
    @test get_visibility(sv, sat_60_deg; fov = (60 + 1e-5)°)
    @test !get_visibility(sv, sat_60_deg; fov = (60 - 1e-5)°)

    # Test a point with manually computed theta angle
    R_e = em.ellipsoid.a
    sin_ρ = R_e / (R_e + sv.lla.alt)
    θ = asin(sin_ρ)
    lat = acos(sin_ρ)
    visible, theta, _ = get_visibility(sv, LLA(lat - 1e-5, 0, 0), ExtraOutput()) # We need 1e-5 to allow for tolerance with earth intersection
    @test isapprox(theta, θ; atol=1e-5)

    # Test extra output
    visible, θ, _ = get_visibility(sv, LLA(0, 0, 0), ExtraOutput())
    @test abs(θ) < 1e-5
    visible, θ, _ = get_visibility(sv, LLA(0, 180°, 0), ExtraOutput())
    @test isnan(θ) # The target is blocked by earth
    visible, θ, _ = get_visibility(sv, LLA(0, 0, 0), ExtraOutput(); boresight=:NegativeZ)
    @test θ ≈ π && !visible # The point is behind the satellite reference face
    visible, θ, _ = get_visibility(sv, LLA(0, 0, 0), ExtraOutput(); boresight=:PositiveX)
    @test θ ≈ π/2 && visible

    # Test error with different earthmodel
    @test_throws "EarthModel" get_visibility(sv, SatView(LLA(0, 0, 500km), EarthModel())) # Different EarthModels
end

@testset "Get Mutual Visibitiliy" begin
    sv1 = SatView(LLA(0, 0, 800km), em)
    sv2 = SatView(LLA(0, 0, 1000km), em)
    visible, (fwd, rtn) = get_mutual_visibility(sv1, sv2, ExtraOutput())
    @test !visible && fwd === rtn # The reference face for sv1 is still nadir, so it doesn't see sv2. If we don't set short_circuit to false, fwd and rtn are equivalent as rtn is not computed
    visible, ((fwd_visible, _), (rtn_visible, _)) = get_mutual_visibility(sv1, sv2, ExtraOutput(); short_circuit = false)
    @test !visible && !fwd_visible && rtn_visible # The reference face for sv1 is still nadir, so it doesn't see sv2. but sv2 sees sv1

    visible = get_mutual_visibility(sv1, sv2; boresights=(:NegativeZ, :PositiveZ))
    @test visible
end