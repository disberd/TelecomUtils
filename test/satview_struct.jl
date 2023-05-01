import Unitful: °, km
import TelecomUtils: wgs84_ellipsoid, ExtraOutput

@testset "Satview Struct" begin
    sp_ell = SphericalEllipsoid()
    lla2ecef = ECEFfromLLA(sp_ell)

    em = EarthModel(sp_ell)
    sat_lla = LLA(0,0,700km)

    sv = SatView(sat_lla, em)

    @testset "SatView Creation" begin
        @test SatView(sat_lla, em; face = :PositiveY).face === TelecomUtils.PositiveY
        @test SatView(sat_lla, em; face = -1).face === TelecomUtils.NegativeX
        @test SatView(sat_lla, em; face = TelecomUtils.PositiveX).face === TelecomUtils.PositiveX
    end

    @testset "Get Range" begin
        # 2D Pointing based
        @test get_range(sv, (0,0)) ≈ sat_lla.alt
        @test get_range(sv, (0,0); h= 1e3) ≈ sat_lla.alt - 1e3
        @test get_range(sv, (0,0); face=:NegativeZ) |> isnan # Default altitude is 0m so there is no intersection from face -Z
        @test get_range(sv, (0,0); h = sv.lla.alt + 100e3, face=:NegativeZ) ≈ 100e3 # When providing an altitude above the satellite, the intersection from -Z is found

        # ECEF/LLA/ReferenceView based
        @test get_range(sv, SatView(LLA(0,0,600km), em)) ≈ 100e3 # 2 ReferenceViews
        @test get_range(sv, SatView(LLA(0,0,800km), em)) |> isnan # The target is above the satellite, so not visible from the reference face
        @test get_range(sv, SatView(LLA(0,0,800km), em); face = :NegativeZ) ≈ 100e3 # The target is above the satellite, so it's visible from -Z

        # Throws with different earth models
        @test_throws "EarthModel" get_range(sv, SatView(LLA(0,0,500km), EarthModel())) # Different EarthModels
    end

    @testset "Get Pointing/LLA/ECEF" begin
        lla_ref = LLA(1°, 1°, 1km)
        ecef_ref = lla2ecef(lla_ref)
        ref_uv = get_pointing(sv, lla_ref)
        @test ref_uv ≈ get_pointing(sv, ecef_ref)
        @test get_lla(sv, ref_uv;h = lla_ref.alt) ≈ lla_ref
        @test get_ecef(sv, ref_uv;h = lla_ref.alt) ≈ ecef_ref     

        above = LLA(sv.lla.lat,sv.lla.lon,sv.lla.alt + 100e3)
        below = LLA(sv.lla.lat,sv.lla.lon,sv.lla.alt - 100e3)
        @test isapprox(get_pointing(sv, below),SA_F64[0,0]; atol = 1e-10)
        @test get_pointing(sv, above) |> isnan
        @test isapprox(get_pointing(sv, above; face = :NegativeZ),SA_F64[0,0]; atol = 1e-10)
        # We test the X and Y faces.
        @test isapprox(get_pointing(sv, below; face = :PositiveX),SA_F64[-1,0]; atol = 1e-10)
        @test isapprox(get_pointing(sv, below; face = :NegativeX),SA_F64[1,0]; atol = 1e-10)
        @test isapprox(get_pointing(sv, below; face = :PositiveY),SA_F64[0,-1]; atol = 1e-10)
        @test isapprox(get_pointing(sv, below; face = :NegativeY),SA_F64[0,1]; atol = 1e-10)

        # Test a point with manually computed theta angle
        R_e = em.ellipsoid.a
        sin_ρ = R_e/(R_e + sv.lla.alt)
        θ = asin(sin_ρ)
        lat = acos(sin_ρ)
        thetaphi = get_pointing(sv, LLA(lat - 1e-5, 0, 0); pointing_type = :thetaphi) # We need 1e-5 to allow for tolerance with earth intersection
        @test isapprox(thetaphi, SA_F64[θ, π/2]; atol = 1e-5)

        # Test extra output
        _, xyz = get_pointing(sv, LLA(0,0,0), ExtraOutput())
        @test isapprox(xyz, SA_F64[0,0,sv.lla.alt]; atol = 1e-10)
        _, xyz = get_pointing(sv, LLA(0,0,0), ExtraOutput(); face = :NegativeZ)
        @test isnan(xyz) # The point is behind the satellite reference face
        _, xyz = get_pointing(sv, LLA(0,0,0), ExtraOutput(); face = :PositiveX)
        @test isapprox(xyz, SA_F64[-sv.lla.alt, 0,0]; atol = 1e-10)

        # Test error with different earthmodel
        @test_throws "EarthModel" get_pointing(sv, SatView(LLA(0,0,500km), EarthModel())) # Different EarthModels

        # Get ECEF
        @test get_ecef(sv, (0,0)) ≈ lla2ecef(LLA(0,0,0))
        @test get_ecef(sv, (0,0); face = :NegativeZ) |> isnan
        @test get_ecef(sv, (-1,0); face = :PositiveX) ≈ lla2ecef(LLA(0,0,0))
        @test get_ecef(sv, (1,0); face = :NegativeX) ≈ lla2ecef(LLA(0,0,0))
        @test get_ecef(sv, (0,-1); face = :PositiveY) ≈ lla2ecef(LLA(0,0,0))
        @test get_ecef(sv, (0,1); face = :NegativeY) ≈ lla2ecef(LLA(0,0,0))
    end

    @testset "Get ERA" begin
        # We test that a non-visible point is NaN
        @test get_era(UserView(LLA(40°, -39°, 0), em), sv) |> isnan
        # We test that a visible point is not NaN, and with the expected value
        @test get_era(UserView(LLA(0,0,500km), em), sv) ≈ ERA(90°, 200km, 0°)
        @test_throws "EarthModel" get_era(UserView(LLA(0,0,0), EarthModel()), sv)
        @test_throws "UserView" get_era(sv, LLA(0,0,0))
    end

    @testset "Get Distance on Earth" begin
        r = 6371e3
        sp = SphericalEllipsoid(r)
        em = EarthModel(sp)
        sv = SatView(LLA(0,0,700km), em)
        target_dist = 2π*r / 360
        @test get_distance_on_earth(LLA(0°, 0°, 0), LLA(1°, 0°, 0); em) ≈ target_dist
    end

    @testset "Get Nadir Beam Diameter" begin
        # Check that the beam diameter does not depend on latitude on a spherical ellipsoid
        @test get_nadir_beam_diameter(SatView(LLA(90°,0°,735km), EarthModel()), 55) ≈ get_nadir_beam_diameter(SatView(LLA(0°,0°,735km), EarthModel()), 55)
        # Test that the wgs84 ellipsoid makes a difference in the beam diameter computation
        @test get_nadir_beam_diameter(SatView(LLA(90°,0°,735km), EarthModel(wgs84_ellipsoid)), 55) ≉ get_nadir_beam_diameter(SatView(LLA(0°,0°,735km), EarthModel(wgs84_ellipsoid)), 55)
    end
end