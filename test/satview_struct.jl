import Unitful: °, km
import TelecomUtils: wgs84_ellipsoid

@testset "Satview Struct" begin
    sp_ell = SphericalEllipsoid()
    lla2ecef = ECEFfromLLA(sp_ell)

    em = EarthModel(sp_ell)
    sat_lla = LLA(0,0,700km)

    sv = SatView(sat_lla, em)

    @testset "Get Range" begin
        @test get_range(sv, (0,0)) ≈ sat_lla.alt
        @test get_range(sv, (0,0), 1e3) ≈ sat_lla.alt - 1e3
    end

    @testset "Get Pointing/LLA/ECEF" begin
        lla_ref = LLA(1°, 1°, 1km)
        ecef_ref = lla2ecef(lla_ref)
        ref_uv = get_pointing(sv, lla_ref)
        @test ref_uv ≈ get_pointing(sv, ecef_ref)
        @test get_lla(sv, ref_uv;h = lla_ref.alt) ≈ lla_ref
        @test get_ecef(sv, ref_uv;h = lla_ref.alt) ≈ ecef_ref     
    end

    @testset "Get ERA" begin
        # We test that a non-visible point is NaN
        @test get_era(sv, LLA(40°, -39°, 0)) |> isnan
        # We test that a visible point is not NaN, and with the expected value
        @test get_era(sv,LLA(0,0,500km)) ≈ ERA(90°, 200km, 0°)
    end

    @testset "Get Distance on Earth" begin
        r = 6371e3
        sp = SphericalEllipsoid(r)
        em = EarthModel(sp)
        sv = SatView(LLA(0,0,700km), em)
        target_dist = 2π*r / 360
        @test get_distance_on_earth(sv, LLA(0°, 0°, 0), LLA(1°, 0°, 0)) ≈ target_dist
    end

    @testset "Get Nadir Beam Diameter" begin
        # Check that the beam diameter does not depend on latitude on a spherical ellipsoid
        @test get_nadir_beam_diameter(SatView(LLA(90°,0°,735km), EarthModel()), 55) ≈ get_nadir_beam_diameter(SatView(LLA(0°,0°,735km), EarthModel()), 55)
        # Test that the wgs84 ellipsoid makes a difference in the beam diameter computation
        @test get_nadir_beam_diameter(SatView(LLA(90°,0°,735km), EarthModel(wgs84_ellipsoid)), 55) ≉ get_nadir_beam_diameter(SatView(LLA(0°,0°,735km), EarthModel(wgs84_ellipsoid)), 55)
    end
end