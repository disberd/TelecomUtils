import Unitful: °, km
import LinearAlgebra: normalize
import TelecomUtils: earth_intersection

@testset "Satview Transformations" begin
    sp_ell = SphericalEllipsoid()
    ### ECEF <-> LLA ###
    lla2ecef = ECEFfromLLA(sp_ell)
    @test lla2ecef(LLA(0,0,0)) ≈ [6371e3, 0, 0]
    lla1 = LLA(30°, 45°, 100km)
    @test lla2ecef(lla1) |> inv(lla2ecef) ≈ lla1 # Test transformation inversion

    ### Earth Intersection ###
    @testset "Earth Intersection" begin
        sat_lla = LLA(0,0,500km)
        sat_ecef = lla2ecef(sat_lla)

        n̂ = normalize(-sat_ecef)
        intersection_ecef = earth_intersection(n̂, sat_ecef, sp_ell.a, sp_ell.b)
        @test intersection_ecef ≈ [sp_ell.a, 0, 0]
    end

    ### Compute sat positions ###
    @testset "Compute Sat Positions" begin
        h = 6371e3 + 735e3
        lla = LLA(10°, 25°, 0km)
        era2ecef = ECEFfromERA(lla)
        el, az = 35°, 80°
        satecef = compute_sat_position(era2ecef, el, az;h)
        invera = inv(era2ecef)(satecef)
        @test invera.el ≈ el && invera.az ≈ az
    end
end
