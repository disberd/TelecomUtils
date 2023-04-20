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

    @testset "UVfromLLA" begin
        sat_lla = LLA(0°, 0°, 600km)
        target_lla = [
            LLA(0°, 1°, 0km), # Right - U Negative, V 0
            LLA(1°, 0°, 0km), # Top - U 0, V Positive
            LLA(0°, -1°, 0km), # Left - U Positive, V 0
            LLA(-1°, 0°, 0km), # Bottom - U 0, V Negative
        ]
        target_uv = map(target_lla) do lla
            UVfromLLA(sat_lla; ellipsoid=SphericalEllipsoid())(lla) |> normalize
        end
        @test target_uv[1] == [-1,0]
        @test target_uv[2] == [0,1]
        @test target_uv[3] == [1,0]
        @test target_uv[4] == [0,-1]


        # We now test that targets behind the reference direction are not visible (NaN)
        target_uv = UVfromLLA(LLA(0,0,600km))(LLA(0,0,610km))
        @test all(isnan.(target_uv))
    end
end
