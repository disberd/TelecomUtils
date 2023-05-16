import Unitful: °, km
using Rotations
import LinearAlgebra: normalize
import TelecomUtils: earth_intersection, ExtraOutput
using TelecomUtils
using StaticArrays
using Test

sp_ell = SphericalEllipsoid()
### ECEF <-> LLA ###
lla2ecef = ECEFfromLLA(sp_ell)
@test lla2ecef(LLA(0, 0, 0)) ≈ [6371e3, 0, 0]
lla1 = LLA(30°, 45°, 100km)
@test lla2ecef(lla1) |> inv(lla2ecef) ≈ lla1 # Test transformation inversion

### Earth Intersection ###
@testset "Earth Intersection" begin
    sat_lla = LLA(0, 0, 500km)
    sat_ecef = lla2ecef(sat_lla)

    n̂ = normalize(-sat_ecef)
    intersection_ecef = earth_intersection(n̂, sat_ecef, sp_ell.a, sp_ell.b)
    @test intersection_ecef ≈ [sp_ell.a, 0, 0]
    # Test that if the direction is opposite of the earth the result is NaN
    @test all(map(isnan, earth_intersection(-n̂, sat_ecef, sp_ell.a, sp_ell.b)))
end

### Compute sat positions ###
@testset "Compute Sat Positions" begin
    h = 6371e3 + 735e3
    lla = LLA(10°, 25°, 0km)
    era2ecef = ECEFfromERA(lla)
    el, az = 35°, 80°
    satecef = compute_sat_position(era2ecef, el, az; h)
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
    @test target_uv[1] == [-1, 0]
    @test target_uv[2] == [0, 1]
    @test target_uv[3] == [1, 0]
    @test target_uv[4] == [0, -1]


    # We now test that targets behind the reference direction are not visible (NaN)
    sat_lla = LLA(0, 0, 600km)
    target_lla = LLA(0, 0, 610km)
    lla2uv = UVfromLLA(sat_lla)
    target_uv = lla2uv(target_lla)
    @test all(isnan.(target_uv))
    # We test that if we rotate by 180 degrees around Y the matrix we actually find a valid uv
    lla2uv = UVfromLLA(lla2uv.origin, RotY(180°) * lla2uv.R, lla2uv.ellipsoid)
    @test all(map(!isnan, lla2uv(target_lla)))
end

@testset "LLAfromUV" begin
    sp = SphericalEllipsoid()
    l2e = ECEFfromLLA(sp)
    # We find the pointing that corresponds to Edge of Earth and check various combination in its vicinity
    sat_lla = LLA(0, 0, 600km)
    sat_ecef = l2e(sat_lla)
    eoe_scan = asin(sp.a / (sp.a + sat_lla.alt))
    u = sin(eoe_scan)
    uv2lla = LLAfromUV(sat_lla; ellipsoid=sp)
    @test !isnan(uv2lla((u * (1 - eps()), 0))) # We should find a solution because we are pointing slightly less than EoE
    @test isnan(uv2lla((u * (1 + eps()), 0))) # We should not find a solution because we are pointing slightly more than EoE
    @test !isnan(uv2lla((u * (1 - eps()), 0), h=100e3)) # We should find a solution because we are looking at 100km above earth
    @test !isnan(uv2lla((u * (1 + eps()), 0), h=100e3)) # We should find a solution because we are looking at 100km above earth
    @test isnan(uv2lla((u * (1 - eps()), 0), h=700e3)) # We should not find a solution because we are looking at 100km above the satellite alitude and with an angle slightly lower than eoe scan, so the corresponding valid point in the pointing direction is located behind earth
    @test !isnan(uv2lla((u * (1 + eps()), 0), h=700e3)) # We should find a solution because we are pointing more than eoe_scan so the earth is not blocking the view of the corresponding point


    lla2uv = inv(uv2lla)
    target_uv = SA_F64[0.1, 0.1]
    target_lla, r = uv2lla(target_uv, ExtraOutput())
    uv2, r2 = lla2uv(target_lla, ExtraOutput())
    @test uv2 ≈ target_uv
    @test r2 ≈ r
end
