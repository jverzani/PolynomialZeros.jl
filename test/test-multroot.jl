using PolynomialZeros
const AGCD = PolynomialZeros.AGCD
const MultRoot = PolynomialZeros.MultRoot
using Polynomials
using Test

multroot=PolynomialZeros.MultRoot.multroot
pejroot=PolynomialZeros.MultRoot.pejroot
agcd = PolynomialZeros.AGCD.agcd
identify_z0s_ls = PolynomialZeros.MultRoot.identify_z0s_ls

x = variable()
_poly(zs,ls) = prod((x-z)^l for (z,l) in zip(zs, ls))

@testset "agcd-Float64" begin

    x = variable()

    p = (x-1)^3
    u,v,w,err = AGCD.agcd(p)
    @test Polynomials.degree(v) == 1

    p = prod(x-i for i in 1:6)
    u,v,w,err = AGCD.agcd(p)
    @test Polynomials.degree(v) == Polynomials.degree(p)

    n = 4
    p = prod((x-i)^i for i in 1:n)
    u,v,w,err = AGCD.agcd(p)
    @test Polynomials.degree(v) == n

    # can fails
    n = 6
    p = prod((x-i)^i for i in 1:n)
    u,v,w,err = AGCD.agcd(p)
    @test Polynomials.degree(v) >= n

    # use big
    n = 6
    p = prod((x-i)^i for i in 1:n)
    u,v,w,err = AGCD.agcd(convert(Poly{BigFloat}, p))
    @test Polynomials.degree(v) == n

    T = Float64
    x = variable(Complex{T})
    p = (x-im)^2 * (x-1)^2 * (x+2im)^2
    u,v,w,err = AGCD.agcd(p)
    @test length(v)-1 == 3


end

@testset "agcd-othertypes" begin


    T = Float32
    x = variable(T)
    p = (x-1)^2 * (x-2)^2 * (x+3)^3
    u,v,w,err = AGCD.agcd(p, θ=1e-4)

    @test degree(v) == 3


    T = BigFloat
    x = variable(T)
    p = (x-1)^2 * (x-2)^2 * (x+3)^3
    u,v,w,err = AGCD.agcd(p^15)

    @test degree(v) == 3 # fails for Float64

end

@testset "pejroot" begin


    zs, ls = [1.0,2,3,4], [4,3,2,1]
    p = _poly(zs, ls)

    delta = 0.1 # works
    z0 = zs + delta*[1,-1,1,-1]
    z1 = pejroot(p, z0, ls)
    @test !all(sort(z0) .== sort(z1))

    delta = 0.2 # fails
    z0 = zs + delta*[1,-1,1,-1]
    z1 = pejroot(p, z0, ls)
    @test all(sort(z0) .== sort(z1))

    ls = [30,20,10,5]
    p = _poly(zs, ls)

    delta = 0.01 # works
    z0 = zs + delta*[1,-1,1,-1]
    z1 = pejroot(p, z0, ls)
    @test !all(sort(z0) .== sort(z1))

end

@testset "identify_ls" begin


    zs, ls = [1.0,2,3,4], [4,3,2,1]
    p = _poly(zs, ls)
    _zs, _ls = identify_z0s_ls(coeffs(p))
    @test all(sort(ls) .== sort(_ls))

    p = _poly(zs, 3*ls)
    _zs, _ls = identify_z0s_ls(coeffs(p))
    @test length(ls) != length(_ls) # fails w/o preconditioning

    n = 4
    zs, ls = cumsum(ones(n)), cumsum(ones(Int, n))
    p = _poly(zs, ls)
    _zs, _ls = identify_z0s_ls(coeffs(p))
    @test all(sort(ls) .== sort(_ls))

    n = 5
    zs, ls = cumsum(ones(n)), cumsum(ones(Int, n))
    p = _poly(zs, ls)
    _zs, _ls = identify_z0s_ls(coeffs(p))
    @test !(length(ls) == length(_ls))

end



@testset "multroot" begin

    zs, ls = [1.0,2,3,4], [4,3,2,1]
    _zs, _ls = multroot(_poly(zs, ls))
    @test all(sort(ls) .== sort(_ls))

    _zs, _ls = multroot(_poly(zs, 3ls))
    @test !all(sort(2ls) .== sort(_ls)) # XXX fails!


    _zs, _ls = multroot(_poly(big.(zs), 3ls))
    @test all(sort(3ls) .== sort(_ls)) # passes


    delta = 0.01
    zs, ls = [1-delta, 1, 1+delta], [5,4,3]
    _zs, _ls = multroot(_poly(zs, ls))
    @test all(sort(ls) .== sort(_ls))

    n = 20
    zs,ls = collect(1.0:n), ones(Int, n)
    _zs, _ls = multroot(_poly(zs, ls))
    @test all(sort(ls) .== sort(_ls))

    _zs, _ls = multroot(_poly(zs, 2ls))
    @test !(length(ls) == length(_ls)) ##XXX fails!

    n = 10
    zs,ls = collect(1.0:n), ones(Int, n)
    _zs, _ls = multroot(_poly(zs, 2ls)) #fails *badly*
    @test !(length(_ls) == length(ls))

    n = 10
    zs,ls = collect(1.0:n), ones(Int, n)
    _zs, _ls = multroot(_poly(big.(zs), 2ls)) # works now
    @test all(sort(ls) .== sort(_ls))

    n = 5
    zs,ls = collect(1.0:n), ones(Int, n)
    _zs, _ls = multroot(_poly(zs, 4ls))
    @test all(sort(4ls) .== sort(_ls))

end
