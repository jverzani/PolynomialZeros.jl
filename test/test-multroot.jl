using PolynomialZeros
const AGCD = PolynomialZeros.AGCD
const MultRoot = PolynomialZeros.MultRoot
using Polynomials
using Test



@testset "agcd-Float64" begin

    x = variable()

    p = prod(x-i for i in 1:20)
    u,v,w,err = AGCD.agcd(p)
    @test Polynomials.degree(v) == Polynomials.degree(p)

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
    u,v,w,err = AGCD.agcd(p)

    @test degree(v) == 3


    T = BigFloat
    x = variable(T)
    p = (x-1)^2 * (x-2)^2 * (x+3)^3
    u,v,w,err = AGCD.agcd(p^15)

    @test degree(v) == 3 # fails for Float64

end
