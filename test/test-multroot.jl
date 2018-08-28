using PolynomialZeros
const AGCD = PolynomialZeros.AGCD
const MultRoot = PolynomialZeros.MultRoot
using Polynomials
using Test


@testset "agcd routines" begin

    x = variable()
    p = (x-1)^2 * (x-2)^2 * (x-3)^2
    q = polyder(p)
    k = 2 # non-singular  (degree(p) - degree(gcd(p,p')) = k -- smallest k with A singular
    psi, xs = AGCD.lemma24(p.a, q.a, k)
    @test psi > 1e-8

    k = 3 # singular 6 - 3 = k = 3
    psi, xs = AGCD.lemma24(p.a, q.a, k)
    @test psi <= 1e-8


    

    p =  (x-1)^3 * (x-2)^3 * (x-3)^3 # n = 9, m = 6; k = 3
    q = polyder(p)
    psis = norm.(AGCD.sylvester_matrix_singular_values(coeffs(p), coeffs(q)))
    k = 3
    @test (psis[k-1] > 1e-8) && (psis[k] < 1e-8)

    p =  (x-1)^10 * (x-2)^10 * (x-3)^10 # n = 30, m = 27; k = 3
    q = polyder(p)
    k  = 3
    psis = norm.(AGCD.sylvester_matrix_singular_values(coeffs(p), coeffs(q)))
    @test ((psis[1:end-1]./psis[2:end])[k-1] > 10 ) &&  ((psis[1:end-1]./psis[2:end])[k] < 10 ) 

    p = prod((x-i)^i for i in 1:6) # n = 21, u = 15, k=6
    q = polyder(p)
    psis = norm.(AGCD.sylvester_matrix_singular_values(coeffs(p), coeffs(q)))
    k = 6
    @test (((psis[1:end-1]./psis[2:end])[k-1] > 10 ) &&  ((psis[1:end-1]./psis[2:end])[k] < 10 ) ) 

end    

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
    p = (x-1)^2 * (x-2)^2 * (x+3)^2
    u,v,w,err = AGCD.agcd(p)
    @test length(v)-1 > 3   # Fails, tolerances are too strict

    u,v,w,err = AGCD.agcd(p, θ=1e-4, ρ=1e-4)
    @test length(v)-1 == 3

end
    
    
