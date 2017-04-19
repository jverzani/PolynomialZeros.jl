using PolynomialZeros
using Polynomials
using Base.Test

# write your own tests here
function Wilkinson(n, T=Float64)
    x = variable(T)
    prod(x-i for i in 1:n)
end

function Chebyshev(n, T=Float64)
    x = variable(T)
    n == 0 && return one(x)
    n == 1 && return x
    2 * x * Chebyshev(n-1, T) - Chebyshev(n-2, T)
end

## only a few
function Cyclotomic(n, T)
    x = variable(T)
    n == 1 && return x-1
    n == 2 && return x + 1
    n == 5 && return x^4 + x^3 + x^2 + x + 1
    n == 10 && return x^4 - x^3 + x^2 -x + 1
    n == 20 && return x^8 - x^6 + x^4 - x^2 + 1
    throw(DomainError)
end
    


@testset "Over.C" begin
    @test length(poly_zeros(Wilkinson(10), Over.C)) == 10
    @test length(poly_zeros(x -> x^5 - x -1, Over.C)) == 5
    @test maximum(norm.(poly_zeros(Chebyshev(5), Over.C))) <= 1
end

@testset "Over.R" begin
#    @test length(poly_zeros(Wilkinson(10), Over.R)) == 10
    @test length(poly_zeros(x -> x^5 - x -1, Over.R)) == 1
    @test maximum(norm.(poly_zeros(Chebyshev(5), Over.R))) <= 1
end 

@testset "Over.Q" begin
    @test length(poly_zeros(Wilkinson(5, Int), Over.Q)) == 5    
    @test length(poly_zeros(x -> x^5 - x -1, Over.Q)) == 0
end


@testset "Over.Zp{q}" begin
    p = x -> x^8 - 1
    @test length(poly_zeros(p, Over.Zp{7})) == 2
    @test length(poly_zeros(p, Over.Zp{17})) == 8
end


@testset "special cases" begin
    x = variable()
    p = (x-1)*(x^2 + 1)
    rts_c = poly_zeros(p, Over.C)
    rts_r = poly_zeros(p, Over.R)
    @test length(rts_c) == 3
    @test length(rts_r) == 1    

    x = variable(Int)
    p = (x-1)*(2x-3)
    rts_q = poly_zeros(p, Over.Q)
    rts_z = poly_zeros(p, Over.Z)
    @test length(rts_q) == 2
    @test length(rts_z) == 1
end
