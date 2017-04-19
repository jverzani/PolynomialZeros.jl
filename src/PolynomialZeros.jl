__precompile__()
module PolynomialZeros


using Polynomials
#import PolynomialRoots
using PolynomialFactors
using Compat

export poly_zeros
export Over 

export agcd, multroot


include("utils.jl")
include("over.jl")
include("special_cases.jl")
include("agcd.jl")
include("multroot.jl")
include("real_roots.jl")
#include("amvw.jl")








"""

`poly_zeros(f, domain)`: Find zeros of the polynomial `f` within the specified domain

* `f` can be an instance of `Poly` type (from `Polynomials.jl`) or a callable object which can be so converted. Will throw a `DomainError` if the object `f` can not be converted to `Poly{T}`.        

* `domain` is one of

    - `over.C` (the default) for solving over complex values (`Complex{Float64}`). Use `Over.CC{T}` to specfy a type `T<: AbstractFloat` other than `Float64`

    - `over.R` for solving over the real line (`Float64`). Use
      `Over.RR{T}` to specify a `T <: Integer` other than
      `Float64`. The algorithm assume the polynomial is square free
      (none of its factors are squares over R). This is important for
      floating point coefficients. Pass the argument
      `square_free=false` to have an *approximate* gcd used to create
      a square free version.

    - `over.Q` for solving over the rational (`Rational{Int}`). Use `Over.RR{T}` to specify a `T <: Integer` other than `Int`

    - `over.Z` for solving over the integers (`Int`). Use `Over.ZZ{T}` to specify a `T` other than `Int`

    - `over.Zp{p}` for solving over the finite field `ZZ_p`, p a prime

Returns an array of zeros, possible empty. May throw error if polynomial type is appropriate for specified domain.


    
"""
poly_zeros(f) = poly_zeros(f, Over.CC{Float64}) # default
poly_zeros(f, ::Type{Over.C}) = poly_zeros(f, Over.CC{Float64})
function poly_zeros{T<:Float64}(f, U::Type{Over.CC{T}})
    
    p = as_poly(Complex{T}, f)
    
    fn = special_case(p.a, U)
    if fn == identity
        PolynomialRoots.roots(p.a, polish=true)
    else
        fn(p.a, U)
    end

end

function poly_zeros{T<:BigFloat}(f, U::Type{Over.CC{T}})
    
    p = as_poly(Complex{T}, f)
    fn = special_case(p.a, U)
    if fn == identity
        PolynomialRoots.roots(p.a)
    else
        fn(p.a, U)
    end
end



poly_zeros(f, ::Type{Over.R};square_free=true) = poly_zeros(f, Over.RR{Float64}, square_free=square_free)
function poly_zeros{T <: AbstractFloat}(f, U::Type{Over.RR{T}}; square_free=true)
    p = as_poly(T, f)
    fn = special_case(p.a, U)
    if fn == identity
        if square_free
            real_roots_sqfree(p)
        else
            real_roots(p)
        end
    else
        fn(p.a, U)
    end
    
end

                                                                                      
poly_zeros(f, ::Type{Over.Q}) = poly_zeros(f, Over.QQ{Int})
function poly_zeros{T <: Integer}(f, U::Type{Over.QQ{T}})
    
    p = as_poly(Rational{T}, f)
    fn = special_case(p.a, U)
    if fn == identity
        d = PolynomialFactors.factor(p)
        d = filter((k,v) -> degree(k) == 1, d)
        vcat([ones(Rational{T},v)*(-p[0]//p[1]) for (p,v) in d]...)
    else
        fn(p.a, U)
    end
    
end


poly_zeros(f, ::Type{Over.Z}) = poly_zeros(f, Over.ZZ{Int64})
function poly_zeros{T <: Integer}(f, U::Type{Over.ZZ{T}})

    p = as_poly(T, f)

    fn = special_case(p.a, U)
    if fn == identity
        d = PolynomialFactors.factor(p)
        d = filter((k,v) -> degree(k) == 1, d)
        d = filter((k,v) -> rem(k[0], k[1]) == 0, d)
        vcat([ones(Rational{T},v)*(-div(p[0], p[1])) for (p,v) in d]...)
    else
        fn(p.a, U)
    end
end



function poly_zeros{q}(f, U::Type{Over.Zp{q}})
    # error if q is not prime?
    
    p = as_poly(BigInt, f)

    fn = special_case(p.a, U)
    if fn == identity
        fs = PolynomialFactors.factormod(p,q)
        ls = filter((r,n) -> degree(r) == 1, fs)
        [mod(-r[0] * invmod(r[1],q), q) for (r,n) in ls]
    else
        fn(p.a, U)
    end
    
end


end # module
