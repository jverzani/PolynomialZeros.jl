module PolynomialZeros


using Polynomials
import PolynomialRoots
using PolynomialFactors
using Compat

export poly_zeros
export agcd, multroot
export Over 


include("utils.jl")
include("over.jl")
include("agcd.jl")
include("multroot.jl")
include("real_roots.jl")
#include("amvw.jl")






immutable RingType{T} end


"""

`poly_zeros(f, domain)`: Find zeros of the polynomial `f` within the specified domain

* `f` can be an instance of `Poly` type (from `Polynomials.jl`) or a callable object which can be so converted. Will throw a `DomainError` if the object `f` can not be converted to `Poly{T}`.        

* `domain` is one of

    - `over.C` (the default) for solving over complex values (`Complex{Float64}`). Use `Over.CC{T}` to specfy a type `T<: AbstractFloat` other than `Float64`

    - `over.R` for solving over the real line (`Float64`). Use `Over.RR{T}` to specify a `T <: Integer` other than `Float64`

    - `over.Q` for solving over the rational (`Rational{Int}`). Use `Over.RR{T}` to specify a `T <: Integer` other than `Int`

    - `over.Z` for solving over the integers (`Int`). Use `Over.ZZ{T}` to specify a `T` other than `Int`

    - `over.Zp{p}` for solving over the finite field `ZZ_p`, p a prime

Returns an array of zeros, possible empty. May throw error if polynomial type is appropriate for specified domain.


    
"""
poly_zeros(f) = poly_zeros(f, Over.CC{Float64}) # default
poly_zeros(f, ::Type{Over.C}) = poly_zeros(f, Over.CC{Float64})
function poly_zeros{T<:Float64}(f, ::Type{Over.CC{T}})
    
    p = as_poly(Complex{T}, f)
    PolynomialRoots.roots(p.a, polish=true)

end

function poly_zeros{T<:BigFloat}(f, ::Type{Over.CC{T}})
    
    p = as_poly(Complex{T}, f)
    PolynomialRoots.roots(p.a)

end



poly_zeros(f, ::Type{Over.R}) = poly_zeros(f, Over.RR{Float64})
function poly_zeros{T <: AbstractFloat}(f, ::Type{Over.RR{T}})
    p = as_poly(T, f)
    real_roots(p)
end

                                                                                      
poly_zeros(f, ::Type{Over.Q}) = poly_zeros(f, Over.QQ{Int64})
function poly_zeros{T <: Integer}(f, ::Type{Over.QQ{T}})
    
    p = as_poly(Rational{T}, f)
    PolynomialFactors.rational_roots(p)
    
end


poly_zeros(f, ::Type{Over.Z}) = poly_zeros(f, Over.ZZ{Int64})
function poly_zeros{T <: Integer}(f, ::Type{Over.ZZ{T}})

    p = as_poly(T, f)
    rts = PolynomialFactors.rational_roots(p)
    T[convert(T, r) for r in  filter(u->u.den==one(T), rts)]
end



function poly_zeros{q <: Int}(f, ::Type{Over.Zp{q}})
    # error if q is not prime?
    
    p = as_poly(Int, f)
    fs = PolynomialFactors.factormod(p,q)
    ls = filter((r,n) -> degree(r) == 1, fs)
    [mod(-r[0] * invmod(r[1],q), q) for (r,n) in ls]
    
end


end # module
