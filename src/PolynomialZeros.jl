__precompile__()
module PolynomialZeros

using Polynomials
import PolynomialRoots

using Iterators
using Combinatorics
using PolynomialFactors

#using AMVW
using Compat


export poly_roots
export Over 

# export agcd, multroot # qualify these





include("utils.jl")
include("over.jl")
include("special_cases.jl")

include("agcd/agcd.jl")        # AGCD.agcd
include("agcd/multroot.jl")    # MultRoot.multroot

include("RealRoots/real_roots.jl")  # RealRoots.real_roots(p::Poly)


using .AGCD
using .MultRoot


# Separate pakage?
#import AMVW
#include("AMVW/AMVW.jl")        # AMVW.amvw or AMVW.amvw_pencil
#using .AMVW



"""

`poly_roots(f, domain)`: Find zeros of the polynomial `f` within the specified domain.

* `f` can be an instance of `Poly` type (from `Polynomials.jl`) or a callable object which can be so converted. Will throw a `DomainError` if the object `f` can not be converted to `Poly{T}`.        

* `domain` is one of

    - `over.C` (the default) for solving over complex values
      (`Complex{Float64}`). Use `Over.CC{T}` to specfy a type `T<:
      AbstractFloat` other than `Float64`. The default method is from
      `PolynomialRoots`. Pass the argument `method=:Roots` to use the
      `roots` function from `Polynomials.jl`. For a degree n
      polynomial over C, all n roots should be returned (including
      multiplicities).

    - `over.R` for solving over the real line (`Float64`). Use
      `Over.RR{T}` to specify a `T <: Integer` other than
      `Float64`. The algorithm assumes the polynomial is square free
      (none of its factors are squares over R). This is important for
      floating point coefficients. Pass the argument
      `square_free=false` to have an *approximate* gcd used to create
      a square-free version. Only unique real roots are returned (no
      multiplicities).

    - `over.Q` for solving over the rational numbers
      (`Rational{Int}`). Use `Over.RR{T}` to specify a `T <: Integer`
      other than `Int`. Only unique rational roots are returned (no
      multiplicities).

    - `over.Z` for solving over the integers (`Int`). Use `Over.ZZ{T}`
      to specify a `T` other than `Int`. Only unique integer roots are
      returned (no multiplicities).

    - `over.Zp{p}` for solving over the finite field `ZZ_p`, p a
      prime.  Only unique roots are returned (no multiplicities).

Returns an array of zeros, possibly empty. May throw error if polynomial type is inappropriate for specified domain.

Examples:

```julia
using Polynomials, PolynomialZeros
x = variable(); p = x^5 - x - 1
poly_roots(p, Over.C)  # 5
poly_roots(p, Over.R)  # 1
poly_roots(p, Over.Q)  # empty

p = x^3 - 1
poly_roots(p, Over.C)  # 3
poly_roots(p, Over.R)  # 1
poly_roots(p, Over.Q)  # 1
poly_roots(p, Over.Z)  # 1
poly_roots(p, Over.Zp{7})  # 3
```
    
"""
function poly_roots(f; method=:PolynomialRoots)
    T = eltype(float(zero(e_type(f))))
    poly_roots(f, Over.CC{T}; method=method) # default
end
poly_roots(f, ::Type{Over.C}; method=:PolynomialRoots) = poly_roots(f, Over.CC{Float64}; method=method)
function poly_roots{T<:AbstractFloat}(f, U::Type{Over.CC{T}}; method=:PolynomialRoots)


    ps = poly_coeffs(T, f)
    fn = special_case(ps, U)
    
    if fn == identity
        if method == :PolynomialRoots
            PolynomialRoots.roots(ps, polish=true)
#        elseif method == :AMVW
#            AMVW.poly_roots(ps)
#        elseif method ==:AMVW_Pencil
#            AMVW.poly_roots(ps)
        else
            convert(Vector{Complex{Float64}}, roots(Poly(ps)))
        end
    else
        fn(ps, U)
    end

end



function poly_roots(f, ::Type{Over.R};square_free=true)
    T = eltype(float(zero(e_type(f))))
    poly_roots(f, Over.RR{T}, square_free=square_free)
end
function poly_roots{T <: Real}(f, U::Type{Over.RR{T}}; square_free::Bool=true)

    ps = convert(Vector{T},poly_coeffs(f))
    fn = special_case(ps, U)
    
    if fn == identity
        RealRoots.real_roots(as_poly(T, f), square_free=square_free)
    else
        fn(ps, U)
    end
    
end
# should I do an alias? Best to add if requestd, keeping name space light for now
## realroots(f; kwargs...) = poly_roots(f, Over.R; kwargs...) #real_roots?
                                                                                      
function poly_roots(f, ::Type{Over.Q})
    T = e_type(f)
    poly_roots(f, Over.QQ{T})
end
function poly_roots{T <: Integer}(f, U::Type{Over.QQ{T}})
    
    p = as_poly(Rational{T}, f)
    fn = special_case(poly_coeffs(p), U)
    if fn == identity
        d = PolynomialFactors.factor(p)
        d = filter((k,v) -> degree(k) == 1, d)
        vcat([-p[0]//p[1] for p in keys(d)]...)
    else
        fn(p.a, U)
    end
    
end


function poly_roots(f, ::Type{Over.Z})
    T = eltype(as_poly(f)(0))
    poly_roots(f, Over.ZZ{T})
end

function poly_roots{T <: Integer}(f, U::Type{Over.ZZ{T}})

    p = as_poly(T, f)

    fn = special_case(poly_coeffs(p), U)
    if fn == identity
        d = PolynomialFactors.factor(p)
        d = filter((k,v) -> degree(k) == 1, d)
        d = filter((k,v) -> rem(k[0], k[1]) == 0, d)
        vcat([-div(p[0], p[1]) for (p,v) in d]...)
    else
        fn(p.a, U)
    end
end



function poly_roots{q}(f, U::Type{Over.Zp{q}})
    # error if q is not prime?
    
    p = as_poly(e_type(f), f)

    fn = special_case(poly_coeffs(p), U)
    if fn == identity
        fs = PolynomialFactors.factormod(p,q)
        ls = filter((r,n) -> degree(r) == 1, fs)
        [mod(-r[0] * invmod(r[1],q), q) for (r,n) in ls]
    else
        fn(p.a, U)
    end
    
end


end # module
