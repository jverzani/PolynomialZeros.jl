function deflate_leading_zeros{T}(ps::Vector{T})
    ## trim any 0s from the end of ps
    N = findlast(!iszero, ps)
    K = findfirst(!iszero, ps)

    N == 0 && return(zeros(T,0), length(ps))
    ps = ps[K:N]
    ps, K-1
end

## take poly [p0, p1, ..., pn] and return
## [q_m-1, q_m-2, ..., q0], k
## where we trim of k roots of 0, and then make p monic, then reverese
## monomial x^5
function reverse_poly{T}(ps::Vector{T})
    # assume we have called deflate_leading_zeros
    qs = reverse(ps./ps[end])[2:end]
    qs
end

#
function quadratic_equation{T <: Real}(a::T, b::T, c::T)   
    qdrtc(a, -(0.5)*b, c)
end

## make more robust
function quadratic_equation{T}(a::Complex{T}, b::Complex{T}, c::Complex{T})
    d = sqrt(b^2 - 4*a*c)
    e1 = (-b + d)/(2a); e2 = (-b-d)/(2a)
    return (real(e1), imag(e1), real(e2), imag(e2))
    
end

## Kahan quadratic equation with fma
##  https://people.eecs.berkeley.edu/~wkahan/Qdrtcs.pdf

## solve ax^2 - 2bx + c
function qdrtc{T <: Real}(a::T, b::T, c::T)
    # z1, z2 roots of ax^2 - 2bx + c
    d = discr(a,b,c)  # (b^2 - a*c), as 2 removes 4
    
    if d <= 0
        r = b/a  # real
        s = sqrt(-d)/a #imag
        return (r,s,r,-s)
    else
        r = sqrt(d) * (sign(b) + iszero(b)) + b
        return (r/a, zero(T), c/r, zero(T))
    end
end

## more work could be done here.
function discr{T}(a::T,b::T,c::T)
    pie = 3.0 # depends on 53 or 64 bit...
    d = b*b - a*c
    e = b*b + a*c

    pie*abs(d) > e && return d

    p = b*b
    dp = muladd(b,b,-p)
    q = a*c
    dq = muladd(a,c,-q)

    (p-q) + (dp - dq)
end

##
## solve degree 2 or less case
function solve_simple_cases{T}(ps::Vector{T})
    S = T <: Complex ? T : Complex{T}
    
    N = length(ps)
    
    if N <= 1
        return S[]
    elseif N == 2
        return S[-ps[1]/ps[2]]
    elseif N == 3
        c,b,a = ps
        r1,i1, r2,i2 = quadratic_equation(a,b,c)
        S[complex(r1, i1), complex(r2, i2)]
    end
end

