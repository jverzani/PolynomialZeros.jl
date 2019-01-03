module MultRoot
using Polynomials
using LinearAlgebra
include("../utils.jl")
## The main function here is `MultRoot.multroot`

## Polynomial root finder for polynomials with multiple roots
##
## Based on "Computing multiple roots of inexact polynomials"
## http://www.neiu.edu/~zzeng/mathcomp/zroot.pdf
## Author: Zhonggang Zeng
## Journal: Math. Comp. 74 (2005), 869-903
##
## Zeng has a MATLAB package `multroot`, from which this name is derived.
## Basic idea is
## 1) for polynomial p we do gcd decomposition p = u * v; p' = u * w. Then roots(v) are the roots without multiplicities.
## 2) can repeat with u to get multiplicities.
##
## This is from Gauss, as explained in paper. Zeng shows how to get u,v,w when the polynomials
## are inexact due to floating point approximations or even model error. This is done in his
## algorithm II.
## 3) Zeng's algorithm I (pejroot) uses the pejorative manifold of Kahan and Gauss-Newton to
## improve the root estimates from algorithm II (roots(v)). The pejorative manifold is defined by
## the multiplicities l and is operationalized in evalG and evalJ from Zeng's paper.

using Polynomials
using ..AGCD
monic(p) = p/p[end]
rcoeffs(p) = reverse(p.a)

## map monic(p) to a point in C^n
## p = 1x^n + a1x^n-1 + ... + an_1 x + an -> (a1,a2,...,an)
function p2a(p::Poly)
    p = monic(p)
    rcoeffs(p)[2:end]
end

## get value of gl(z). From p16
function evalG(zs::Vector, ls::Vector)
    length(zs) == length(ls) || throw("Length mismatch")

    s = prod([poly([z])^l for (z,l) in zip(zs, ls)])  # \prod (x-z_i)^l_i
    p2a(s)
#    rcoeffs(s)[2:end]
end

## get jacobian J_l(z), p16

## get jacobian J_l(z), p16
function evalJ!(J, zs::Vector{T}, ls::Vector) where {T}
    length(zs) == length(ls) || throw("Length mismatch")

    n, m = sum(ls), length(zs)
    u = evalG(zs, ls .- 1)
    pushfirst!(u, one(T))

#    x = Polynomials.variable(T)
#    u = prod((x-z)^(l-1) for (z,l) in zip(zs, ls))  # \prod (x-z_i)^l_i


    for j in 1:m
        s = -ls[j] * u

        for (l, zl) in zip(1:m, zs)
            l == j && continue
            s = AGCD._polymul(s, (one(T), -zl))
        end
    J[:,j] = s#rceffs(s)
    end
    J
end
function evalJ(zs::Vector{T}, ls) where {T}
    n, m = sum(ls), length(zs)
    J = zeros(T, n, m)
    evalJ!(J, zs, ls)
    J
end

## Gauss-Newton iteration to solve weighted least squares problem
## G_l(z) = a, where a is related to monic version of polynomial p
## l is known multiplicity structure of polynomial p = (x-z1)^l1 * (x-z2)^l2 * ... * (x-zn)^ln
## Algorithm I, p17
function pejroot(p::Poly, z0::Vector, l::Vector{Int};
                 wts::Union{Vector, Nothing}=nothing, # weight vector
                 tol = 1e-8,
                 maxsteps = 100
                      )

    a = p2a(p) #rcoeffs(monic(p))[2:end] # an_1, an_2, ..., a2, a1, a0

    if wts == nothing
        wts = map(u -> min(1, 1/abs.(u)), a)
    end

    ## Solve WJ Δz = W(Gl(z) - a) in algorithm I
#    G(z) = (evalG(z, l) - a)
#    update(z, l) = z -  AGCD.weighted_least_square(evalJ(z,l), G(z), wts)

    G = z -> evalG(z,l) - a

    zk = copy(z0);
    zk1 = zk - AGCD.weighted_least_square(evalJ(zk,l), G(zk), wts)
    deltaold = norm(zk1 - zk,2)
    zk = zk1

    cvg = false
    for ctr in 1:maxsteps

        zk1 = zk - AGCD.weighted_least_square(evalJ(zk,l), G(zk), wts)
        delta = norm(zk1 - zk, 2)

        if delta > deltaold
            @info "Growing delta. Best guess is being returned."
            break
        end

        ## add extra abs(delta) < 100*eps() condition
        if delta^2 / (deltaold - delta) < tol || abs(delta) < 100*eps()
            cvg = true
            break
        end

        deltaold = delta
        zk=zk1
    end

    if !cvg @info ("""
Returning the initial estimates, as the
algorithm failed to improve estimates for the roots on the given
pejorative manifold.
""")
        return(z0)
    end
    return(zk1)
end

"""
    multroot(p; [θ, ρ, ϕ, δ])

Find roots of polynomial `p`,

The `multroot` function returns the roots and their multiplicities
for `Poly` objects. It performs better than `roots` when the
polynomial has multiplicities.

Based on "Computing multiple roots of inexact polynomials"
Zhonggang Zeng
Journal: Math. Comp. 74 (2005), 869-903
http://www.neiu.edu/~zzeng/mathcomp/zroot.pdf

Zeng has a MATLAB package `multroot`, from which this name is derived.

Basic idea is:
* for polynomial p we do gcd decomposition p = u * v; p' = u * w. Then roots(v) are the roots without multiplicities.
* can repeat with u to get multiplicities.

The basic idea is from Gauss, as explained in paper. Zeng shows how to get u,v,w when the polynomials
are inexact due to floating point approximations or even model error. This is done in his
algorithm II.

Zeng's algorithm I (pejroot) uses the pejorative manifold of Kahan and Gauss-Newton to
improve the root estimates from algorithm II (roots(v)). The pejorative manifold is defined by
the multiplicities l and is operationalized in `evalG` and `evalJ` from Zeng's paper.

Examples:
```
x = poly([0.0]);
p = (x-1)^4 * (x-2)^3 * (x-3)^3 * (x-4)
multroot(p) # ([4.0, 3.0, 2.0, 1.0], [1, 3, 3, 4])

## For "prettier" printing, results can be coerced to a dict
Dict(k=>v for (k, v) in zip(multroot(p)...))
## Dict{Float64,Int64} with 4 entries:
##   4.0 => 1
##   2.0 => 3
##   1.0 => 4
##   3.0 => 3

## compare to
roots(p)
# 11-element Array{Complex{Float64},1}:
#   4.000000000049711 + 0.0im
#   3.000340992482669 + 0.0005901104690775108im
#   3.000340992482669 - 0.0005901104690775108im
#  2.9993180137669726 + 0.0im
#  2.0006129094631833 + 0.0im
#  1.9996935459268157 + 0.0005303204210266187im
#  1.9996935459268157 - 0.0005303204210266187im
#  1.0006651061397798 + 0.00066713538542472im
#  1.0006651061397798 - 0.00066713538542472im
#  0.9993348938107887 + 0.0006630950464016094im
#  0.9993348938107887 - 0.0006630950464016094im

## Large order polynomials prove difficult. We can't match the claims in Zeng's paper
## as we don't get the pejorative manifold structure right.
p = poly(1.0:7.0));
multroot(p^2) ## should be 1,2,3,4,...,7 all with multplicity 2, but
## ([7.00028, 6.99972, 6.00088, 5.99912, 5.00102, 4.99898, 4.00055, 3.99945, 3.00014, 2.99986, 2.00002, 1.99998, 1.0, 0.999999], [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])

## nearby roots can be an issue
delta = 0.00001  ## delta = 0.0001 works as desired.
p = (x-1 - delta)*(x-1)*(x-1 + delta)
multroot(p)
## ([1.0], [3])
```
"""
function multroot(p::Poly;
                  θ::Real=1e-8,  #
                  ρ::Real=1e-10, # initial residual tolerance
                  ϕ::Real=1e2,   # residual tolerance growth factor
                  δ::Real=1e-8   # passed to solve y sigma

                  )

    Polynomials.degree(p) == 0 && error("Degree of `p` must be atleast 1")

    if Polynomials.degree(p) == 1
        return (roots(p), [1])
    end

    p = Poly(float(coeffs(p)))  # floats, not Int

    u_j, v_j, w_j, residual= AGCD.agcd(p, polyder(p), θ=θ,  ρ=ρ)
    ρ = max(ρ, ϕ * residual)

    ## bookkeeping
    zs = roots(v_j)
    ls = ones(Int, length(zs))

    p0 = u_j

    while Polynomials.degree(p0) > 0
        if Polynomials.degree(p0) == 1
            z = roots(p0)[1]
            tmp, ind = findmin(abs.(zs .- z))
            ls[ind] = ls[ind] + 1
            break
        end

        u_j, v_j, w_j, residual= AGCD.agcd(p0, polyder(p0), θ=θ, ρ=ρ)

        ## need to worry about residual between
        ## u0 * v0 - monic(p0) and u0 * w0 - monic(Polynomials.polyder(p0))
        ## resiudal tolerance grows with m, here it depends on
        ## initial value and previous residual times a growth tolerance, ϕ
        ρ = max(ρ, ϕ * residual)

        ## update multiplicities
        for z in roots(v_j)
            tmp, ind = findmin(abs.(zs .- z))
            ls[ind] = ls[ind] + 1
        end

        ## rename
        p0 = u_j
    end


    if maximum(ls) == 1
        return (zs, ls)
    else
        zs = pejroot(p, zs, ls)
        return (zs, ls)
    end
end

## Different interfaces

## can pass in vector too
multroot(p::Vector{T}; kwargs...) where {T <: Real} = multroot(Poly(p); kwargs...)

## Can pass in function
function multroot(f::Function; kwargs...)
    p = as_poly(Float64, f)
    multroot(p; kwargs...)

end


end
