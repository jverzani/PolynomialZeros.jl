module AGCD
using Polynomials
using LinearAlgebra

## This provide AGCD.agcd for finding an *approximate* GCD of two polynomials. The most common use in this package
## is to reduce a polynomial `p` to a square free polynomial `q=p/gcd(p, p')`.


## This isn't too accurate for higher order polys, but
## seems to be okay now for degree 20 or less
## Can we do better?
## Follow up on these:
# * [Winkler&Hasan](http://ac.els-cdn.com/S0377042712003135/1-s2.0-S0377042712003135-main.pdf?_tid=126266d0-ea31-11e6-a1a7-00000aab0f6b&acdnat=1486140896_2119601f37cf58ff0d3eefbc05fecbb4)
# * [LiLiuZhi](http://www.sciencedirect.com/science/article/pii/S0377042707000271?np=y&npKey=f5e3b7a1bc9583214348aad066170ea389212ebba443454938cacfdbf96177f8)
# * [Winkler](http://ml.dcs.shef.ac.uk/summer_school/ln.pdf)
# * [Matonardi&...](http://www.sciencedirect.com/science/article/pii/S0024379515005649)
# * [FarmerPhD](http://www.learninglink.bbk.ac.uk/site/assets/files/1025/farmer.pdf)
# * [WinklerLaoHasan](http://www.sciencedirect.com/science/article/pii/S0377042712000751)
# * [Sagraloff&...](https://arxiv.org/pdf/1308.4088v2.pdf)
# * [](https://arxiv.org/pdf/1207.0630v3.pdf)

## ------------------------

# we use vectors, not polynomials. Some of this if copied and modified from Polynomials.jl
# originally by @vtjnash
_degree(p::Vector) = length(p) - 1
function _polyval(p::Vector{T}, x::S) where {T, S}
    R = promote_type(T,S)
    y = convert(R, p[end])
    for i in lastindex(p)-1:-1:1
        y = p[i] + x * y
    end
    y
end
function _monic!(ps::Vector{T}) where T
    lambda = iszero(ps[end]) ? one(T) : 1 / ps[end]
    ps .= ps .* lambda
end
function _monic(p::Vector)
    p1 = copy(p)
    _monic!(p1)
    p1
end
_float(p::Vector) = float(p)


function _polymul(p::Vector{T}, q) where {T}
    R = promote_type(T,eltype(q))
    n = length(p)-1
    m = length(q)-1
    a = zeros(R, m+n+1)

    for i = 1:n+1
        for j = 1:m+1
            a[i+j-1] += p[i] * q[j]
        end
    end
    a
end

function _polyder(p::Vector{T}) where {T}
    S = eltype(float(p[1]))
    n = length(p)
    a2 = Vector{S}(undef, n-1)
    for i = 1:n-1
        a2[i] = p[i+1] * i
    end
    a2
end

## ------------------------

## Special matrices
## p = [p_0, p_1, ..., p_n]
function cauchy_matrix_size(p::Vector, k::Integer)
    n = length(p)
    (n + k - 1, k)
end
function cauchy_matrix!(M, p::Vector{T}, k::Integer) where {T}
    n = length(p)
    for i in 1:k
        M[(1:n) .+ (i-1), i] = p
    end
end
function cauchy_matrix(p::Vector{T}, k::Integer) where {T}
    n,m = cauchy_matrix_size(p, k)
    out = zeros(T, n, m)
    cauchy_matrix!(out, p, k)
    out
end

## Jacobian F(u,v,w) = [p,p'] is J(u,v,w)
function JF_size(u, v, w)

    m, k, j = _degree(u), _degree(v), _degree(w)
    n, l = m + k, m + j

    ai, aj = cauchy_matrix_size(v, m + 1)
    bi, bj = cauchy_matrix_size(u, k + 1)
    di, dj = cauchy_matrix_size(w, m + 1)
    fi, fj = cauchy_matrix_size(u, j + 1)
    ci, cj = ai, fj
    ei, ej = di, bj




    (1 + ai + di, aj + bj + cj)



end

function JF!(M, u::Vector{T}, v, w) where {T}
    m, k, j = _degree(u), _degree(v), _degree(w)
    n, l = m + k, m + j

    # JF_size should return these
    ai, aj = cauchy_matrix_size(v, m + 1)
    bi, bj = cauchy_matrix_size(u, k + 1)
    di, dj = cauchy_matrix_size(w, m + 1)
    fi, fj = cauchy_matrix_size(u, j + 1)
    ci, cj = ai, fj
    ei, ej = di, bj


    M[1,1] = one(T)
    cauchy_matrix!(view(M, 1 .+ (1:ai), 1:aj),        v, m + 1)
    cauchy_matrix!(view(M, 1 .+ (1:bi), aj .+ (1:bj)), u, k + 1)
#    M[1 + (1:ci), aj + bj + (1:cj)] = zeros(T, n+1, m+1-j)
    cauchy_matrix!(view(M, 1 + ai .+ (1:di), 1:dj),   w, m + 1)
#    M[1 + bi + (1:ei), dj + (1:ej)] = zeros(T, m + 1, n + 1 - j)
    cauchy_matrix!(view(M, 1 + ci .+ (1:fi), dj + ej .+ (1:fj)), u, j + 1)


end
function JF(u::Vector{U}, v::Vector{V}, w::Vector{W}) where {U,V,W}
    R = promote_type(U,V, W)
    M = zeros(R, JF_size(u, v, w)...)
    JF!(M, u, v, w)
    M
end


# eqn (12) to get weights
function _weights(p::Vector{T},q) where {T}
    n, m = length(p), length(q)
    wts = ones(T, 1 + n + m)
    for (i,pj) in enumerate(p)
        wts[1+i] = min(one(real(T)), inv(abs(pj)))
    end
    for (i,qj) in enumerate(q)
        wts[1 + n + i] = min(one(real(T)), inv(abs(qj)))
    end
    wts
end

# reduce residual by Gauss-Newton algorithm
# updates u,v,w; returns residual error and flag for convergence
function reduce_residual!(u,v,w, p::Vector{T}, q, wts, ρ) where {T}
    m, n, l = map(_degree, (u, v, w))
    # preallocate
    A = zeros(T, JF_size(u, v, w)...)
    b = zeros(T, 1 + length(p) + length(q))
    up,vp,wp = copy(u),copy(v),copy(w)

    inc = zeros(T, m + n + l + 3) #weighted_least_square(A, b, wts)

    ρm_1, ρm = Inf, Inf
    MAXSTEPS = 5
    ctr = 0
    flag = :not_converged



    while ctr <= MAXSTEPS
        ρm =  agcd_update!(p, q, u,v,w, m, n, A, b, inc, up, vp, wp, wts, ρm_1)

        if ρm > 1.5 * ρm_1
            return (ρm, flag)
        end

        # threshold here is really important, but not clear what it should be
        # we empirically get better results---it seems---with sqrt(norm(p,2))
        # though recomendation if norm(u,2); Might also make sense to include
        # q in the computation
        if ρm <= ρ * norm(u,2)
            flag = :converged
            return (ρm, flag)
        end
        ctr += 1
        ρm_1 = ρm
    end

    ρm, flag
end



## Lemma 4.1, (25)
## solve weighted least squares problem W*(Ax - b) = 0
function weighted_least_square(A, b, w)
    W = diagm(0 => w)
    (W * A) \ (W * b)
end

function weighted_least_square!(M, A, b, w)
    W = diagm(0 => w)
    M[:] .= (W * A) \ (W * b)
end



## compute F(u,v,w) - [p, p'] = [u*v, u*w] - [p, p']
Fmp(p,q,u,v,w) =  [u[end]; _polymul(u,v); _polymul(u,w)] .- [1; p; q]
function Fmp!(b, p,q,u,v,w)
    b[1] = u[end] - 1
    for (i,v) in enumerate(_polymul(u,v) .- p)
        b[i] = v
    end
    offset = 1 + length(p)
    for (i,v1) in enumerate(_polymul(u,w) .- q)
        b[offset + i] = v1
    end
end

residual_error(p,q,u,v,w, wts=ones(1 + length(p) + length(q))) = norm( ([_polymul(u,v); _polymul(u,w)] .- [p; q]) .* wts[2:end], 2)

### refine u,v,w estimate using
## Newton's method and weighted least squares
# updates in place u,v,w
function agcd_update!(p, q, u,v,w,m,n, A, b, inc, up, vp, wp, wts, err0)

    JF!(A, u,v,w)
    Fmp!(b, p,q,u,v,w)
    weighted_least_square!(inc, A, b, wts)

    Δu = inc[1:(1+m)]
    Δv = inc[(m+2):(m+n+2)]
    Δw = inc[(m+n+3):end]

    up .- u; vp .= v; wp .= w

    up .-=  Δu; _monic!(up)
    vp .-=  Δv; _monic!(vp)
    wp .-=  Δw; _monic!(wp)

    err = residual_error(p,q,up,vp,wp, wts)

    # update if we improve(ish)
    if err  < 1.1*err0
        u .= up; v .= vp; w .= wp
    end


    # return error estimate


    err
end




function _mult(Gs, A)
    for G in Gs
        A = G * A
    end
    A
end

function qr_sylvester!(Gs, A0, k)

        A = copy(A0)
        m,n = size(A)
        for i in m:-1:2
            Gi, r = givens(A[1,1],A[i,1], 1, i)
            A = Gi * A
            push!(Gs, Gi)  # [G1, G2, ..., Gn]
        end

        for i in m:-1:3
            Gi, r = givens(A[2,2],A[i,2], 2, i)
            A = Gi * A
            push!(Gs, Gi)
        end
        return A[1:2,1:2]
end

function qr_sylvester!(Gs, A0, k, R)
    T = eltype(A0)
    N = (k-1)
    Zs = zeros(T, N)
    B = vcat(hcat(Zs, Zs), A0)
    B = _mult(Gs, B)
    m,n = size(B)
    for i in m:-1:(2k)
        Gi, r = givens(B[2(k-1)+1, 1],   B[i, 1], 2(k-1)+1, i)
        B = Gi * B
        push!(Gs, Gi)
    end
    for i in m:-1:(2k+1)
            Gi, r = givens(B[2k, 2], B[i, 2], 2k, i)
        B = Gi * B
        push!(Gs, Gi)
    end
    return [vcat(R, zeros(T, 2, 2*(k-1))) B[1:(2k),:]]
end


# return sigma, x
# converge on smallest eigenvalue of (A'*A) using power rule
# A=QR; (A'*A)^-1) = (R'*Q'*Q*R)^(-1) = (R'*R)^(-1)
# instead of computing x_{i+1} = (R'*R)^{-1} xi; we solve R'y=xi; Rz=y;x=z/|z|
function smallest_eigval(R::LinearAlgebra.UpperTriangular{T}, thresh=sqrt(eps())) where {T}

    if iszero(det(R))
        return (:iszero, zero(T), T[])
    end

    m,n = size(R)

    x = ones(T, n)
    y = zeros(T, m)
    z = zeros(T, n)
    σ, σ1 = Inf*one(real(T)), Inf*one(real(T))

    flag = :ispositive
    for cnt in 1:10
        y .= R' \ x
        z .= R  \ y
        nz = 1/norm(z,2)
        x .= z .* nz
        sigma = norm(R * x, 2)
        σ1 = abs(sigma)

        if σ1 < 1/10 * σ
            σ = σ1
            continue
        end

        if σ1 < thresh
            flag = :ispossible
            break
        end
        if  abs(σ - σ1)  < 1.1 * σ1
            flag = :ispossible
            break
        end

        σ = σ1
    end

    return (flag, σ1, x)
end


"""

    `agcd(p, q, θ=sqrt(eps()), ρ=1e-2*θ)`

Find an approximate GCD for polynomials `p` and `q` using an algorithm of [Zeng](https://doi.org/10.1090/S0025-5718-04-01692-8).


Returns u,v,w, err where:

* `u*v ≈ monic(p)`
* `u*w ≈ monic(q)`
* The total residual error in these approximations is bounded by `err`.

Further,
* `v` and `w` should share no common roots (`u` is a gcd of `u*v` and `u*w`)
* the roots of `v` should exhaust unique values of roots of `p`.

(If `p` and `q` are specified as vectors, the returned values will be
vectors. If specified as objects of `Polynomial` type then the
returned values will be as well.)

The tolerances are:

* θ: singular value threshold. Used to identify if smallest singular value of Sylvester matrix is ≈ 0
* ρ: initial residual tolerance. If we can get (u,v,w) error less than this, we stop

The algorithm looks for the first `k` for which the corresponding
Sylvester matrix is rank deficient. This follows Lemma 2.4 of the paper, which
finds the smallest singular value. In the process, an estimate for
`u`,`v`, and `w` is produced. If this singular value is approximately
0 (as determined with the parameter θ) *and* the estimates of the
decomposition can be refined via a Gauss-Jordan process to a close
enough approximation (as determined with the parameter ρ) then the
desired values are found.

The act of identifying the value of `k` depends on determining if the
smallest singular value is 0 which in practice is quite hard. The
simple threshold of `norm(p,2) ⋅ θ` will fail for higher-degree
polynomials. For manually identifying this value, the complete set of
estimated singular values for different values of `k` are returned by
`sylvester_matrix_singular_values`. For a given `k` the refined values
of `u`, `v`, and `w` are returned by `rank_k_agcd`.


"""
function agcd(ps::Vector{T}, qs::Vector{S}=_polyder(ps);
              θ=sqrt(eps()),  # no T dependence
              ρ=1e-2*θ, maxk=length(qs)) where {T,S}

    _monic!(ps); _monic!(qs)
    n, m = length(ps), length(qs)
    wts = _weights(ps,qs)

    R = promote_type(T,S)
    A0::Matrix{T} = R[ps vcat(zeros(T, n-m),qs)]

    nm = norm(ps, 2)
    nm = sqrt(nm)  # paper has theta ||f|!_2; but we  we use (||f||_2)^(1/2)
    thresh = nm * θ ## this is sensitive

    Gs = LinearAlgebra.Givens{T}[]
    k = 1
    R = qr_sylvester!(Gs, A0, k)

    local x::Vector{T}

    while k <=  maxk
        V = UpperTriangular(R)
        flag, sigma, x = smallest_eigval(V)

        if flag != :iszero &&  sigma < thresh

            v = x[2:2:end]; _monic!(v)
            w = x[(1+2(n-m)):2:end]; _monic!(w)
            A =  cauchy_matrix(v, n - length(v) + 1)
            u = A \ ps
            _monic!(u)


            ρm, flag = reduce_residual!(u, v, w, ps, qs, wts, ρ)

            if flag == :converged
                return (u, v, w, ρm)
            end

        end

        if flag == :iszero

            u = qs
            w = ones(T, 1)
            v = cauchy_matrix(u, length(ps) - m + 1) \ ps
            _monic!(v)

            return (u,v,w, zero(T))

        end

        # else bumpup k
        k += 1
        if k > m
            u = ones(T,1)
            w = qs
            v = ps
            return (u, v, w, zero(T))
        end
        R =  qr_sylvester!(Gs, A0, k, R)
    end

    # Okay, we gave it a good shot, but didn't find a perfect candidate
    # let's move on...
    v = x[2:2:end]; _monic!(v)
    w = x[3:2:end]; _monic!(w)
    A =  cauchy_matrix(v, n - length(v) + 1)
    u = A \ ps
    _monic!(u)

    wts = _weights(ps,qs)
    ρm, flag = reduce_residual!(u,v,w, ps, qs, wts, ρ)
    return (u, v, w, ρm)

end







# for polys, we precondition

monic(p) = p/p[end]
_float(p::Poly{T}) where {T <: AbstractFloat} = p
_float(p::Poly{T}) where {T} = Poly(float.(p.a), p.var)

## --------------
## preconditioning code
## taken from https://who.rocq.inria.fr/Jan.Elias/pdf/je_masterthesis.pdf
function geometric_mean(a::Vector{T})::T where {T}
    b = filter(!iszero, a)
    n = length(b)
    prod(abs(u)^(one(T)/n) for u in b)
end

function ratio(p::Poly,q::Poly, atol=Base.eps(), rtol=Base.eps())
    as = norm.(filter(!iszero, p.a))
    length(as) == 0 && return Inf

    bs = norm.(filter(!iszero, q.a))
    length(bs) == 0 && return Inf

    max(maximum(as), maximum(bs)) / min(minimum(as), minimum(bs))
end

## find alpha, gamma that minimize ratio of max coefficients to min coefficients, for getting zeros
## 1.12 writes this as a linear programming problem, we just ...
function precondition(p::Poly{T}, q::Poly) where {T}

    m = ratio(p,q)

    alphas = [(2*one(T))^i for i in -5:5]
    phis = [(2*one(T))^i for i in -5:5]

    out = ones(eltype(p),2)

    for α in alphas
        for ϕ in phis
            r = ratio(polyval(p, ϕ * variable(p)), α * polyval(q, ϕ * variable(q)))
            if r < m
                out = [α, ϕ]
            end
        end
    end

    α, ϕ = out

    p = polyval(p, ϕ * variable(p))
    q = α * polyval(q, ϕ * variable(q))

    p = p * (1/geometric_mean(coeffs(p)))
    q = q * (1/geometric_mean(coeffs(q)))

    p, q, ϕ, α

end


##

monic(p) = p/p[end]
_float(p::Poly{T}) where {T <: AbstractFloat} = p
_float(p::Poly{T}) where {T} = Poly(float.(p.a), p.var)

## --------------
## preconditioning code
## taken from https://who.rocq.inria.fr/Jan.Elias/pdf/je_masterthesis.pdf

## compute q(x) = p(phi*x)
function Tphi(ps::Vector{T}, phi) where {T}
    ps .* (phi^(i-1) for i in eachindex(ps))
end
## compute q(x) = p(x-alpha)
function T_alpha(ps::Vector{T}, alpha) where {T}
end

function geometric_mean(a::Vector{T})::T where {T}
    b = filter(!iszero, a)
    n = length(b)
    prod(abs(u)^(one(T)/n) for u in b)
end

function ratio(p::Vector,q::Vector, atol=Base.eps(), rtol=Base.eps())
    as = norm.(filter(!iszero, p))
    length(as) == 0 && return Inf

    bs = norm.(filter(!iszero, q))
    length(bs) == 0 && return Inf

    max(maximum(as), maximum(bs)) / min(minimum(as), minimum(bs))
end

## find alpha, gamma that minimize ratio of max coefficients to min coefficients, for getting zeros
## 1.12 writes this as a linear programming problem, we just ...
function precondition(p::Vector{T}, q::Vector) where {T}

    m = ratio(p,q)

    alphas = [(2*one(T))^i for i in -5:5]
    phis = [(2*one(T))^i for i in -5:5]

    out = ones(eltype(p),2)

    for α in alphas
        for ϕ in phis
            r = ratio(Tphi(p, ϕ), α * Tphi(q, ϕ))
            if r < m
                out = [α, ϕ]
            end
        end
    end

    α, ϕ = out



    p = Tphi(p, ϕ)
    q = α * Tphi(q, ϕ)

    p = p * (1/geometric_mean(p))
    q = q * (1/geometric_mean(q))

    p, q, ϕ, α

end


## -----------------

function agcd(p0::Poly{T}, q0::Poly{S}=polyder(p0); kwargs...) where {T <: Number,S <: Number}

    p, q = _float(p0), _float(q0)

    p, q, phi, alpha = precondition(p, q)

    u0,v0,w0,err = agcd(coeffs(p), coeffs(q); kwargs...)
    u,v,w = Poly(u0), Poly(v0), Poly(w0)

    x = (1/phi) * variable(p) # reverse preconditioning
    u, v, w = map(monic, (polyval(u, x), polyval(v, x), polyval(w, x)))

    u,v,w,err
end






end
