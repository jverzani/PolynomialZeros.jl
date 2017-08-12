module AGCD
using Polynomials

using Compat
import Compat.iszero

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
_polyval(p, x) = "XXX"
_monic(p::Vector) = p./p[end]
_monic!(p::Vector) = p ./= p[end]
_float(p::Vector) = float(p)

function _polymul{T,S}(p::Vector{T}, q::Vector{S})
    R = promote_type(T,S)
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
function _polyder{T}(p::Vector{T})
    S = eltype(float(p[1]))
    n = length(p)
    a2 = Vector{S}(n-1)
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
function cauchy_matrix!{T}(M, p::Vector{T}, k::Integer)
    n = length(p)
    for i in 1:k
        M[(1:n) + (i-1), i] = p
    end
end
function cauchy_matrix{T}(p::Vector{T}, k::Integer)
    n,m = cauchy_matrix_size(p, k)
    out = zeros(T, n, m)
    cauchy_matrix!(out, p, k)
    out
end

function sylvester_matrix_size(p, q, k::Int=0)
    n, m = length(p), length(q)
    
    if n < m
        n,m = m,n
    end
    Δ = n - m


    i,j = Δ + k, k
    a,b = cauchy_matrix_size(q,i)
    c,d = cauchy_matrix_size(p, j)
    (a, b + d)
end

function sylvester_matrix!(M, p::Vector, q::Vector, k::Int=0)
    @assert k >= 0
    n,m = length(p)-1, length(q) - 1
    if n < m
        p,q = q,p
        n,m = m,n
    end

    # M is u × w = (v1+v2) × w
    Δ = n - m
    i,j = k + Δ, k
    
    v,w = cauchy_matrix_size(p, j)
    u = size(M)[2]

    
    cauchy_matrix!(view(M, :, 1:(u-w)),   q, i)
    cauchy_matrix!(view(M, :, (u-w+1):u), p, j)
    
end

function sylvester_matrix{T,S}(p::Vector{T}, q::Vector{S}, k::Int=0)
    R = promote_type(T,S)
    M = zeros(R, sylvester_matrix_size(p,q,k)...)
    sylvester_matrix!(M, p, q, k)
    M
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
    
function JF!{T}(M, u::Vector{T}, v, w)
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
    cauchy_matrix!(view(M, 1 + (1:ai), 1:aj),        v, m + 1)
    cauchy_matrix!(view(M, 1 + (1:bi), aj + (1:bj)), u, k + 1)
#    M[1 + (1:ci), aj + bj + (1:cj)] = zeros(T, n+1, m+1-j)
    cauchy_matrix!(view(M, 1 + ai + (1:di), 1:dj),   w, m + 1)
#    M[1 + bi + (1:ei), dj + (1:ej)] = zeros(T, m + 1, n + 1 - j)
    cauchy_matrix!(view(M, 1 + ci + (1:fi), dj + ej + (1:fj)), u, j + 1)

    
end
function JF{U,V,W}(u::Vector{U}, v::Vector{V}, w::Vector{W})
    R = promote_type(U,V, W)
    M = zeros(R, JF_size(u, v, w)...)
    JF!(M, u, v, w)
    M
end

## ----------------------------------------

## converge on right singular vector of A associated with singular value sigma (not a double value)
## we return if sigma < tol; delta \appro 0;
## how to @assert that sigma_2 > sigma_1?
## claim is that this could be improved were A recycled
function lemma24{T}(p::Vector{T}, q::Vector{T}, k::Int, θ=1e-8) 

    A = sylvester_matrix(p, q , k)
    const MAXSTEPS = 100
    Q,R = Base.qr(A)

    ## if R is rank deficient, we can have issues solving R' \ x below
    ## What to do in this case with Big values???

    if iszero(det(R))
        u = q
        w = [one(T)]
        v =  cauchy_matrix(u, length(p) - length(u) + 1) \ p; _monic(v)
     
        return (zero(T), vcat(v,w))
    end
   
    #    x = rand(size(A)[2]);
    #    x ./= norm(x,2)   # use random initial guess
    x = ones(T, size(A)[2]) # use fixed initial guess
    σ, σ1 = 1e8, Inf

    ## how long do we update? Size is only issue, not accuracy so we iterate until we stop changing much
    m, n = size(R)
    y = zeros(T, m)
    z = zeros(T, n)
    flag = :maxteps
    for cnt in 1:MAXSTEPS
        y[:] = R' \ x
        z[:] = R \ y
        x[:] = z/norm(z,2)
        sigma = norm(R * x, 2)  # y/norm(z,2)
        σ1 = abs(sigma)

        if σ1 < θ
            flag = :threshhold
            break
        end
        if  (abs((σ - σ1) / σ1) < 1.1)
            flag = :pause
            break
        end
        
        σ = σ1
    end

    return (σ1, x)
end



# eqn (12) to get weights
function _weights{T}(p::Vector{T},q)
    n, m = length(p), length(q)
    wts = ones(T, 1 + n + m)
    for (i,pj) in enumerate(p)
        wts[1+i] = 1/max(1, abs(pj))
    end
    for (i,qj) in enumerate(q)
        wts[1 + n+i] = 1/max(1, abs(qj))
    end    
    wts
end

# reduce residual by Gauss-Newton algorithm
# updates u,v,w; returns residual error and flag for convergence
function reduce_residual!{T}(u,v,w, p::Vector{T}, q, wts, ρ)
    m, n, l = map(_degree, (u, v, w))
    # preallocate
    A = zeros(T, JF_size(u, v, w)...)
    b = zeros(T, 1 + length(p) + length(q))
    inc = zeros(T, m + n + l + 3) #weighted_least_square(A, b, wts)

    ρm_1, ρm = Inf, Inf
    const MAXSTEPS = 20
    ctr = 0
    flag = :not_converged

    
    
    while ctr <= MAXSTEPS
        ρm =  agcd_update!(p, q, A, b, inc, u, v, w, m, n, wts)

        np2 = norm(p,2)
        nu2 = norm(u,2)
        
        if ρm > 1.5 * ρm_1
            return (ρm, flag)
        end
        # threshold here is really important, but not clear what it should be
        # we empirically get better results---it seems---with sqrt(norm(p,2))
        # though recomendation if norm(u,2); Might also make sense to include
        # q in the computation
        if ρm <= ρ * norm(u,2) #sqrt(norm(p,2)) 
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
    W = diagm(w)
    (W * A) \ (W * b)
end
    
function weighted_least_square!(M, A, b, w)
    W = diagm(w)
    M[:] = (W * A) \ (W * b)
end



## compute F(u,v,w) - [p, p'] = [u*v, u*w] - [p, p']
Fmp(p,q,u,v,w) =  [u[end]; _polymul(u,v); _polymul(u,w)] .- [1; p; q]
function Fmp!(b, p,q,u,v,w)
    b[1] = u[end] - 1
    for (i,v) in enumerate(_polymul(u,v) .- p)
        b[i] = v
    end
    offset = 1 + length(p)
    for (i,v) in enumerate(_polymul(u,w) .- q)
        b[offset + i] = v
    end
end

residual_error(p,q,u,v,w, wts=ones(1 + length(p) + length(q))) = norm( ([_polymul(u,v); _polymul(u,w)] .- [p; q]) .* wts[2:end], 2)

### refine u,v,w estimate using
## Newton's method and weighted least squares
# updates in place u,v,w
function agcd_update!(p, q, A, b, inc, u, v, w, m, n, wts)
    JF!(A, u, v, w)
    Fmp!(b, p,q,u,v,w)
    weighted_least_square!(inc, A, b, wts)

    Δu = inc[1:(1+m)]
    Δv = inc[(m+2):(m+n+2)]
    Δw = inc[(m+n+3):end]

    u .-= Δu; _monic!(u)
    v .-= Δv; _monic!(v)
    w .-= Δw; _monic!(w)

    
    # return error estimate
    err = residual_error(p,q,u,v,w, wts)
#    println("err=$err")
    err
end


## ----------------------

## some diagnostic functions
function sylvester_matrix_singular_values{T}(p0::Vector{T}, q0::Vector=_polyder(p0))
    p, q = _float(p0), _float(q0)
    p, q = _monic(p), _monic(q)
    n, m = _degree(p), _degree(q)

    n < m && return sylvester_matrix_ranks(q, p)

    
    ψs = zeros(T, m-1)
    for k in 1:(m-1)
        ψ, x = lemma24(p, q, k)
        ψs[k] = ψ
    end
    ψs
end

# find u,v,w,err assuming gcd is rank k
function rank_k_agcd(k, p::Vector, q::Vector=_polyder(p), θ=1e-8, ρ=1e-10)
    ψ, x = lemma24(_monic(p), _monic(q), k, θ)
    v = x[1:(length(x)-k)]; _monic(v)
    w = x[(length(x)-k+1):end]; _monic(w)
    # solve for u using least-squares, not division
    u = cauchy_matrix(v, length(p) - length(v) + 1) \ p; _monic(u)

    ρm, flag = reduce_residual!(u,v,w, p, q, _weights(p, q), ρ)

    (u, v, w, ρm)
end
    

## ----------------------


"""
Return k, u,v,w where k reveals rank; u*v ≈ p; u*w ≈ q; v & w coprime

Following Zeng, section 4, this could be made more efficient by recycling QR decomposition

Default parameter choice comes from Zeng

Can we engineer around the minimum?
"""
function reveal_rank{T}(p::Vector{T}, q::Vector{T}, θ=1e-8, ρ=1e-10)
    # we check k until 1) lemma24 is small *and* we can refine residual error
    n,m = length(p), length(q)  # assume n >= m

    θ = θ * norm(p,2)  # p23
    wts = _weights(p,q)

    # umin = T[]; vmin=T[]; wmin=T[];
    # ρm_min1 = Inf
    # ρm_min = Inf
    # mk = 0

    # psis = zeros(T, m-1)
    
    for k in 1:m-1
        psi, x = lemma24(p,q,k, θ)
        # psis[k] = psi
        
        if abs(psi) < norm(p,2) * θ #  norm(p,2)? norm(A,2)?

            
            ## degree u = n - k; degree(v) = k
            v = x[1:(length(x)-k)]; _monic(v)
            w = x[(length(x)-k+1):end]; _monic(w)
            # solve for u using least-squares, not division
            u = cauchy_matrix(v, length(p) - length(v) + 1) \ p; _monic(u)

            ρm, flag = reduce_residual!(u,v,w, p, q, wts, ρ)
            if flag == :converged
                return (k, u, v, w)
            # elseif ρm < ρm_min 
            #    umin = u; vmin=v; wmin = w; ρm_min = ρm; ρm_min1 = ρm_min; mk=k
            end
        end
    end

    ## This is tricky to get correct! The psis should be 0, but what is 0 is unclear.
    ## Following Zeng, we have norm(p,2)*θ *and* the residual error less than sqrt(norm(p,2))*ρ
    ## Zeng uses norm(U,2) there.
    ## is by no means perfect.
    
    # if ρm_min < 2/3 * ρm_min1
    #     println("Is k=$mk a winner? ρm = $ρm_min next is $ρm_min1")
    # end
    
    return (n, ones(T,1), p, q)


end 


"""

    `agcd(p, q, θ=1e-8, ρ=1e-10)`

Find an approximate GCD for polynomials `p` and `q` using an algorithm of [Zeng](http://www.ams.org/journals/mcom/2005-74-250/S0025-5718-04-01692-8/home.html). 


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
sylvester matrix is rank deficient. This is done in `lemma24` which
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
function agcd{T <: Number,S <: Number}(p::Vector{T}, q::Vector{S}=_polyder(p);
                                       θ = 1e-8,  
                                       ρ = 1e-10   
              )

    n, m = length(p), length(q)
    if m > n
        p, q=q, p
    end

    R = promote_type(T,S)
    p0 = float.(convert(Vector{R}, p)); _monic!(p0)
    q0 = float.(convert(Vector{R}, q)); _monic!(q0)
    

    if m == 0
        return (ones(R,1), p0, q0, zero(R))
    end

    k, u, v, w = reveal_rank(p0, q0, θ, ρ)
    
    (u, v, w, residual_error(p0, q0, u, v, w))
end


# for polys, we precondition

monic(p) = p/p[end]
_float{T <: AbstractFloat}(p::Poly{T}) = p
_float{T}(p::Poly{T}) = Poly(float.(p.a), p.var)

## --------------
## preconditioning code
## taken from https://who.rocq.inria.fr/Jan.Elias/pdf/je_masterthesis.pdf
function geometric_mean{T}(a::Vector{T})::T
    b = filter(!iszero, a)
    n = length(b)
    prod(abs(u)^(one(T)/n) for u in b)  
end

function ratio(p,q, atol=Base.eps(), rtol=Base.eps())
    as = norm.(filter(!iszero, p.a))
    length(as) == 0 && return Inf

    bs = norm.(filter(!iszero, q.a))
    length(bs) == 0 && return Inf

    max(maximum(as), maximum(bs)) / min(minimum(as), minimum(bs))
end

## find alpha, gamma that minimize ratio of max coefficients to min coefficients, for getting zeros
## 1.12 writes this as a linear programming problem, we just ...
function precondition{T}(p::Poly{T}, q::Poly)

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

## -----------------

function agcd{T <: Number,S <: Number}(p0::Poly{T}, q0::Poly{S}=polyder(p0);
                                   θ = 1e-8,  
                                   ρ = 1e-10)

    p, q = _float(p0), _float(q0)
    
    p, q, phi, alpha = precondition(p, q) 
    
    u0,v0,w0,err = agcd(p.a, q.a; θ=θ, ρ=ρ)
    u,v,w = Poly(u0), Poly(v0), Poly(w0)

    x = (1/phi) * variable(p) # reverse preconditioning
    u, v, w = map(monic, (polyval(u, x), polyval(v, x), polyval(w, x)))

    u,v,w,err
end
              

    



end
