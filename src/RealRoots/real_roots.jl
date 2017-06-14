module RealRoots

## This is pretty slow:
## * using `Poly` class is slow as new polys allocate
## * the poly translation: Tλ, Hλ, and R are slow. Could rewrite using arrays too
## This is not robust
## * the newton search is more primitive than suggested, as it cuts down on calls to admissible points
## * we don'tuse interval arithmetic, as it take about 3 times as long per computaion. This would be worth it, but...
## * we aren't careful with the floating point size.


using Polynomials
using Compat
using ..AGCD
import Roots # fzero

## Find real roots of a polynomial using a DesCartes Method
##
## State of the art for this approach is summarized here:
## Cf. "Computing Real Roots of Real Polynomials ... and now For Real!" by Alexander Kobel, Fabrice Rouillier, Michael Rouillier
## (https://arxiv.org/abs/1605.00410) 
## 
## Earlier work of a Descartes-like method  [Vincent's theorem](http://en.wikipedia.org/wiki/Vincent's_theorem)
## are here "Efficient isolation of polynomial’s real roots" by
## Fabrice Rouillier; Paul Zimmermann
##
## This implementation doesn't take nearly enough care with the details, but takes some ideas
## to implement a means to find real roots of non-pathological polynomials (lowish degree, roots separated)
##
## implementation using interval arithmetic seemed to be slower, but no more succesful in avoiding numeric issues

## Polynomial transformations
##
## The Taylor shift here is the most expensive operation
## https://arxiv.org/pdf/1605.00410.pdf has a better strategy
## of partial Taylor shifts, using just the nearby roots
## 
## We did not implement the full Newton test, but this would need to
## be implemented too.

" `p(x + λ)`: Translate polynomial left by λ "
Tλ(p, λ=1)   = polyval(p, variable(p) + λ)

" `R(p)` finds  `x^n p(1/x)` which is a reversal of coefficients "
R(p) = Poly(reverse(p.a), p.var)

" `p(λ x)`: scale x axis by λ "
Hλ(p, λ=1//2) = polyval(p, λ * variable(p))




## Upper bound on size of real roots that is tighter than cauchy
## titan.princeton.edu/papers/claire/hertz-etal-99.ps
function upperbound{T}(p::Poly{T})::T
    descartes_bound(p) == 0 && return zero(T)
    
    q, d = p/p[end], degree(p)
    
    d == 0 && error("degree 0 is a constant")
    d == 1 && abs(q[0])


    a1 = abs(q[d-1])
    B = maximum([abs(q[i]) for i in 0:(d-2)])

    a,b,c = 1, -(1+a1), a1-B
    (-b + sqrt(b^2 - 4a*c))/2
end

function lowerbound{T <: Real}(p::Poly{T})::T
    q = p(-variable(p))
    -upperbound(q)
end


"""

Use Descarte's rule of signs to compute a bound on the number of *postive* roots of p: n-2k

"""
## count upper bound on number of positive roots
## use -1 to signal an issue
function descartes_bound{T}(p::Poly{T})
    sgs = Int[]
    n = length(p.a)
    for i in n:-1:1
        s = sign(p.a[i])
        if s < 0
            push!(sgs, -1)
        elseif s > 0
            push!(sgs, 1)
        elseif iszero(s)
            # do not push zero
        end
    end
    count_sign_changes(sgs)
end

# function descartes_bound{T}(p::Poly{ValidatedNumerics.Interval{T}})
#     sgs = Int[]
#     n = length(p.a)
#     for i in n:-1:1
#         s = sign(p.a[i])
#         if s < 0
#             push!(sgs, -1)
#         elseif s > 0
#             push!(sgs, 1)
#         elseif iszero(s.lo) && iszero(s.hi)
#             # do not push zero
#         else
#             return -1 # inconclusive
#         end
#     end
#     count_sign_changes(sgs)
# end



function count_sign_changes(as::Vector)
    length(as) == 0 && return -1
    cnt, sgn = 0, sign(as[1])
    for i in 2:length(as)
        nsgn = sign(as[i])
        if nsgn * sgn < 0
            sgn = nsgn
            cnt += 1
        end
    end
    cnt
end


## Descarte's upperbound on possible zeros in (0,1)
function DesBound(p::Poly)
    q = Tλ(R(p),1)
    descartes_bound(q)
end



## Interval with [a,b], N
## N is for Newton-Test
#mutable struct Intervalab{T,S}
type Intervalab{T,S}
    a::T
    b::S
    N::Int
    bnd::Int
end
width(I::Intervalab) = I.b-I.a
midpoint(I::Intervalab) = I.a + (0.5) * width(I)


## how we store the state of our algorithm to find zero in [0,1]
#@compat struct State{T}
immutable State{T}
    Internal::Vector{Intervalab{T,T}}                    # DesBound > 1
    Isol::Vector{Intervalab{T,T}}                        # DesBound == 1
    Unresolved::Vector{Intervalab{T,T}}
    p::Poly{T}
end

State{T <: AbstractFloat}(p::Poly{T}) = State(Intervalab{T,T}[], Intervalab{T,T}[], Intervalab{T,T}[], p)


DesBound(p, a::Real, b::Real) = DesBound(Pab(p, a, b))
DesBound(p, node::Intervalab) = DesBound(p, node.a, node.b)

function DescartesBound_ab{T}(st::State{T}, node)

    a, b = node.a, node.b
    
    p1 = Hλ(Tλ(st.p, a), b-a)
    p2 = Tλ(R(p1),1)
    cnt = descartes_bound(p2)

    return cnt
end

## Tests

## Return true or false
zero_test(st::State, node)  =   DescartesBound_ab(st, node) == 0 
one_test(st::State, node)  =   DescartesBound_ab(st, node) == 1  

## return count -1 (can't tell), 0, 1, or more
zero_one_test(st::State, node) = DescartesBound_ab(st, node)  


# find admissible point
# XXX improve me
function find_admissible_point{T}(st::State{T},  I::Intervalab, m=midpoint(I), Ni::T=one(T), c::T=one(T))::T
    N = ceil(Int, c * degree(st.p)/2)
    ep = min(m-I.a, I.b - m) / (4*Ni)
    mis = [m + i/N * ep for i in -N:N]
    for m in shuffle(mis)
        (m < I.b || m > I.a) || continue
        norm(st.p(m)) > 0 && return m
    end
#    mx, i = findmax(norm.(st.p.(mis)))
#    mis[i]
end


## split an interval up into 2^n pieces. This is useful at the outset
function break_up_interval{T,R}(st::State{T}, I::Intervalab{R}, n::Int=5)
    ointervals = [I]
    N = I.N
    while n > 0
        intervals = Any[]
        for node in ointervals
            mi = find_admissible_point(st, node)
            u,v=promote(node.a, mi)
            push!(intervals, Intervalab(u,v, N, -2))
            v,w = promote(mi, node.b)
            push!(intervals, Intervalab(v, w, N, -2))
        end
        ointervals = copy(intervals)
        n -= 1
    end
    ointervals
end

# find splitting point
# find an admissible point that can be used to split interval. Needs to have all
# coefficients not straddle 0
# return (logical, two intervals)
function split_interval{T}(st::State{T},I::Intervalab,  m=midpoint(I), Ni=one(T), c=one(T))
    N = ceil(Int, c * degree(st.p)/2)
    ep = min(1, width(I)) / (16*Ni)
    mis = T[m + i/N * ep for i in -N:N]
    mis = filter(m -> m > I.a && m < I.b, mis)
    mx, i = findmax(norm.(st.p.(mis)))

    ## Now, we need apoint that is bigger than max and leaves conclusive
    for i in eachindex(mis)
        mi = mis[i]
        abs(st.p(mi)) >= min(mx/4, one(T)) || continue
        ileft = Intervalab(I.a, mi, I.N, -2)
        nl = DescartesBound_ab(st, ileft)        
        nl == -1 && continue
        iright = Intervalab(mi, I.b, I.N, -2)
        nr = DescartesBound_ab(st, iright)        
        nr == -1 && continue

        # XXX improve this XXX
        # identify degenerate cases here. This is not good.
        # if nl + nr < I.bnd
        #     if nl == 0 || nr == 0
        #         if (nl == 0 && rem(I.bnd - nr,2) == 1) || ( nl == 0 && rem(I.bnd - nr,2) == 1)
        #             println("nl=$nl, nr=$nr, bnd=$(I.bnd) -- is this an error")
        #             return(false, I, I)
        #         end
        #     elseif nl + nr < I.bnd
        #         println("nl=$nl, nr=$nr, bnd=$(I.bnd) -- is this an error")
        #         return(false, I, I)
        #     end
        # end
        
        ileft.bnd = nl; iright.bnd=nr

##        println("Split $(I.a) < $mi < $(I.b): $(nl) & $(nr)")
        
        
        return (true, ileft, iright)
    end

    return (false, I, I)
end


## split interval
## adds to intervals if successful
## adds node to Unresolved if not
function linear_step{T}(st::State{T}, node)
    succ, I1, I2 = split_interval(st, node)

    if succ
        push!(st.Internal, I1)
        push!(st.Internal, I2)
    else
        push!(st.Unresolved, node)
    end
    return true
    
end


## return (true, I), (false, node)
## I will be smaller interval containing all roots in node
function newton_test{T}(st::State{T}, node)
    const NMAX = 1024  # 2147483648 = 2^31
    (node.N > NMAX) && return (false, node) 
    (zero_one_test(st, node) in (0,1)) && return(false, node)

    a, b, w, N, bnd  = node.a, node.b, width(node), node.N, node.bnd

    a1 = a - st.p(a) / polyder(st.p)(a)
    b1 = b - st.p(b) / polyder(st.p)(b)

    if a < a1 && zero_test(st, Intervalab(a, a1, N, 0))
        if b1 < b && zero_test(st, Intervalab(b1, b, N, 0))
            return (true, Intervalab(a1, b1, N*N, -2))
        else
            return(true, Intervalab(b1, b, N*N, -2))
        end
    elseif b1 < b && zero_test(st, Intervalab(b1, b, N, 0))
        return (true, Intervalab(a, b, N*N, -2))
    end

    ## boundary test?

    mlstar::T = find_admissible_point(st, node, a + w/(2N))
    if mlstar > a && zero_test(st, Intervalab(mlstar, b, N,-2))
        return (true, Intervalab(a, mlstar, N,-2))
    end

    mrstar::T = find_admissible_point(st, node, b - w/(2N))
    if mrstar < b && zero_test(st, Intervalab(a, mrstar, N,-2))
        return (true, Intervalab(mrstar, b, N,-2))
    end
    return (false, node)
end

## This is the newton test from the paper. The above one avoids doing more checks which are costly
function newton_test_XX{T}(st::State{T}, node)

    const NMAX = 2147483648  # 2^31
    (node.N > NMAX) && return (false, node)  
    (zero_one_test in (0,1)) && return(false, node)

    a, b, w, N, bnd  = node.a, node.b, width(node), node.N, node.bnd
    psis = T[a + (j/4) * w for j in 1:3]
    psi_stars = T[find_admissible_point(st, node, psi) for psi in psis]
    vs = st.p.(psi_stars) ./ polyder(st.p).(psi_stars)

    lambdas = zeros(T, 3, 3)
    for i in 1:3
        for j in (i+1):3
            
            lambdas[i,j] = psi_stars[i] - (psi_stars[j] - psi_stars[i]) / (vs[j] -vs[i]) * vs[i]
            # compute approximations
            if  lambdas[i,j]  < a ||  lambdas[i,j]  > b
                ## println("outside [a.b]")
                continue
            end
            lij = floor((lambdas[i,j] - a) * 4* N / w)

            aij = a + max(zero(T), lij - 1) * w / (4N)
            ## println("aij=$aij, admiss=$(find_admissible_point(st, node, aij))")
            isnan(aij) && continue
            aij = find_admissible_point(st, Intervalab(a, b, node.N, -2), aij)

            bij = a + min(4N, lij + 2) * w / (4N)
            isnan(bij) && continue            
            bij = find_admissible_point(st, Intervalab(a, b, node.N, -2), bij)

            println("Okay, zero test?")
            println(zero_test(st, Intervalab(a, aij, N,0)))
            println(zero_test(st, Intervalab(bij, b, N,0)))
            println("$a < $aij < $bij < $b")
            println("===")
            ## okay, how did we do?
            if (aij >= a && bij <= b) && zero_test(st, Intervalab(a, aij, N,0)) && zero_test(st, Intervalab(bij, b, N,0)) 
#                println("squeeze: a=$aij, b = $bij")
                return (true, Intervalab(aij, bij, N^2, node.bnd))
            end
        end
    end
    # boundary test
    mlstar::T = find_admissible_point(st, node, a + w/(2N))
    mrstar::T = find_admissible_point(st, node, b - w/(2N))
    al = Intervalab(a, mlstar, N,-2)
    br = Intervalab(mrstar, b, N,-2)

    if mlstar > a && zero_test(st, Intervalab(mlstar, b, N,-2))
        return (true, Intervalab(a, mlstar, N,-2))
    elseif mrstar < b && zero_test(st, Intervalab(a, mrstar, N,-2))
        return (true, Intervalab(mrstar, b, N,-2))
    else
        return (false, node)
    end
end



## Add successors to I
## We have
function addSucc{T}(st::State{T}, node)

    val, I = newton_test(st, node)

    if val
        ## a,b = node.a, node.b
        ## println("Newton test: a=$a <= $(I.a) <= $(I.b) <= $b=b")
        push!(st.Internal, I)
    else
        ##        println("linear step")
        succ = linear_step(st, node)
        if !succ
            warn("node $node was a failure")
        end
    end
    true
end


## m, M should bound the roots
## essentially algorithm 4
function ANewDsc{T <: Real}(p::Poly{T}, m = lowerbound(p), M=upperbound(p))

    st = State(p)
    base_node = Intervalab(m, M, 4, -2)    
    base_node.bnd = DescartesBound_ab(st, base_node)


    if base_node.bnd == -1
        append!(st.Internal, break_up_interval(st, base_node, 4))
    else
        push!(st.Internal, base_node)
    end
    
    while length(st.Internal) > 0
        node = pop!(st.Internal)

        bnd = node.bnd
        if bnd == -2
            bnd = DescartesBound_ab(st, node) 
            node.bnd = bnd
        end

        if bnd < 0
            # this is a bad node!
            warn("Bad node, how did it get here: $node")
        elseif bnd == 0
            next
        elseif bnd == 1
            push!(st.Isol, node)
            next
        else
            addSucc(st, node)
        end
    end
    st
end

# populate `Isol`
# p must not have any roots with even degree. (e.g. no (x-c)^(2k) exists as a factor for any c,k
# assumed square free (at least no roots of even multiplicity)
function isolate_roots{T <: Real}(p::Poly{T}, m, M)

    try
        st = ANewDsc(p, m, M)
        return st
    catch err
        if  !(T <: BigFloat)
            try
                st = ANewDsc(convert(Poly{BigFloat}, p), m, M)
                return st
            catch err
                rethrow(err)
            end
        end
    end
        
end



function real_roots{T <: Real}(p::Poly{T}, m = lowerbound(p), M=upperbound(p); square_free::Bool=true)

    # deflate zero
    nzroots = findfirst(p.a) - 1
    if nzroots > 0
        p = Poly(p.a[nzroots+1:end])
    end
    
    
    
    if !square_free
        u,v,w,err = AGCD.agcd(p, polyder(p))
    else
        v = p
    end

    st = isolate_roots(v, m, M)

    if length(st.Unresolved) > 0
        println("Some intervals are unresolved:")
        println("------------------------------")
        for node in st.Unresolved
            @printf "* There may be up to %d roots in (%0.16f, %0.16f).\n" node.bnd node.a node.b
        end
        println("------------------------------")                
    end

    rts = zeros(T, length(st.Isol))
    for i in eachindex(st.Isol)
        node = st.Isol[i]
        ## println("find root in $(node.a), $(node.b), $(p(node.a)), $(p(node.b))")
        rt = try
            Roots.fzero(x -> v(x), node.a, node.b)
        catch err
            Roots.fzero(x -> v(x), big(node.a), big(node.b))
        end
        rts[i] = rt
    end

    if nzroots > 0
        rts = push!(rts, zero(T))
    end
    
    rts
end
        



end



