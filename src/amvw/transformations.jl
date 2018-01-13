##
##################################################

"""
 rotations; find values
 Real Givens
 This subroutine computes c and s such that,

 [c -s] * [a, b] = [r,0]; c^2 + s^2 = 1

 and 

 r = sqrt(|a|^2 + |b|^2).

  XXX seems faster to just return r, then not
"""
function givensrot(a::T,b::T) where {T <: Real}
    iszero(b) && return (sign(a) * one(T), zero(T), abs(a))
    iszero(a) && return(zero(T), -sign(b) * one(T), abs(b))

    r = hypot(a,b)
    return(a/r,-b/r,r)
end

## givens rotation
##################################################
# Compute Givens rotation zeroing b
#
# G1 [ ar + i*ai ] = [ nrm ]
# G1 [    b      ] = [     ]
#
# all variables real (nrm complex)
# returns (copmlex(cr, ci), s) with
# u=complex(cr,ci), v=s; then [u -v; v conj(u)] * [complex(ar, ai), s] has 0
#
# XXX: Could hand write this if needed, here we just use `givens` with a flip
# to get what we want, a real s part, not s part
function givensrot(a::Complex{T},b::Complex{T}) where {T <: Real}
    G, r = givens(b, a, 1, 2)
    G.s, -G.c, r
end
givensrot(a::Complex{T},b::T) where {T <: Real} = givensrot(a, complex(b, zero(T)))


####   Operations on [,[ terms

## The zero_index and stop_index+1 point at "D" matrices
## 
## Let a D matrix be one of [1 0; 0 1] or [-1 0; 0 1] (D^2 = I). Then we have this move
## D    --->   D  (we update the rotator)
##   [       [
##
## this is `dflip`
function dflip(a::RealRotator{T}, d=one(T)) where {T}
    a.s = sign(d)*a.s
end

# get d from rotator which is RR(1,0) or RR(-1, 0)
function getd(a::RealRotator{T}) where {T}
    c, s = vals(a)
    norm(s) <= 4eps(T) || error("a is not a diagonal rotator")
    sign(c)
end

## This is main case
#  Q           D Q
#     D --> D 
    

"""
   D  --> D
U           V
"""
function Dflip(r::ComplexComplexRotator{T}, d::ComplexComplexRotator{T}) where {T}
    !is_diagonal(d) && error("d must be diagonal rotator")

    # D is fixed,
    alpha = d.c
    r.s = r.s * conj(alpha)
end

##   U --> U Da
## D           Da
## (not the reverse!)
function Dflip(d::ComplexRealRotator{T}, r::ComplexRealRotator{T}) where {T}
#    !is_diagonal(d) && error("d must be diagonal rotator")

    alpha = d.c
    c,s = vals(r)
    vals!(r, c*conj(alpha), s)
end

## We have this for left fuse and for deflation
#
#  Di                         Di                     
#    Qi+1             Si+1     Di+1       Si+1   
#        Qi+2    -->    Si+2    Di+2    =    Si+2  * diagm([alpha, I, conj(alpha)])
#           ...           ...     ...          ...
#             Qj            Sj      Dj          Sj
function cascade(Qs, D, alpha, i, j)
    # Q = CR(c,s) -> S = CR(c*conj(alpha), s)
    for k in (i+1):j
        c,s = vals(Qs[k])
        vals!(Qs[k], c*conj(alpha), s)
    end
    D[i] *= alpha
    D[j+1] *= conj(alpha)
end




## Fuse
## fuse combines two rotations, a and b, into one,


## For ComplexRealRotator, the result of a*b will not have a real sign
## we output by rotating by alpha.
## return alpha so a*b = f(ab) * [alpha 0; 0 conj(alpha)]
## for left with have uv -> (u') Di
function fuse(a::ComplexRealRotator{T}, b::ComplexRealRotator{T},::Type{Val{:left}}) where {T}
    #    idx(a) == idx(b) || error("can't fuse")
    u = a.c * b.c - conj(a.s) * b.s
    v = conj(a.c) * b.s + a.s * b.c

    alpha =  conj(v)/norm(v)

    a.c = u * alpha
    a.s = norm(v)

    conj(alpha)
end

# for right we have uv -> (v') Di
function fuse(a::ComplexRealRotator{T}, b::ComplexRealRotator{T}, ::Type{Val{:right}}) where {T}
#    idx(a) == idx(b) || error("can't fuse")
    u = a.c * b.c - conj(a.s) * b.s
    v = conj(a.c) * b.s + a.s * b.c

    
    alpha =  conj(v)/norm(v)

    b.c = u * alpha
    b.s = norm(v)

    conj(alpha)
end


## Fuse for general rotation
## We have two functions as it seems a bit faster
fuse(a::Rotator{T}, b::Rotator{T}, dir, d) where {T} = fuse(a,b,dir)

function fuse(a::Rotator{T}, b::Rotator{T},::Type{Val{:left}}) where {T}
    #    idx(a) == idx(b) || error("can't fuse")
    u = a.c * b.c - conj(a.s) * b.s
    a.s = conj(a.c) * b.s + a.s * b.c
    a.c = u

    one(T)
end
function fuse(a::Rotator{T}, b::Rotator{T}, ::Type{Val{:right}}) where {T}
#    idx(a) == idx(b) || error("can't fuse")
    u = a.c * b.c - conj(a.s) * b.s
    b.s = conj(a.c) * b.s + a.s * b.c
    b.c = u

    one(T)
end


# Turnover: Q1    Q3   | x x x |      Q1
#              Q2    = | x x x | = Q3    Q2  <-- misfit=3 Q1, Q2 shift; 
#                      | x x x |
#
# misfit is Val{:right} for <-- (right to left turnover), Val{:left} for -->
#
# This is the key computation once matrices are written as rotators
# We wrote this for complex rotators where sine part may be complex
# so we make use of alpha and beta, which isn't otherwise needed
# could streamline, but doesn't seem to incur an expense

function _turnover(Q1::Rotator{T}, Q2::Rotator{T}, Q3::Rotator{T}) where {T}    
#    i,j,k = idx(Q1), idx(Q2), idx(Q3)
#    (i == k) || error("Need to have a turnover up down up or down up down: have i=$i, j=$j, k=$k")
#    abs(j-i) == 1 || error("Need to have |i-j| == 1")
    
    c1,s1 = vals(Q1)
    c2,s2 = vals(Q2)
    c3,s3 = vals(Q3)

    # key is to find U1,U2,U3 with U2'*U1'*U3' * (Q1*Q2*Q3) = I
    # do so by three Givens rotations to make (Q1*Q2*Q3) upper triangular
    
    # initialize c4 and s4
    a = conj(c1)*c2*s3 + s1*c3 
    b = s2*s3
    # check norm([a,b]) \approx 1    
    c4, s4, temp = givensrot(a,b)#, Val{true})

    # initialize c5 and s5

    a = c1*c3 - conj(s1)*c2*s3
    b = temp
    # check norm([a,b]) \approx 1
    c5, s5, alpha = givensrot(a, b)

    alpha = alpha/norm(alpha)
    c5 *= conj(alpha) # make diagonal elements 1
    s5 *= alpha
    
    # second column
    u = -c1*conj(s3) - conj(s1)*c2*conj(c3)
    v = conj(c1)*c2*conj(c3) - s1*conj(s3)
    w = s2 * conj(c3)

    a = c4*conj(c5)*v - conj(s4)*conj(c5)*w + s5*u
    b = conj(c4)*w + s4*v

    c6, s6, beta = givensrot(a,b)

    beta = beta/norm(beta)
    c6 *= conj(beta) # make diagonal elements 1
    s6 *= beta
    
    (c4, s4, c5, s5, c6, s6)
end




function turnover(Q1::Rotator{T}, Q2::Rotator{T}, Q3::Rotator{T}, ::Type{Val{:right}}) where {T}

    c4,s4,c5,s5,c6,s6 = _turnover(Q1,Q2,Q3)
    vals!(Q3, conj(c4), -s4)
    vals!(Q1, conj(c5), -s5)
    vals!(Q2, conj(c6), -s6)
    idx!(Q3, idx(Q2))  # misfit is right one
    
end

turnover(Q1::Rotator{T}, Q2::Rotator{T}, Q3::Rotator{T}) where {T} = turnover(Q1, Q2, Q3, Val{:right})

function turnover(Q1::Rotator{T}, Q2::Rotator{T}, Q3::Rotator{T}, ::Type{Val{:left}}) where {T}
    
    c4,s4,c5,s5,c6,s6 = _turnover(Q1,Q2,Q3)
    vals!(Q2, conj(c4), -s4)
    vals!(Q3, conj(c5), -s5)
    vals!(Q1, conj(c6), -s6)
    idx!(Q1, idx(Q2))   # misfit is left one

end





## passthrough
## Pass a rotator through a diagonal matrix with phase shifts
## D U -> U' D'
## Here D[i] = D[1]
## usually call with view(state.d, idx(U):idx(U)+1)
function passthrough(D, U::ComplexRealRotator{T}) where {T}
    i = idx(U)
    alpha, beta = D[i], D[i+1]

    c, s = vals(U)
    u = c * alpha * conj(beta)
    v = s
    vals!(U, u, v)
    D[i], D[i+1] = beta, alpha
end

function passthrough(D::ComplexRealRotator{T}, U::ComplexRealRotator{T}) where {T}
    norm(D.s) <= 1e2*eps(T) || error("D not diagonal")
    alpha, ds = vals(D)
    c,s = vals(U)
    vals!(U, c*alpha*alpha, s)
    vals!(D, conj(alpha), ds)
end
