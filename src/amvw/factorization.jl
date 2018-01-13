## Factorization code

## Factor
##
##  1 0 -p1                   1 0 0 0 0       -p1
##  0 1 -p2  --> Yn + X -->   0 1 0 0 0       -p2
##  0 0 -p3                   0 0 1 0 -1  +   -p3
##                            0 0 0 1 0       - 1
## we pass in ps = [-p1, -p2, ..., -pn, -1] the whol column including 1
## We can leave as Ct * B * D or D * Ct * B (the default)
##
function R_factorization(ps::Vector{Complex{T}}, Ct, B, side=Val{:left}) where {T}

    N = length(ps) - 1 # ps has 1 in it?
    c,s,tmp = givensrot(conj(ps[N]), -one(Complex{T}))

    nrm = norm(c)
    alpha = c/nrm
    
    alpha = alpha
    
    vals!(Ct[N], conj(c), -s);
    idx!(Ct[N], N)
    
    vals!(B[N], -s*alpha, norm(c))
    idx!(B[N], N)

    ## do we leave on left? If so, it must pass through Cts and B
    gamma = side == Val{:left} ? alpha : one(Complex{T})

    for i in (N-1):-1:1
        c,s,tmp = givensrot(ps[i], tmp)
        vals!(Ct[i], conj(c*gamma), -s)
        idx!(Ct[i], i)

        vals!(B[i], c*gamma, s)
        idx!(B[i], i)
    end

    side == Val{:left} ? alpha : conj(alpha)
end


## Factor R =
# 1 0 0 p1     1 0 0 0       p1
# 0 1 0 p2 --> 0 1 0 0  -->  p2
# 0 0 1 p3     0 0 0 -1      p3
#              0 0 1 0       -1
function R_factorization(ps::Vector{T}, Ct, B, side=Val{:not_applicable}) where {T <: Real}
    N = length(ps) - 1
    c,s,tmp = givensrot(ps[N], - one(T))

    vals!(Ct[N], c, -s);
    idx!(Ct[N], N)
    vals!(B[N], -s, c)
    idx!(B[N], N)
    
   
    for i in (N-1):-1:1
        c,s,tmp = givensrot(ps[i], tmp)
        vals!(Ct[i], c, -s)
        idx!(Ct[i], i)

        vals!(B[i], c, s)
        idx!(B[i], i)
    end
end


function Q_factorization(state::FactorizationType{T, St, P, Tw}) where {T, St, P, Tw}
    N = state.N
    Q = state.Q
    S = St == Val{:SingleShift} ? Complex{T} : T # type of cosine term in rotator
    for ii = 1:(N-1)
        vals!(Q[ii], zero(S), one(T))
        idx!(Q[ii], ii)
    end
    vals!(Q[N], one(S), zero(T)) #I
    idx!(Q[N], N)
end


## Init State uses factorization of R and Q above
function init_state(state::FactorizationType{T, St, P, Val{:NotTwisted}}, decompose) where {T, St, P}
    # (Q,D, sigma)
    # Ct, B [Ct1, B1]
    init_triu(state, decompose)
    Q_factorization(state)
end


function init_state(state::FactorizationType{T, St, P, Val{:IsTwisted}}, decompose) where {T, St, P}
    # (Q,D, sigma)
    # Ct, B [Ct1, B1]
    alpha = init_triu(state, decompose)
    Q_factorization(state)
    ## XXX ... Need to twist the Q here .... XXXX

end

## Upper triangular piecies R and (V,W)
## populate (D,Ct,B)
function init_triu(state::FactorizationType{T, Val{:SingleShift}, Val{:NoPencil}, Tw}, decompose) where {T, Tw}
    ps = decompose(state.POLY)
    alpha =  R_factorization(ps, state.Ct, state.B)
    state.D[state.N] = alpha
    state.D[state.N+1] = conj(alpha)
end

function init_triu(state::FactorizationType{T, Val{:DoubleShift}, Val{:NoPencil}, Tw}, decompose) where {T, Tw}
    ps = decompose(state.POLY)
    alpha =  R_factorization(ps, state.Ct, state.B)
end

# Here we have V = D Ct B; W = D1 Ct1 B1
function init_triu(state::FactorizationType{T, Val{:SingleShift}, Val{:HasPencil}, Tw}, decompose) where {T, Tw}

    ps, qs = decompose(state.POLY)
    beta =  R_factorization(qs, state.Ct1, state.B1)
    state.D1[state.N] = beta
    state.D1[state.N+1] = conj(beta)
    
    alpha =  R_factorization(ps, state.Ct, state.B)
    state.D[state.N] = alpha
    state.D[state.N+1] = conj(alpha)
end

# Here we have V = D Ct B; W = D1 Ct1 B1
function init_triu(state::FactorizationType{T, Val{:DoubleShift}, Val{:HasPencil}, Tw}, decompose) where {T, Tw}
    ps, qs = decompose(state.POLY)
    beta =  R_factorization(qs, state.Ct1, state.B1)
    alpha =  R_factorization(ps, state.Ct, state.B)
end





# # If there is an issue, this function can be used to resetart the algorithm
# # could be merged with init_state?
# function restart(state::ShiftType{T}) where {T}
#     # try again
#     init_state(state)
    
#     for i in 1:state.N
#         state.REIGS[i] = state.IEIGS[i] = zero(T)
#     end
#     state.ctrs.zero_index = 0
#     state.ctrs.start_index = 1
#     state.ctrs.stop_index = state.N - 1
#     state.ctrs.it_count = 0
#     state.ctrs.tr = state.N - 2
# end




