## Types

## A container for our counters
#mutable struct AMVW_Counte
type AMVW_Counter
    zero_index::Int
    start_index::Int
    stop_index::Int
    it_count::Int
    tr::Int
end


## Rotators

## Our rotators have field names c, s where c nad s are either T or Complex{T}
@compat abstract type CoreTransform{T} end
@compat abstract type Rotator{T} <: CoreTransform{T} end

is_diagonal{T}(r::Rotator{T}) = norm(r.s) <= eps(T)


Base.copy(a::Rotator) = Rotator(a.c, a.s, a.i)
function Base.copy!(a::Rotator, b::Rotator)
    vals!(a, vals(b)...)
    idx!(a, idx(b))
end

## set values
vals{T}(r::Rotator{T}) = (r.c, r.s)
idx(r::Rotator) = r.i
idx!(r::Rotator, i::Int) = r.i = i




#the index is superflous for now, and a bit of a hassle to keep immutable
#but might be of help later if twisting is approached. Shouldn't effect speed, but does mean 3N storage (Q, Ct, B)
#so may be
#
#mutable struct RealRotator{T} <: Rotator{T}
type RealRotator{T} <: Rotator{T}
c::T
s::T
i::Int
end

function Base.ctranspose(r::RealRotator)
    RealRotator(r.c, -r.s, r.i)
end


Base.one{T}(::Type{RealRotator{T}}) = RealRotator(one(T), zero(T), 0)
#Base.ones{T}(S::Type{RealRotator{T}}, N) = [one(S) for i in 1:N]
_ones{T}(S::Type{RealRotator{T}}, N) = [one(S) for i in 1:N]

## set values
function vals!{T}(r::RealRotator, c::T, s::T)
    # normalize in case of roundoff errors
    # but, using hueristic of 6.3 on square roots
    
    nrmi = sqrt(c^2 + s^2 )
    nrmi = norm(nrmi - one(T)) >= 1e2*eps(T) ? inv(nrmi) : one(T)
    r.c = c * nrmi
    r.s = s * nrmi
end

##################################################
### Okay, now try with complex C, real s
    
#mutable struct ComplexRealRotator{T} <: Rotator{T}
type ComplexRealRotator{T} <: Rotator{T}    
c::Complex{T}
s::T
i::Int
end

function Base.ctranspose(r::ComplexRealRotator)
    ComplexRealRotator(conj(r.c), -r.s, r.i)
end

function vals!{T}(r::ComplexRealRotator, c::Complex{T}, s::T)
    # normalize in case of roundoff errors
    # but, using hueristic on 6.3 on square roots
    
    nrmi = sqrt(abs(c * conj(c) + s^2))
    nrmi = norm(nrmi - one(T)) >= eps(T) ? inv(nrmi) : one(T)
    r.c = c * nrmi
    r.s = s * nrmi
end
function vals!{T}(r::ComplexRealRotator, c::Complex{T}, s::Complex{T})
##    abs(imag(s)) < 4eps(T) || error("setting vals needs real s, got $s")
    vals!(r, c, real(s))
end
vals!{T}(r::ComplexRealRotator{T}, c::T, s::T) = vals!(r, complex(c, zero(T)), s)

Base.one{T}(::Type{ComplexRealRotator{T}}) = ComplexRealRotator(complex(one(T), zero(T)), zero(T), 0)
#Base.ones{T}(S::Type{ComplexRealRotator{T}}, N) = [one(S) for i in 1:N]
_ones{T}(S::Type{ComplexRealRotator{T}}, N) = [one(S) for i in 1:N]


Base.copy(a::ComplexRealRotator) = ComplexRealRotator(a.c, a.s, a.i)
function Base.copy!(a::ComplexRealRotator, b::ComplexRealRotator)
    vals!(a, vals(b)...)
    idx!(a, idx(b))
end
  



##################################################
## We use two complex, rather than 3 reals here.
## Will be basically the ame storage, as we don't need to include a D, but not quite (12N, not 11N)
    
#mutable struct ComplexComplexRotator{T} <: Rotator{T}
type ComplexComplexRotator{T} <: Rotator{T}
c::Complex{T}
s::Complex{T}
i::Int
end

function Base.ctranspose(r::ComplexComplexRotator)
    ComplexComplexRotator(conj(r.c), -r.s, r.i)
end


Base.one{T}(::Type{ComplexComplexRotator{T}}) = ComplexComplexRotator(complex(one(T), zero(T)), complex(zero(T), zero(T)), 0)
#Base.ones{T}(S::Type{ComplexComplexRotator{T}}, N) = [one(S) for i in 1:N]
_ones{T}(S::Type{ComplexComplexRotator{T}}, N) = [one(S) for i in 1:N]

## set values
function vals!{T}(r::ComplexComplexRotator, c::Complex{T}, s::Complex{T})
    # normalize in case of roundoff errors
    # but, using hueristic on 6.3 on square roots
    
    nrmi = sqrt(abs(c * conj(c) + s * conj(s)))
    nrmi = norm(nrmi - one(T)) >= eps(T) ? inv(nrmi) : one(T)
    r.c = c * nrmi
    r.s = s * nrmi
end
vals!{T}(r::ComplexComplexRotator, c::Complex{T}, s::T) = vals!(r, c, complex(s,zero(T)))
vals!{T}(r::ComplexComplexRotator{T}, c::T, s::T) = vals!(r, complex(c, zero(T)), complex(s, zero(T)))



##################################################
## Factorization Types ##
##################################################


# N -- degree of poly
# POLY -- [p0, p1, p2, ..., pn] (not reversed, same as p.a for Poly type
# reverse in the init_state

@compat abstract type FactorizationType{T, ShiftType, Pencil, Twisted} end

## SingleShift(:SingleShift) -- Q, D, U, Ut
## DoubleShift(:DoubleShift) -- Q, U, V, Ut, Vt, W

## :NoPencil -- Ct, B
## :HasPencil -- [D], Ct, B, Ct1, B1, [D1]

## :NotTwisted --
## :IsTwisted -- sigma

## All have: N, POLY, REIGS, IEIGS,
## A, Bk, R, e1, e2, ctr

############### No Pencil, Not Twisted ###################################

## RDS, no pencil, not twisted
#struct Real_DoubleShift_NoPencil_NotTwisted{T} <: FactorizationType{T, Val{:DoubleShift}, Val{:NoPencil}, Val{:NotTwisted}
immutable Real_DoubleShift_NoPencil_NotTwisted{T} <: FactorizationType{T, Val{:DoubleShift}, Val{:NoPencil}, Val{:NotTwisted}}

N::Int
POLY::Vector{T}
##
Q::Vector{RealRotator{T}}
Ct::Vector{RealRotator{T}}  # We use C', not C here
B::Vector{RealRotator{T}}
##
REIGS::Vector{T}
IEIGS::Vector{T}
## reusable storage
U::RealRotator{T}
Ut::RealRotator{T}
V::RealRotator{T}
Vt::RealRotator{T}
W::RealRotator{T}
A::Matrix{T}    # for parts of A = QR
R::Matrix{T}    # temp storage, sometimes R part of QR
e1::Vector{T}   # eigen values e1, e2
e2::Vector{T}
ctrs::AMVW_Counter
end

function Base.convert{T <: AbstractFloat}(::Type{FactorizationType{T, Val{:DoubleShift}, Val{:NoPencil}, Val{:NotTwisted}}}, ps::Vector{T})

    N = length(ps) - 1

    Real_DoubleShift_NoPencil_NotTwisted(N, ps,
                                         _ones(RealRotator{T}, N), #Q
                                         _ones(RealRotator{T}, N), #Ct
                                         _ones(RealRotator{T}, N), #B
                                         zeros(T, N),  zeros(T, N), #EIGS
                                         one(RealRotator{T}), one(RealRotator{T}), #U, Ut
                                         one(RealRotator{T}), one(RealRotator{T}), #V, Vt
    one(RealRotator{T}), # W
    zeros(T, 2, 2),zeros(T, 3, 2), # A, R
    zeros(T,2), zeros(T,2),
    AMVW_Counter(0,1,N-1, 0, N-2)
    )
end


# ComplexReal Double Shift, no pencil, not twisted

#mutable struct ComplexReal_SingleShift_NoPencil_NotTwisted{T} <: FactorizationType{T, Val{:SingleShift}, Val{:NoPencil}, Val{:NotTwisted}}
type ComplexReal_SingleShift_NoPencil_NotTwisted{T} <: FactorizationType{T, Val{:SingleShift}, Val{:NoPencil}, Val{:NotTwisted}}

N::Int
POLY::Vector{Complex{T}}
Q::Vector{ComplexRealRotator{T}}
D::Vector{Complex{T}}
Ct::Vector{ComplexRealRotator{T}}  # We use C', not C here
B::Vector{ComplexRealRotator{T}}
#
REIGS::Vector{T}
IEIGS::Vector{T}
# reusable storage
U::ComplexRealRotator{T}
Ut::ComplexRealRotator{T}
A::Matrix{Complex{T}}    # for parts of A = QR
R::Matrix{Complex{T}}    # temp storage, sometimes R part of QR
e1::Vector{T}   # eigen values e1, e2, store as (re,imag)
e2::Vector{T}
ray::Bool
ctrs::AMVW_Counter
end

function Base.convert{T}(::Type{FactorizationType{T, Val{:SingleShift}, Val{:NoPencil}, Val{:NotTwisted}}}, ps::Vector{Complex{T}}; ray=true)

    N = length(ps) - 1
    
    ComplexReal_SingleShift_NoPencil_NotTwisted(N, ps,
                           _ones(ComplexRealRotator{T}, N), #Q
                           ones(Complex{T}, N+1), # D
                           _ones(ComplexRealRotator{T}, N), #Ct
                           _ones(ComplexRealRotator{T}, N), #B
                           zeros(T, N),  zeros(T, N), #EIGS
    one(ComplexRealRotator{T}), one(ComplexRealRotator{T}), #U, Ut
    zeros(Complex{T}, 2, 2),zeros(Complex{T}, 3, 2), # A R
    zeros(T,2), zeros(T,2),
    #    true,  # true for Wilkinson, false for Rayleigh.Make adjustable!
    ray,
    AMVW_Counter(0,1,N-1, 0, N-2)
    )
end



############## Has Pencil, Not twisted ####################################

## RDS, no pencil, not twisted
#struct Real_DoubleShift_HasPencil_NotTwisted{T} <: FactorizationType{T, Val{:DoubleShift}, Val{:HasPencil}, Val{:NotTwisted}}
immutable Real_DoubleShift_HasPencil_NotTwisted{T} <: FactorizationType{T, Val{:DoubleShift}, Val{:HasPencil}, Val{:NotTwisted}}

N::Int
POLY::Vector{T}
##
Q::Vector{RealRotator{T}}
Ct::Vector{RealRotator{T}}  # We use C', not C here
B::Vector{RealRotator{T}}
#
Ct1::Vector{RealRotator{T}}  # W = Ct * B
B1::Vector{RealRotator{T}}   # 
##
REIGS::Vector{T}
IEIGS::Vector{T}
## reusable storage
U::RealRotator{T}
Ut::RealRotator{T}
V::RealRotator{T}
Vt::RealRotator{T}
W::RealRotator{T}
#
A::Matrix{T}    # for parts of A = QR
R::Matrix{T}    # temp storage, sometimes R part of QR
e1::Vector{T}   # eigen values e1, e2
e2::Vector{T}
ctrs::AMVW_Counter
end

function Base.convert{T}(::Type{FactorizationType{T, Val{:DoubleShift}, Val{:HasPencil}, Val{:NotTwisted}}}, ps::Vector{T})

    N = length(ps) - 1

    Real_DoubleShift_HasPencil_NotTwisted(N, ps,
                                          _ones(RealRotator{T}, N), #Q
                                          _ones(RealRotator{T}, N), #Ct
                                          _ones(RealRotator{T}, N), #B
                                          _ones(RealRotator{T}, N), #Ct1
                                          _ones(RealRotator{T}, N), #B1                                         
    zeros(T, N),  zeros(T, N), #EIGS
    one(RealRotator{T}), one(RealRotator{T}),
    one(RealRotator{T}), one(RealRotator{T}),
    one(RealRotator{T}), #U,U',V,V',W
    zeros(T, 2, 2), zeros(T, 3, 2), # A R
    zeros(T,2), zeros(T,2),
    AMVW_Counter(0,1,N-1, 0, N-2)
    )
end


# ComplexReal Double Shift, has pencil, not twisted

#mutable struct ComplexReal_SingleShift_HasPencil_NotTwisted{T} <: FactorizationType{T, Val{:SingleShift}, Val{:HasPencil}, Val{:NotTwisted}}
type ComplexReal_SingleShift_HasPencil_NotTwisted{T} <: FactorizationType{T, Val{:SingleShift}, Val{:HasPencil}, Val{:NotTwisted}}

N::Int 
POLY::Vector{Complex{T}}
Q::Vector{ComplexRealRotator{T}}
D::Vector{Complex{T}}
Ct::Vector{ComplexRealRotator{T}}  # We use C', not C here
B::Vector{ComplexRealRotator{T}}
# pencil part
D1::Vector{Complex{T}}
Ct1::Vector{ComplexRealRotator{T}}  # not inverses
B1::Vector{ComplexRealRotator{T}}
#
REIGS::Vector{T}
IEIGS::Vector{T}
## reusable storage
U::ComplexRealRotator{T}
Ut::ComplexRealRotator{T}
#
A::Matrix{Complex{T}}    # for parts of A = QR
R::Matrix{Complex{T}}    # temp storage, sometimes R part of QR
e1::Vector{T}   # eigen values e1, e2, store as (re,imag)
e2::Vector{T}
#
ray::Bool
ctrs::AMVW_Counter
end

function Base.convert{T}(::Type{FactorizationType{T, Val{:SingleShift}, Val{:HasPencil}, Val{:NotTwisted}}}, ps::Vector{Complex{T}}; ray=true)

    N = length(ps) - 1
    
    ComplexReal_SingleShift_HasPencil_NotTwisted(N, ps,
                                                 _ones(ComplexRealRotator{T}, N), #Q
                                                 ones(Complex{T}, N+1), # D
                                                 _ones(ComplexRealRotator{T}, N), #Ct
                                                 _ones(ComplexRealRotator{T}, N), #B
                                                 ones(Complex{T}, N+1), # D1
    _ones(ComplexRealRotator{T}, N), #Ct1
    _ones(ComplexRealRotator{T}, N), #B1
    zeros(T, N),  zeros(T, N), #EIGS
    one(ComplexRealRotator{T}), one(ComplexRealRotator{T}), #U, Ut
    zeros(Complex{T}, 2, 2),zeros(Complex{T}, 3, 2), # A , R
    zeros(T,2), zeros(T,2),
    #    true,  # true for Wilkinson, false for Rayleigh.Make adjustable!
    ray,
    AMVW_Counter(0,1,N-1, 0, N-2)
    )
end

