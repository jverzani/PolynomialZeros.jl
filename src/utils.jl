#
"""
Try to convert a callable object into a polynomial, `Poly{T}`. `T` can be specified, or guessed from calling `f(0)`.
"""
as_poly(f::Poly) = f
    
function as_poly(f)
    T = typeof(f(0))
    as_poly(T, f)
end

function as_poly(T, f::Poly)
    convert(Poly{T}, f)
end


"""
as_poly{T}(T, f)

Convert `f` to a polynomial of type `Poly{T}`, or throw a `DomainError`.
"""
function as_poly(T, f)
    p = try
        x = variable(T)
        convert(Poly{T}, f(x))
    catch err
        throw(DomainError)
    end
    p
end
    

"""
reverse coefficients of a polynomial
"""
rcoeffs(p::Poly) = reverse(coeffs(p))

"""
monic
"""
monic(p::Poly) = p[degree(p)] != 0 ? Poly(p.a * inv(p[degree(p)]), p.var) : p




## replace with more robust method, but don't want to load Roots here
function _bisection_method(f, a, b)
    u,v = float(a), float(b)
    if u > v
        u,v = v,u
    end
    T = eltype(u)
    tol = 4*eps(T)

    fu, fv =  f(u), f(v)
    sign(fu) * sign(fv) > 0 && throw(DomainError) # not a bracket
    
    w = u + (v-u) * 0.5
    ctr = 1
    while norm(v - u) > tol
        fw = f(w)

        fw  == zero(T) && return w

        if sign(fu) * sign(fw) < 1
            u,v,fu,fv = u,w,fu,fw
        else
            u,v,fu,fv = w,v,fw,fv
        end
        w = u + (v-u) * 0.5
        ctr = ctr + 1; ctr > 100 && return w
    end
    w
end


