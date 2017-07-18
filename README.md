# PolynomialZeros

Methods to find zeros (roots) of polynomials over given domains

[![Build Status](https://travis-ci.org/jverzani/PolynomialZeros.jl.svg?branch=master)](https://travis-ci.org/jverzani/PolynomialZeros.jl)

[![Coverage Status](https://coveralls.io/repos/jverzani/PolynomialZeros.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/jverzani/PolynomialZeros.jl?branch=master)

[![codecov.io](http://codecov.io/github/jverzani/PolynomialZeros.jl/coverage.svg?branch=master)](http://codecov.io/github/jverzani/PolynomialZeros.jl?branch=master)


This package provides the method `poly_roots` to find roots of
univariate polynomial functions over the complex numbers, the real
numbers, the rationals, the integers, and $Z_p$. (A "root" is the name
for a "zero" of a polynomial.) The package takes advantage of many
other root-finding packages for polynomials within Julia (e.g.,
`PolynomialRoots` for numeric solutions over the complex numbers and
`PolynomialFactors` or exact solutions over the intergers).

The basic interface is

```
poly_roots(f, domain)
```

Where `f` is in `Poly{T}` (from the `Polynomials.jl` package) or can be converted into `Poly{T}`. The domain is specified by `Over.C`, `Over.R`, `Over.Q`, `Over.Z`, or `over.Zp{p}`. Not all polynomials will have such a factorization.


Examples:

```
julia> poly_roots(x -> x^4 - 1, Over.C)  # uses `roots` from `PolynomialRoots.jl`
4-element Array{Complex{Float64},1}:
  0.0+1.0im
  1.0-0.0im
  0.0-1.0im
 -1.0+0.0im


julia> poly_roots(x -> x^4 - 1, Over.R)  
2-element Array{Float64,1}:
  1.0
  -1.0
  
julia> poly_roots(x -> x^4 - 1, Over.Q) # uses `PolynomialFactors.jl`
2-element Array{Rational{Int64},1}:
 -1//1
  1//1

julia> poly_roots(x -> x^4 - 1, Over.Z) # uses `PolynomialFactors.jl`
2-element Array{Int64,1}:
 -1
  1

julia> poly_roots(x -> x^4 - 1, Over.Zp{5}) # uses `PolynomialFactors.jl`
4-element Array{Int64,1}:
 4
 1
 3
 2
```

Domains can also have their types specified. For example, to solve
over the `BigFloat` type, we have:

```julia
poly_roots(x -> x^4 - 1, Over.CC{BigFloat})  # `CC{BigFloat}` not just `C`
```

## Details


This package uses:

* The `PolynomialRoots` package to find roots over the
complex numbers. The `Roots` package can also be used.

* The `PolynomialFactors` package to return roots over the
rationals, integers, and integers module a prime.

* As well, it provides an algorithm to find the real
roots of polynomials that originally was in the `Roots` package.


The main motivation for this package was to move the polynomial
specific code out of the `Roots` package. This makes the `Roots`
package have fewer dependencies and a more focused task. In addition,
the polynomial specific code could use some better implementations of
the underlying algorithms.

In the process of doing this, making a common interface to the other
root-finding packages seemed to make sense.

### Other possibly useful methods

The package also provides

* `PolynomialZeros.agcd` for computing an *approximate* GCD of
  polynomials `p` and `q` over `Poly{Float64}`. (This is used to
  reduce a polynomial over the reals to a square-free
  polynomial. Square-free polynomials are needed for the
  algorithm used. However, the `agcd` function begins to be inaccurate for
  higher degree polynomials, so this is not the default behaviour, rather it must be requested.)

* `PolynomialZeros.multroot` for finding roots of `p` in
  `Poly{Float64}` over `Complex{Float64}` which has some advantage if
  `p` has high multiplicities. The `roots` function from the
  `Polynomials` package will find all the roots of a polynomial. Its
  performance degrades when the polynomial has high
  multiplicities. The multroot function is provided to handle this
  case a bit better. The function follows algorithms due to Zeng,
  "Computing multiple roots of inexact polynomials", Math. Comp. 74
  (2005), 869-903.

