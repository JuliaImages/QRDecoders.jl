# Syndrome decoding
## outline
### 1. Compute the syndromes polynomial S(x) = ∑ᵢsᵢxⁱ
### 2. Compute the erasure/error locator polynomial Λ(x) = ∏ₖ(2^{iₖ}⋅xᵏ - 1) = ∑ᵢλᵢxⁱ
### 3. Compute the erasure/error evaluator polynomial Ω(x) = S(x)Λ(x) mod xⁿ
### 4. Compute the erasure/error magnitude polynomial eₖ = (2^{i_k}⋅Ω(2^{-i_k})) / (Λ'(2^{-i_k}))
### 5. Repair the input message simply by subtracting the magnitude polynomial from the input message.

module Syndrome

using QRCoders.Polynomial: mult, Poly, gfpow2, gflog2, gfinv, divide, iszeropoly, unit, rstripzeros
using QRDecoders: ReedSolomonError

"""
    polynomial_eval(p::Poly, x::Int)

Evaluates the polynomial p(x) at x in GF256.
"""
function polynomial_eval(p::Poly, x::Int)
    val = last(p.coeff) ## leading term
    coeff = @view(p.coeff[end-1:-1:1])
    return foldl((i, j) -> mult(i, x) ⊻ j, coeff; init=val)
end

"""
    derivative_polynomial(p::Poly)

Computes the derivative of the polynomial p.
"""
function derivative_polynomial(p::Poly)
    lp = length(p)
    lp == 1 && return zero(Poly) ## constant polynomial
    coeff = Vector{Int}(undef, lp - 1)
    for (i, c) in enumerate(@view(p.coeff[2:end]))
        coeff[i] =  isodd(i) ? c : 0 ## i⋅x = 0 if i is even
    end
    return Poly(coeff)
end

"""
    syndrome_polynomial(message::Poly, n::Int)

Computes the syndrome polynomial S(x) for the message polynomial where n is the number of syndromes.
"""
function syndrome_polynomial(msg::Poly, n::Int)
    ## S(x) = ∑ᵢsᵢxⁱ where sᵢ = R(2ⁱ) = E(2ⁱ)
    syndromes = polynomial_eval.(Ref(msg), gfpow2.(0:(n-1)))
    return Poly(syndromes)
end

"""
    haserrors(msg::Poly, n::Int)

Returns true if the message polynomial has errors. (may go undetected when the number of errors exceeds n)
"""
haserrors(msg::Poly, n::Int) = !iszeropoly(syndrome_polynomial(msg, n))

"""
    erratalocator_polynomial(errpos::AbstractVector)

Compute the erasures/error locator polynomial Λ(x) from the erasures/errors positions.
"""
function erratalocator_polynomial(errpos::AbstractVector)
    isempty(errpos) && return unit(Poly)
    return reduce(*, Poly([1, gfpow2(i)]) for i in errpos)
end

"""
    evaluator_polynomial(syd::Poly, errloc::Poly, n::Int)

Return the evaluator polynomial Ω(x) where Ω(x)≡S(x)Λ(x) mod xⁿ.
"""
evaluator_polynomial(syd::Poly, errloc::Poly, n::Int) = Poly((syd * errloc).coeff[1:n])

"""
    fillearsed(msg::Poly, errpos::AbstractVector, n::Int)

Forney algorithm, computes the values (error magnitude) to correct the input message.
"""
fillearsed(recieved::Poly, errpos::AbstractVector, n::Int) = fillearsed!(copy(recieved), errpos, n)
function fillearsed!(recieved::Poly, errpos::AbstractVector, n::Int)
    ## number of errors exceeds limitation of the RS-code
    length(errpos) > n && throw(ReedSolomonError())

    ## syndrome polynomial S(x)
    sydpoly = syndrome_polynomial(recieved, n) 
    
    ## error locator polynomial
    errloc = erratalocator_polynomial(errpos)
    ## derivative of error locator polynomial
    errderi = derivative_polynomial(errloc)
    
    ## evaluator polynomial Ω(x)≡S(x)Λ(x) mod xⁿ
    evlpoly = evaluator_polynomial(sydpoly, errloc, n)
    
    ## computes error magnitudes using Forney algorithm
    ### e_k = \frac{2^{i_k}⋅Ω(2^{-i_k})} {Λ'(2^{-i_k})}
    forneynum(k::Int) = mult(gfpow2(k), polynomial_eval(evlpoly, gfpow2(-k))) ## numerator
    forneyden(k::Int) = polynomial_eval(errderi, gfpow2(-k)) ## denominator
    errvals = @. divide(forneynum(errpos), forneyden(errpos))
    ### correct errors
    recieved.coeff[1 .+ errpos] .⊻= errvals
    return recieved
end

"""
    reducebyHorner(p::Poly, a::Int)

Horner's rule, find the polynomial q(x) such that p(x)-(x-a)q(x) is a constant polynomial.
"""
function reducebyHorner(p::Poly, a::Int)
    reduced = Vector{Int}(undef, length(p))
    reduced[1] = val = last(p.coeff)
    coeff = @view(p.coeff[end-1:-1:1])
    for (i, c) in enumerate(coeff)
        val = mult(val, a) ⊻ c
        reduced[i + 1] = val
    end
    return Poly(reverse!(reduced))
end

"""
    findroots(p::Poly)

Computes the roots of the polynomial p using Horner's method. Returns a empty list if p(x) contains duplicate roots or roots not in GF(256).
"""
function findroots(p::Poly)
    n = length(p) - 1
    roots = Vector{Int}(undef, n)
    for r in 0:255
        reducepoly = reducebyHorner(p, r)
        ## decrease the degree of the polynomial by 1
        popfirst!(reducepoly.coeff) == 0 || continue
        roots[n] = r
        n, p = n - 1, reducepoly
        n == 0 && return reverse!(roots)
    end
    ## p(x) contains duplicate roots or roots not in GF(256)
    return Int[] 
end

"""
    getposition(Λx::Poly)

Recover positions of errors from the error locator polynomial. Note that Λx = \\Prod_{k=1}^{v}(2^{i_k}⋅x + 1).
"""
getposition(Λx::Poly) = mod.(-gflog2.(findroots(Λx)), 255)

"""
    erratalocator_polynomial(recieved::Poly, n::Int)

Berlekamp-Massey algorithm, compute the error locator polynomial Λ(x).
"""
function erratalocator_polynomial(recieved::Poly, nsym::Int)
    ## syndromes
    S = syndrome_polynomial(recieved, nsym).coeff

    ## initialize
    L = 0 # number of errors
    x = Poly([0, 1])
    Λx = unit(Poly) # error locator polynomial
    Bx = copy(Λx) # copy of the last Λx since L is updated

    ## discrepancy Δᵣ = Λ₀Sᵣ₋₁ + Λ₁Sᵣ₋₂ + ⋯ 
    getdelta(r) = @views reduce(⊻, mult(i, j) for (i, j) in zip(S[r:-1:1], Λx.coeff))

    ## iteration
    for r in 1:nsym
        Δ = getdelta(r)
        if Δ == 0 || 2 * L > r - 1 # δ = 0
            Λx, Bx = Λx + Δ * x * Bx, x * Bx
        else # δ = 1
            L = r - L
            Λx, Bx = Λx + Δ * x * Bx, gfinv(Δ) * Λx
        end
    end
    Λx = rstripzeros(Λx)

    ## number of errors exceeds limitation of the RS-code
    if iszeropoly(Λx) || length(Λx) - 1 > nsym ÷ 2 
        throw(ReedSolomonError())
    end
    return Λx
end

"""
    erratalocator_polynomial(recieved::Poly, erasures::AbstractVector, n::Int)

Berlekamp-Massey algorithm, compute the error locator polynomial Λ(x) given the erasures positions.
"""
# function erratalocator_polynomial(recieved::Poly, erasures::AbstractVector, nsym::Int)
    
# end

end