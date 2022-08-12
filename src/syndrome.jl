# Syndrome decoding
## outline
### 1. Compute the syndromes polynomial S(x) = ∑ᵢsᵢxⁱ
### 2. Compute the erasure/error locator polynomial Λ(x) = ∏ₖ(2^{iₖ}⋅xᵏ - 1) = ∑ᵢλᵢxⁱ
### 3. Compute the erasure/error evaluator polynomial Ω(x) = S(x)Λ(x) mod xⁿ
### 4. Compute the erasure/error magnitude polynomial eₖ = (2^{i_k}⋅Ω(2^{-i_k})) / (Λ'(2^{-i_k}))
### 5. Repair the input message simply by subtracting the magnitude polynomial from the input message.

module Syndrome

using QRCoders.Polynomial: mult, geterrorcorrection, Poly, gfpow2, divide
using QRDecoders: ReedSolomonError

"""
    polynomial_eval(p::Poly, x::Int)

Evaluates the polynomial p(x) at x in GF256.
"""
function polynomial_eval(p::Poly, x::Int)
    coeff = reverse(p.coeff)
    val = popfirst!(coeff)
    for c in coeff
        val = mult(val, x) ⊻ c
    end
    return val
end

"""
    derivative_polynomial(p::Poly)

Computes the derivative of the polynomial p.
"""
function derivative_polynomial(p::Poly)
    length(p) == 1 && return Poly([0]) ## constant polynomial
    coeff = Int[]
    for (i, c) in enumerate(@view(p.coeff[2:end]))
        push!(coeff, isodd(i) ? c : 0) ## i⋅x = 0 if i is even
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
haserrors(msg::Poly, n::Int) = !all(iszero, syndrome_polynomial(msg, n))

"""
    erratalocator_polynomial(errpos::AbstractVector)

Compute the erasures/error locator polynomial Λ(x) from the erasures/errors positions.
"""
function erratalocator_polynomial(errpos::AbstractVector)
    isempty(errpos) && return Poly([1])
    return reduce(*, Poly([1, gfpow2(i)]) for i in errpos)
end

"""
    evaluator_polynomial(syd::Poly, errloc::Poly, n::Int)

Return the evaluator polynomial Ω(x) where Ω(x)≡S(x)Λ(x) mod xⁿ.
"""
evaluator_polynomial(syd::Poly, errloc::Poly, n::Int) = Poly((syd * errloc).coeff[1:n])

"""
    syndrome_decoder(msg::Poly, errpos::AbstractVector, n::Int)

Forney algorithm, computes the values (error magnitude) to correct the input message.
"""
syndrome_decoder(recieved::Poly, errpos::AbstractVector, n::Int) = syndrome_decoder!(copy(recieved), errpos, n)
function syndrome_decoder!(recieved::Poly, errpos::AbstractVector, n::Int)
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
    coeff = reverse(p.coeff)
    val = popfirst!(coeff)
    reduced = [val]
    for c in coeff
        val = mult(val, a) ⊻ c
        push!(reduced, val)
    end
    return Poly(reverse!(reduced))
end

"""
    findroots(p::Poly)

Computes the roots of the polynomial p using Horner's method. Returns a empty list if p(x) contains duplicate roots or roots not in GF(256).
"""
function findroots(p::Poly)
    n = length(p) - 1
    roots = Int[]
    for r in 0:255
        reducepoly = reducebyHorner(p, r)
        popfirst!(reducepoly.coeff) == 0 || continue
        push!(roots, r)
        n, p = n - 1, reducepoly
        n == 0 && return roots
    end
    ## p(x) contains duplicate roots or roots not in GF(256)
    return Int[] 
end

end