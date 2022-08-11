# Syndrome decoding
## outline
### 1. Compute the syndromes polynomial
### 2. Compute the erasure/error locator polynomial (from the syndromes)
### 3. Compute the erasure/error evaluator polynomial (from the syndromes and erasure/error locator polynomial)
### 4. Compute the erasure/error magnitude polynomial (from all 3 polynomials above)
### 5. Repair the input message simply by subtracting the magnitude polynomial from the input message.

module Syndrome

using QRCoders.Polynomial: mult, geterrorcorrection, Poly, gfpow2

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

Computes the syndrome polynomial S(x) for the message polynomial where the generator polynomial is of degree n.
"""
function syndrome_polynomial(msg::Poly, n::Int)
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
    ### e_k = \frac{2^{i_k}⋅Ω(2^{-i_k})}{Λ'(2^{-i_k})}
    errvals = Int[]
    for k in errpos
        num = mult(gfpow2(k), polynomial_eval(evlpoly, gfpow2(-k))) ## numerator
        den = polynomial_eval(errderi, gfpow2(-k)) ## denominator
        push!(errvals, divide(num, den))
    end
    ### correct errors
    recieved.coeff[1 .+ errpos] .⊻= errvals
    return recieved
end

end