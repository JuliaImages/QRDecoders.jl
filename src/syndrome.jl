# Syndrome decoding
## outline
### 1. Compute the syndromes polynomial
### 2. Compute the erasure/error locator polynomial (from the syndromes)
### 3. Compute the erasure/error evaluator polynomial (from the syndromes and erasure/error locator polynomial)
### 4. Compute the erasure/error magnitude polynomial (from all 3 polynomials above)
### 5. Repair the input message simply by subtracting the magnitude polynomial from the input message.

module Syndrome

import QRCoders.Polynomial: mult, geterrorcorrection, Poly, logtable

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
    syndrome_polynomial(message::Poly, n::Int)

Computes the syndrome polynomial S(x) for the message polynomial where the generator polynomial is of degree n.
"""
function syndrome_polynomial(msg::Poly, n::Int)
    syndromes = polynomial_eval.(Ref(msg), getindex.(Ref(logtable), 0:(n-1)))
    return Poly(syndromes)
end

"""
    haserrors(msg::Poly, n::Int)

Returns true if the message polynomial has errors. (may go undetected when the number of errors exceeds n)
"""
haserrors(msg::Poly, n::Int) = !all(iszero, syndrome_polynomial(msg, n))

"""
    erratalocator_polynomial(errpos::AbstractVector)

Compute the erasures/error locator polynomial λ(x) from the erasures/errors positions.
"""
function erratalocator_polynomial(errpos::AbstractVector)
    isempty(errpos) && return Poly([1])
    return reduce(*, Poly([1, logtable[i]]) for i in errpos)
end

"""
    evaluator_polynomial(syd::Poly, errloc::Poly, n::Int)

Return the evaluator polynomial Ω(x) where Ω(x)≡S(x)λ(x) mod xⁿ.
"""
error_evaluator(syd::Poly, errloc::Poly, n::Int) = syd * errloc % Poly(push!(zeros(Int, n), 1))

"""
    correct_erased_errors(msg::Poly, errpos::AbstractVector, n::Int)

Forney algorithm, computes the values (error magnitude) to correct the input message.
"""
# function correct_erased_errors(msg::Poly, errpos::AbstractVector, n::Int)
#     
# end

end