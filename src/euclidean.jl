# Extended Euclidean division algorithm
using QRDecoders: ReedSolomonError
using QRCoders.Polynomial: iszeropoly, degree, euclidean_divide, rstripzeros
using QRDecoders.Syndrome: erratalocator_polynomial, syndrome_polynomial, getpositions, polynomial_eval, forney_algorithm

"""
    extended_euclidean_divide(r₁::Poly, r₂::Poly)

Return polynomials u(x) and v(x) such that u(x)r₁(x) + v(x)r₂(x) = gcd(r₁(x), r₂(x)).

## illustration
Let rₖ = uₖr₁ + vₖr₂, k ≥ 2\\
    r₀ = q₀⋅r₁ + r₂ => r₂ = r₀ - q₀r₁, where r₀=r₂, q₀=0\\
    r₁ = q₁⋅r₂ + r₃ => r₃ = r₁ - q₁⋅r₂\\
    r₂ = q₂⋅r₃ + r₄ => r₄ = r₂ - q₂⋅r₃\\
    rₖ = qₖ⋅r_{k+1} + r_{k+2} => r_{k+2} = rₖ - qₖ⋅r_{k+1}\\
                              => r_{k+2} = uₖr₁ + vₖr₂ - qₖ(u_{k+1}r₁ + v_{k+1}r₂)\\
                                         = (uₖ - qₖ(u_{k+1})r₁ + (vₖ - qₖv_{k+1})r₂\\
                              => u_{k+1} = uₖ - qₖu_{k+1}\\
                                 v_{k+1} = vₖ - qₖv_{k+1}\\
    u₂, v₂ = 0, 1\\
    u₃, v₃ = 1, -q₁\\
    loop until\\
    rₜ = qₜ⋅r_{t+1} + 0 => r_{t+1} is the great common factor of r₁ and r₂
"""
function extended_euclidean_divide(r₁::Poly, r₂::Poly)
    u₁, v₁, u₂, v₂ = Poly.([[1], [0], [0], [1]])
    iszeropoly(r₂) && return u₁, v₁, r₁
    q, r₃ = euclidean_divide(r₁, r₂)
    # case r₁ == 0 => r₁ = 0⋅r₂ + 0 
    #              => r₃ == 0 => covered by the following loop
    while !iszeropoly(r₃) ## if rₖ == 0, then r_{k-1} = gcd(r₁, r₂)
        # @assert r₁ = u₁raw₁ + v₁raw₂
        # @assert r₂ = u₂raw₁ + v₂raw₂
        u₁, v₁, u₂, v₂ = u₂, v₂, u₁ + q * u₂, v₁ + q * v₂
        # @assert r₂ = u₁raw₁ + v₁raw₂
        # @assert r₃ = u₂raw₁ + v₂raw₂
        r₁, r₂ = r₂, r₃
        q, r₃ = euclidean_divide(r₁, r₂)
    end
    return u₂, v₂, r₂
end

"""
    Sugiyama_euclidean_divide(r₁::Poly, r₂::Poly, upperdeg::Int)

Yasuo Sugiyama's adaptation of the Extended Euclidean algorithm.
Find u(x), v(x) and r(x) s.t. r(x) = u(x)r₁(x) + v(x)r₂(x) where r(x) = gcd(r₁(x), r₂(x)) or deg(r(x)) ≤ upperdeg.
"""
function Sugiyama_euclidean_divide(r₁::Poly, r₂::Poly, upperdeg::Int)
    u₁, v₁, u₂, v₂ = Poly.([[1], [0], [0], [1]])
    iszeropoly(r₂) && return u₁, v₁, r₁
    q, r₃ = euclidean_divide(r₁, r₂)
    while degree(r₂) > upperdeg && !iszeropoly(r₃)
        # @assert r₁ = u₁raw₁ + v₁raw₂
        # @assert r₂ = u₂raw₁ + v₂raw₂
        u₁, v₁, u₂, v₂ = u₂, v₂, u₁ + q * u₂, v₁ + q * v₂
        # @assert r₂ = u₁raw₁ + v₁raw₂
        # @assert r₃ = u₂raw₁ + v₂raw₂
        r₁, r₂ = r₂, r₃
        q, r₃ = euclidean_divide(r₁, r₂)
    end
    return rstripzeros.((u₂, v₂, r₂))
end

"""
    euclidean_decoder(received::Poly, erasures::AbstractVector, nsym::Int)

Decode the received polynomial using the Euclidean algorithm(with erasures).
"""
euclidean_decoder(received::Poly, erasures::AbstractVector, nsym::Int) = euclidean_decoder!(copy(received), erasures, nsym)
function euclidean_decoder!(received::Poly, erasures::AbstractVector, nsym::Int)
    ## check data
    length(received) > 255 && throw(DomainError(received, "length of received polynomial must be less than 256"))
    length(erasures) > nsym && throw(ReedSolomonError())

    ## syndrome polynomial
    sydpoly = syndrome_polynomial(received, nsym)
    iszeropoly(sydpoly) && return received
    
    ## erasures locator polynomial Γx
    Γx = erratalocator_polynomial(erasures)
    xn = Poly(push!(zeros(Int, nsym), 1))
    
    ## deg(Ω(x))  ≤ ⌊(nsym + length(erasures)) / 2⌋ - 1
    upperdeg = (nsym + length(erasures)) ÷ 2 - 1
    ## error locator polynomial Λx
    ## evaluator polynomial evlpoly
    Λx, _, Ωx = Sugiyama_euclidean_divide(sydpoly * Γx, xn, upperdeg)
    
    ## errata locator polynomial
    errataloc = Λx * Γx
    
    ## errors positions + erasures positions
    errpos = vcat(getpositions(Λx), erasures)
    length(errpos) != length(errataloc) - 1 && throw(ReedSolomonError)
    
    ## Forney algorithm
    received.coeff[1 .+ errpos] .⊻= forney_algorithm(errataloc, Ωx, errpos)
    return received    
end

"""
    euclidean_decoder(received::Poly, erasures::AbstractVector, nsym::Int)

Decode the received polynomial using the Euclidean algorithm(without erasures).
"""
euclidean_decoder(received::Poly, nsym::Int) = euclidean_decoder!(copy(received), Int[], nsym)
