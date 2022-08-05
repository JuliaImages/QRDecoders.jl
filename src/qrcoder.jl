# Add Functions to QRCode
## Functions about polynomial
import Base: %, ÷

"""
    %(f::Poly, g::Poly)

Remainder of Euclidean division.
"""
%(f::Poly, g::Poly) = Poly((f + (f ÷ g) * g).coeff[1:(length(g) - 1)])

"""
    ÷(f::Poly, g::Poly)

Quotient of Euclidean division.
"""
function ÷(f::Poly, g::Poly)
    quodeg = length(f) - length(g) ## degree of the quotient polynomial
    quodeg < 0 && return Poly([0])
    g <<= quodeg
    gn = lead(g) ## leading term of g(x)
    quocoef = Vector{Int}(undef, quodeg + 1)
    for i in 1:quodeg
        quocoef[i] = divide(lead(f), gn)
        f = init!(quocoef[i] * g + f)
        tail!(g)
    end
    quocoef[end] = divide(lead(f), gn)
    return Poly(reverse!(quocoef))
end

"""
    divide(a::Int, b::Int)

Division of intergers in GF(256).
"""
function divide(a::Int, b::Int)
    b == 0 && throw(DivideError())
    a == 0 && return 0
    xa, xb = antilogtable[a], antilogtable[b]
    ## use mod instead of % to avoid negative numbers
    return logtable[mod(xa - xb, 255)]
end


## Functions about encoding & decoding
"""
    hamming_weight(x::Poly)

Calculate the Hamming weight of a polynomial.
"""
hamming_weight(x::Poly) = hamming_weight(x.coeff)
hamming_weight(x::Int) = hamming_weight(digits(x; base=2))
hamming_weight(x::AbstractVector) = count(!=(0), x)


"""
    hamming_distance(x::Poly, y::Poly)

Calculate the Hamming distance between two polynomials.
"""
hamming_distance(x::Poly, y::Poly) = hamming_distance(x.coeff, y.coeff)
hamming_distance(x::Int, y::Int) = hamming_distance(digits(x; base=2), digits(y; base=2))
hamming_distance(x::AbstractVector, y::AbstractVector) = hamming_weight(x - y)


"""
    qrformat(fmt::Int)

Generate standard format infomation (format + error correction + mask).
"""
function qrformat(fmt::Int)
    err = fmt << 10
    g = 0x537 ## generator polynomial(= 0b10100110111 in binary)
    for i in 4:-1:0
        if !iszero(err & (1 << (i + 10)))
            err ⊻= g << i
        end
    end
    fmt << 10 ⊻ err ⊻ 0x5412 ## mask(= 0b101010000010010 in binary)
end
qrformat(fmt::Integer) = qrformat(Int(fmt)) ## to avoid integer overflow

"""
Binary code of information of quality.
"""
quality2binary = Dict(
    Low() => 0b01,
    Medium() => 0b00,
    Quartile() => 0b11,
    High() => 0b10)
binary2quality = Dict(val=>key for (key, val) in quality2binary)


"""
    qrdecode_format(fmt::Int)

Decode format information.
"""
function qrdecode_format(fmt::Int)
    best_fmt, best_dist = -1, 15
    for test_fmt in 0:31
        test_code = qrformat(test_fmt)
        test_dist = hamming_weight(fmt ⊻ test_code)
        if test_dist < best_dist
            best_dist = test_dist
            best_fmt = test_fmt
        end
    end
    ## best_dist > 5 && return -1
    return best_fmt
end