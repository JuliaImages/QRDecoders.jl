# infomation of QR-code

## Hamming distance
"""
    hamming_weight(x::Poly)

Calculate the Hamming weight of a polynomial.
"""
hamming_weight(x::Poly) = hamming_weight(x.coeff)
hamming_weight(x::AbstractVector) = count(!iszero, x)
hamming_weight(x::Integer) = hamming_weight(digits(x; base=2))

"""
    hamming_distance(x, y)

Calculate the Hamming distance between two elements.
"""
hamming_distance(x::AbstractVector, y::AbstractVector) = hamming_weight(x - y)
hamming_distance(x::Integer, y::Integer) = hamming_weight(x ⊻ y)
hamming_distance(x::Poly, y::Poly) = hamming_weight(x.coeff - y.coeff)

## Functions about encoding & decoding
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
Bit data of information of quality.
"""
quality2binary = Dict(
    Low() => 0b01,
    Medium() => 0b00,
    Quartile() => 0b11,
    High() => 0b10)
binary2quality = Dict(val=>key for (key, val) in quality2binary)


"""
    qrdecode_format(fmt::Int)::Int

Decode format information.
"""
function qrdecode_format(fmt_code::Int)::Int
    fmt_info, best_dist = -1, 15
    for test_info in 0:31
        test_code = qrformat(test_info)
        test_dist = hamming_distance(fmt_code, test_code)
        if test_dist < best_dist
            best_dist = test_dist
            fmt_info = test_info
        elseif test_dist == best_dist
            fmt_info = -1
        end
    end
    return fmt_info
end