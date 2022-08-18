# Test for functions in QRDecoders.jl/src/qrinfo.jl

@testset "Format information" begin
    ## format encoding
    tag = true
    for ((qua, mask), val) in formatinfo
        fmt = qrformat(quality2binary[qua] << 3 ⊻ mask)
        val = parse(Int, join(Int.(val)); base=2)
        if fmt != val
            tag = false
            break
        end
    end
    @test tag

    ## fromat decoding
    fmt_info = rand(0:31)
    fmt_code = qrformat(fmt_info)
    fmt_decode_info = qrdecode_format(fmt_code)
    @test fmt_info == fmt_decode_info

    ## random disturbance
    ### errors with dectection capacity(≤5)
    fmt_code_disturb = fmt_code
    errors = unique!(rand(0:14, 5))
    for e in errors
        fmt_code_disturb ⊻= 1 << e
    end
    fmt_dist = hamming_distance(fmt_code, fmt_code_disturb)
    @test fmt_dist == length(errors)

    ## errors within correction capacity(≤2)
    fmt_code_disturb = fmt_code
    for _ in 1:2
        fmt_code_disturb ⊻= 1 << rand(0:14)
    end
    fmt_decode_info = qrdecode_format(fmt_code_disturb)
    @test fmt_info == fmt_decode_info
end

@testset "Hamming distance" begin
    ## distance of vectors
    message = rand(Int, 16)
    received = copy(message)
    errors = unique!(rand(1:16, 8))
    for e in errors
        received[e] += 1
    end
    @test hamming_distance(message, received) == length(errors)

    ## distance of binary numbers
    received = message = rand(1:2^16-1) ## length(message) <= 16
    erros = unique!(rand(0:15, 8))
    for e in errors
        received ⊻= 1 << e
    end
    @test hamming_distance(received, message) == length(errors)

    ## distance of polynomial
    msgpoly = randpoly(16)
    errpoly = randpoly(4)
    received = msgpoly + errpoly
    @test hamming_distance(msgpoly, received) == count(!iszero, errpoly.coeff)
end

