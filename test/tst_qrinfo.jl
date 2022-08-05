# Test for functions in QRDecoders.jl/src/qrcoder.jl

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
    fmt = rand(0:31)
    fmt_code = qrformat(fmt)
    de_fmt = qrdecode_format(fmt_code)
    @test fmt == de_fmt

    ## random disturbance
    disturb_fmt = fmt_code
    for _ in 1:5
        disturb_fmt ⊻= 1 << rand(0:14)
    end
    de_fmt = qrdecode_format(disturb_fmt)
    fmt_dist = hamming_distance(fmt_code, disturb_fmt)
    @test fmt_dist ≤ 5 && fmt == de_fmt
end

@testset "Hamming distance" begin
    ## distance of vectors
    message = rand(Int, 16)
    recieved = copy(message)
    errors = unique!(rand(1:16, 8))
    for e in errors
        recieved[e] += 1
    end
    @test hamming_distance(message, recieved) == length(errors)

    ## distance of binary numbers
    recieved = message = rand(1:2^16-1) ## length(message) <= 16
    erros = unique!(rand(0:15, 8))
    for e in errors
        recieved ⊻= 1 << e
    end
    @test hamming_distance(recieved, message) == length(errors)
end

