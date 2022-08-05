# The following tests are for functions in QRDecoders.jl/src/qrcoder.jl

@testset "Polynomials" begin    
    ## test Euclidean division of polynomials
    fdeg = rand(1:31) ## degree of message polynomial
    f = Poly([rand(1:255) for _ in 1:(fdeg+1)])
    n = rand(1:10) ## degree of generator polynomial
    g = generator(n)
    err = geterrorcorrection(f, n) ## error correction code
    rem = f << n % g ## remainder of euclidean division
    @test rem == err
end

@testset "Format decoding" begin
    ## test format encoding
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

    ## test fromat decoding
    fmt = rand(0:31)
    fmt_code = qrformat(fmt)
    de_fmt = qrdecode_format(fmt_code)
    @test fmt == de_fmt
    ### random disturbance
    for _ in 1:5
        fmt_code ⊻= 1 << rand(0:14)
    end
    de_fmt = qrdecode_format(fmt_code)
    fmt_dist = hamming_distance(fmt_code, qrformat(fmt))
    @test fmt_dist <= 5 && fmt == de_fmt
end


