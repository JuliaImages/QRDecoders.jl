# Test for functions in QRDecoders.jl/src/qrinfo.jl

@testset "Version information" begin
    ## version encoding
    tag = true
    for (v, val) in enumerate(versioninfo)
        v < 7 && continue
        ver_code = qrversion(v)
        val = parse(Int, reverse(join(Int.(val))); base=2)
        if ver_code != val || qrdecode_version(ver_code) != v
            tag = false
            break
        end
    end
    @test tag
    @test qrversion.(7:40) == qrversion.(0x07:0x28)
    @test_throws InfoError qrversion(6)
    @test_throws InfoError qrversion(41)
    ## distance of the Version-Code
    @test hamming_distance(qrversion.(7:40)) == 8

    ## version decoding -- with disturbance
    ver_info = rand(7:40)
    ver_code = qrversion(ver_info)
    ver_code_disturb = ver_code
    errors = sample(7:40, 8; replace=false)
    for e in errors
        ver_code_disturb ⊻= 1 << e
    end
    ver_dist = hamming_distance(ver_code, ver_code_disturb)
    @test ver_dist == length(errors)
    ## errors within correction capacity(≤4)
    ver_code_disturb = ver_code
    for _ in 1:4
        ver_code_disturb ⊻= 1 << rand(7:40)
    end
    ver_decode_info = qrdecode_version(ver_code_disturb)
    @test ver_info == ver_decode_info

    ## get version from QR-matrix
    ## invalid matrix size
    mat = rand(Bool, 21, 25)
    @test_throws DimensionMismatch qrdecode_version(mat)
    ## version information with too much errors
    mat = qrcode("HELLO WORLD"; version=7)
    mat[end-11:end-9, 1:6] .= 0
    @test_throws InfoError qrdecode_version(mat)
    mat = qrcode("HELLO WORLD"; version=7)
    mat[1:6, end-11:end-9] .= 0
    @test_throws InfoError qrdecode_version(mat)

    mat = qrcode("HELLO WORLD", Low(); compact=true)
    @test qrdecode_version(mat) == 1
    ## low degree
    tag = true
    for v in 1:6
        mat = qrcode("HELLO WORLD", Low(); compact=true, version=v)
        if qrdecode_version(mat;noerror=true) != v
            tag = false
            break
        end
    end
    @test tag
    ## higher degree
    tag = true
    for v in 7:40
        mat = qrcode("HELLO WORLD", High(); compact=true, version=v)
        if qrdecode_version(mat;noerror=true) != v
            tag = false
            break
        end
    end
    @test tag
end

@testset "Format information" begin
    ## format encoding
    tag = true
    for ((qua, mask), val) in formatinfo
        fmt = quality2binary[qua] << 3 ⊻ mask
        fmt_code = qrformat(fmt)
        val = parse(Int, join(Int.(val)); base=2)
        if qrdecode_format(fmt_code) != fmt || fmt_code != val
            tag = false
            break
        end
    end
    @test tag
    @test qrformat.(0:31) == qrformat.(0x00:0x1f)
    @test_throws InfoError qrformat(-1)
    @test_throws InfoError qrformat(32)
    ## distance of the Format-Code
    @test hamming_distance(qrformat.(0:31)) == 5

    ## fromat decoding -- with disturbance
    ### errors with dectection capacity(≤5)
    fmt_info = rand(0:7)
    fmt_code = qrformat(fmt_info)
    fmt_code_disturb = fmt_code
    errors = sample(0:14, 5; replace=false)
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

    ## get format from QR-matrix
    mat = qrcode("HELLO WORLD", Low(); compact=true)
    ec, mask = qrdecode_format(mat; noerror=true)
    @test ec == Low() && mask == 0

    mat = qrcode("HELLO WORLD", Medium(); compact=true)
    ec, mask = qrdecode_format(mat; noerror=true)
    @test ec == Medium() && mask == 0

    mat = qrcode("HELLO WORLD", High(); compact=true)
    ec, mask = qrdecode_format(mat; noerror=true)
    @test ec == High() && mask == 5

    mat = qrcode("HELLO WORLD", Quartile(); compact=true)
    ec, mask = qrdecode_format(mat; noerror=true)
    @test ec == Quartile() && mask == 3
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