# Test for functions in QRDecoders.jl/src/qrinfo.jl
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

    mat = qrcode("HELLO WORLD", eclevel = Low(), compact=true)
    @test qrdecode_version(mat) == 1
    ## low degree
    tag = true
    for v in 1:6
        mat = qrcode("HELLO WORLD", eclevel = Low(), compact=true, version=v)
        if qrdecode_version(mat;noerror=true) != v
            tag = false
            break
        end
    end
    @test tag
    ## higher degree
    tag = true
    for v in 7:40
        mat = qrcode("HELLO WORLD", eclevel = High(), compact=true, version=v)
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
        fmt = mode2bin[qua] << 3 ⊻ mask
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
    mat = qrcode("HELLO WORLD", eclevel = Low(), compact=true)
    ec, mask = qrdecode_format(mat; noerror=true)
    @test ec == Low() && mask == 0

    mat = qrcode("HELLO WORLD", eclevel = Medium(), compact=true)
    ec, mask = qrdecode_format(mat; noerror=true)
    @test ec == Medium() && mask == 0

    mat = qrcode("HELLO WORLD", eclevel = Quartile(), compact=true)
    ec, mask = qrdecode_format(mat; noerror=true)
    @test ec == Quartile() && mask == 3

    mat = qrcode("HELLO WORLD", eclevel = High(), compact=true)
    ec, mask = qrdecode_format(mat; noerror=true)
    @test ec == High() && mask == 5
end

@testset "Decompose data" begin
    ## test -- HELLO WORLD
    ### QR-Code with different error correction level
    msg = "HELLO WORLD"
    eclevel = Low()
    msgbits = encodemessage(msg, Alphanumeric(), eclevel, 1)
    mat = qrcode(msg, eclevel = eclevel, compact=true)
    v, ec, mask, databits = qrdecompose(mat; noerror=true)
    @test databits == msgbits && v == 1 && ec == eclevel && mask == 0

    eclevel = Medium()
    msgbits = encodemessage(msg, Alphanumeric(), eclevel, 1)
    mat = qrcode(msg, eclevel = eclevel, compact=true)
    v, ec, mask, databits = qrdecompose(mat; noerror=true)
    @test databits == msgbits && v == 1 && ec == eclevel && mask == 0

    eclevel = Quartile()
    msgbits = encodemessage(msg, Alphanumeric(), eclevel, 1)
    mat = qrcode(msg, eclevel = eclevel, compact=true)
    v, ec, mask, databits = qrdecompose(mat; noerror=true)
    @test databits == msgbits && v == 1 && ec == eclevel && mask == 3

    eclevel = High()
    msgbits = encodemessage(msg, Alphanumeric(), eclevel, 2)
    mat = qrcode(msg, eclevel = eclevel, compact=true)
    v, ec, mask, databits = qrdecompose(mat; noerror=true)
    @test databits == msgbits && v == 2 && ec == eclevel && mask == 5

    ### QR-Code with different version
    msg = "HELLO WORLD"
    tag = true
    for version in 1:40
        msgbits = encodemessage(msg, Alphanumeric(), Medium(), version)
        mat = qrcode(msg; compact=true, version=version)
        v, ec, _, databits = qrdecompose(mat; noerror=true)
        if !(databits == msgbits && v == version && ec == Medium())
            tag = false
            break
        end
    end
    @test tag

    ## random test -- need more tests
    
end

@testset "Integrity of the QRCode -- modify one bit" begin
    ## test matrix
    msg = "HELLO WORLD"
    mat = qrcode(msg; compact=true)
    n = size(mat, 1)

    ## Finder pattern + Separators
    ### left top block
    i, j = rand(1:8, 2)
    mat[i, j] = !mat[i, j]
    @test_throws InfoError qrdecompose(mat)
    mat[i, j] = !mat[i, j]

    ### right top block
    i, j = rand(1:8), rand(n-7:n)
    mat[i, j] = !mat[i, j]
    @test_throws InfoError qrdecompose(mat)
    mat[i, j] = !mat[i, j]

    ### left bottom block
    i, j = rand(n - 7: n), rand(1:8)
    mat[i, j] = !mat[i, j]
    @test_throws InfoError qrdecompose(mat)
    mat[i, j] = !mat[i, j]

    ## Time patterns
    ### horizontal
    i, j = 7, rand(10:n - 8)
    mat[i, j] = !mat[i, j]
    @test_throws InfoError qrdecompose(mat)
    mat[i, j] = !mat[i, j]

    ### vertical
    i, j = rand(10:n - 8), 7
    mat[i, j] = !mat[i, j]
    @test_throws InfoError qrdecompose(mat)
    mat[i, j] = !mat[i, j]

    ## Dark mode
    mat[n-7, 9] = !mat[n-7, 9]
    @test_throws InfoError qrdecompose(mat)
    mat[n-7, 9] = !mat[n-7, 9]

    ## Alignment pattern -- need some tests

end

@testset "Corrupted version/format information" begin
    ## version information with too much errors -- need more tests
    mat = qrcode("HELLO WORLD"; version=7)
    mat[end-11:end-9, 1:6] .= 0
    @test_throws InfoError qrdecode_version(mat)

    mat = qrcode("HELLO WORLD"; version=7)
    mat[1:6, end-11:end-9] .= 0
    @test_throws InfoError qrdecode_version(mat)

    ## format information with too much errors -- need some tests

end