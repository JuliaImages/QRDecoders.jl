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
    @test ver_dist ≤ length(errors)

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
    for eclevel in eclevels, mask in 0:7
        mat = qrcode("HELLO WORLD", eclevel=eclevel, compact=true, mask=mask)
        ec, m = qrdecode_format(mat; noerror=true) 
        @test ec == eclevel && m == mask
    end
end

@testset "Decompose data" begin
    ## test -- HELLO WORLD
    ### QR-Code with different error correction level and masks
    msg = "HELLO WORLD"
    for eclevel in eclevels
        mask = rand(0:7)
        mat = qrcode(msg, eclevel=eclevel, compact=true, mask=mask)
        v, ec, m, databits = qrdecompose(mat; noerror=true)
        msgbits = encodemessage(msg, Alphanumeric(), eclevel, v)
        @test ec == eclevel && m == mask && databits == msgbits
    end

    ### QR-Code with different versions and masks
    msg = "HELLO WORLD"
    for version in 1:40
        mask = rand(0:7)
        mat = qrcode(msg; compact=true, version=version, mask=mask)
        v, ec, m, databits = qrdecompose(mat; noerror=true)
        msgbits = encodemessage(msg, Alphanumeric(), Medium(), version)
        @test databits == msgbits && v == version && ec == Medium() && m == mask
    end
        
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
    ## version information with too much errors
    mat = qrcode("HELLO WORLD"; version=7)
    mat[end-11:end-9, 1:6] .= 0
    @test_throws InfoError qrdecode_version(mat)

    mat = qrcode("HELLO WORLD"; version=7)
    mat[1:6, end-11:end-9] .= 0
    @test_throws InfoError qrdecode_version(mat)

    ## format information with too much errors
    mat = qrcode("HELLO WORLD")
    mat[9, vcat(1:6, 8:9)] .= 0 # left-bar
    mat[vcat(8, 6:-1:1), 9] .= 0 # top-bar
    mat[end:-1:end-6, 9] .= 0 # bottom-bar
    mat[9, end-7:end] .= 0 # right-bar
    @test_throws InfoError qrdecode_format(mat)
end