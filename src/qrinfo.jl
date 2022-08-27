# information of QR-code
## outline
### 1. Version of the QRCode
### 2. integrity of the QR-Code
### 3. error correction level and mask pattern
### 4. extract data bits

using QRCoders: finderpattern, alignmentpattern, alignmentlocation, makemasks, emptymatrix

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

"""
    hamming_distance(C::AbstractVector)

Hamming distance of a linear code.
"""
hamming_distance(C::AbstractVector) = minimum(hamming_weight, C)

### --- division line --- ###
## Version of the QR-Code

"""
    qrversion(fmt::Int)

Encode version information.
"""
function qrversion(ver::Int)
    7 ≤ ver ≤ 40 || throw(InfoError("version code $ver should be no less than 7 and no greater than 40"))
    # error correction code
    err = ver << 12
    # generator polynomial(= 0b1111100100101 in binary)
    g = Int(0x1f25) # use Int(0x1f25) to avoid overflow
    for i in 5:-1:0
        if !iszero(err & (1 << (i + 12)))
            err ⊻= g << i
        end
    end
    return ver << 12 ⊻ err
end
qrversion(ver::Integer) = qrversion(Int(ver))

"""
    qrdecode_version(version::Int)

Decode version information.
"""
function qrdecode_version(version_code::Int)
    ver_info, best_dist = -1, 18
    for test_info in 7:40
        test_code = qrversion(test_info)
        test_dist = hamming_distance(version_code, test_code)
        if test_dist < best_dist
            best_dist = test_dist
            ver_info = test_info
        elseif test_dist == best_dist
            ver_info = -1
        end
    end
    return ver_info
end

"""
    qrdecode_version(mat::AbstractMatrix; noerror)

Return the version of the QR-Code.
"""
function qrdecode_version(mat::AbstractMatrix; noerror=false)
    ## get version from size of the QR-matrix
    m, n = size(mat)
    m != n && throw(DimensionMismatch("matrix must be square"))
    (iszero((m - 17) % 4) && 21 ≤ m ≤ 177) || throw(InfoError("Invalid matrix size"))
    v = (m - 17) ÷ 4
    v < 7 && return v
    ## v ≥ 7 => 6x3 rectangular blocks that contain the version Information
    ## get version from the version string
    leftdown = mat[m-10:m-8, 1:6][:] # left-down block
    righttop = mat[1:6, m-10:m-8]'[:] # right-top block
    ldint = parse(Int, reverse(join(Int.(leftdown))); base=2) # bit string => integer
    rtint = parse(Int, reverse(join(Int.(righttop))); base=2) # bit string => integer
    v == qrdecode_version(ldint) || throw(InfoError("Version information(leftdown) not match"))
    v == qrdecode_version(rtint) || throw(InfoError("Version information(righttop) not match"))
    noerror && (ldint != rtint || qrversion(v) != ldint) && throw(InfoError("The QR-Code contains errors"))
    return v
end

### --- division line --- ###
## Format of the QRCode(EC + Mask)

"""
    qrformat(fmt::Int)

Generate standard format information (format + error correction + mask).
"""
function qrformat(fmt::Int)
    0 ≤ fmt ≤ 31 || throw(InfoError("format code $fmt should be no less than 0 and no greater than 31"))
    err = fmt << 10 # error correction code
    g = 0x537 # generator polynomial(= 0b10100110111 in binary)
    for i in 4:-1:0
        if !iszero(err & (1 << (i + 10)))
            err ⊻= g << i
        end
    end
    fmt << 10 ⊻ err ⊻ 0x5412 # mask(= 0b101010000010010 in binary)
end
qrformat(fmt::Integer) = qrformat(Int(fmt)) # to avoid integer overflow

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


"""
Bit modes of the quality.
"""
quality2binary = Dict(
    Low() => 0b01,
    Medium() => 0b00,
    Quartile() => 0b11,
    High() => 0b10)
binary2quality = Dict(val=>key for (key, val) in quality2binary)

"""
    qrdecode_format(mat::AbstractMatrix)

Return the format of the QR-Code(ErrCorrLevel + mask).
"""
function qrdecode_format(mat::AbstractMatrix; noerror=false)
    ## lefttop bar
    lefttop = vcat(mat[9, 1:6], mat[9, 8:9], mat[8,9], mat[6:-1:1, 9])
    ## down bar and right bar
    rightdown = vcat(mat[end:-1:end-6, 9], mat[9, end-7:end])
    ltint = foldl((i, j) -> (i << 1 ⊻ j), lefttop) # bit array => integer
    rdint = foldl((i, j) -> (i << 1 ⊻ j), rightdown) # bit array => integer
    ltfmt = qrdecode_format(ltint)
    rdfmt = qrdecode_format(rdint)
    ltfmt != rdfmt && throw(InfoError("Invalid format information"))
    noerror && (ltint != rdint || qrformat(ltfmt) != ltint) && throw(InfoError("The QR-Code contains errors"))
    ec, mask = ltfmt >> 3, ltfmt % 8
    return binary2quality[ec], mask
end

### --- division line --- ###
## Data Bits

"""
    extract_databits(mat::AbstractMatrix, datapos::AbstractMatrix)

Extract data bits from the QR-Code.(Inverse procedure of placedata!)
"""
function extract_databits(mat::AbstractMatrix, datapos::AbstractMatrix)
    n = size(mat, 1)
    data = BitArray{1}(count(datapos))
    col, row, cur = n, n + 1, 1
    while col > 0
        # Skip the column with the timing pattern
        if col == 7
            col -= 1
            continue
        end

        # path goes up and down
        if row > n
            row, δrow = n, -1
        else
            row, δrow = 1, 1
        end

        # extract data if the position is masked
        for _ in 1:n
            if datapos[row, col]
                data[cur] = mat[row, col]
                cur += 1
            end
            if datapos[row, col - 1]
                data[cur] = mat[row, col - 1]
                cur += 1
            end
            row += δrow
        end
        # go left
        col -= 2
    end
    return data
end

### --- division line --- ###
## decomposition of the QR-Code
"""
    qrdecompose(mat::AbstractMatrix)

Decompose the QR-Code into its constituent parts.
"""
qrdecompose(mat::AbstractMatrix; noerror=false) = qrdecompose(copy(mat); noerror=noerror)
function qrdecompose!(mat::AbstractMatrix; noerror=false)
    ## Version Information
    v = qrdecode_version(mat;noerror=noerror)
    n = v * 4 + 17

    ## Format Information
    ec, mask = qrdecode_format(mat;noerror=noerror)

    ## get mask matrix
    emptymat = emptymatrix(version) # place nothing on data-bits
    datapos = (emptymat .== nothing) # positions of the data-bits
    maskmat = makemasks(emptymat)[mask] # mask matrix that will be applied to data-bits
    @assert !any(maskmat[(!).(datapos)]) "mask-bit should be zero at non-data position"

    ## check Patterns: Finder Pattern, Alignment Pattern, Timing Pattern and Dark Mode
    ## set nothing on Version Information and Format Information
    if v >= 7
        emptymat[1:6, n - 10:n - 8] .= nothing
        emptymat[n - 10:n - 8, 1:6] .= nothing
    end
    emptymat[9, vcat(1:6, 8:9)] .= nothing # left-bar
    emptymat[vcat(8, 6:-1:1), 9] .= nothing # top-bar
    emptymat[n:-1:n-6, 9] .= nothing # bottom-bar
    emptymat[9, n-7:n] .= nothing # right-bar
    ## locations of the Patterns(Exclude Version Information and Format Information)
    patternpos = (emptymat .!= nothing)
    mat[patternpos] == emptymat[patternpos] || throw(InfoError())

    ## apply mask
    mat .⊻= maskmat
    
    ## extract data bits
    databits = extract_databits(mat, datapos)
    return v, ec, mask, databits
end

"""
    correct_messages(databits::AbstractMatrix, ec::ErrCorrLevel)

Correct the data bits.
"""
# function correct_message(databits::AbstractMatrix, ec::ErrCorrLevel)
#    
# end