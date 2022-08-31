# decode message from the QRCode
## outline
### 1. de-interleave the message
### 2. correct message
### 3. decode message
### 4. pack up the informations

using QRCoders: ecblockinfo, int2bitarray, getecblock, remainderbits, 
                charactercountlength, antialphanumeric, antikanji
using .Syndrome: euclidean_decoder

"""
    deinterleave(bytes::AbstractVector, ec::ErrCorrLevel, version::Int)

De-interleave the message. Sperate msgblocks and ecblocks from data bits.
"""
function deinterleave(bytes::AbstractVector, ec::ErrCorrLevel, version::Int)
    ncodewords, nb1, nc1, nb2, nc2 = ecblockinfo[ec][version, :]
    return deinterleave(bytes, ncodewords, nb1, nc1, nb2, nc2)
end

function deinterleave(bytes::AbstractVector, ncodewords::Int, nb1::Int, nc1::Int, nb2::Int, nc2::Int)
    blocks = vcat([Vector{Int}(undef, nc1) for _ in 1:nb1], [Vector{Int}(undef, nc2) for _ in 1:nb2])
    ecblocks = [Vector{Int}(undef, ncodewords) for _ in 1:(nb1 + nb2)]
    ind = length(bytes) # index start from the end of the message

    ## Error correction block
    for i in ncodewords:-1:1, j in (nb1 + nb2):-1:1
        ecblocks[j][i] = bytes[ind]
        ind -= 1
    end
    
    ## Encoded block
    for i in nc2:-1:(1 + nc1), j in (nb1 + nb2):-1:(nb1 + 1)
        blocks[j][i] = bytes[ind]
        ind -= 1
    end
    for i in nc1:-1:1, j in (nb1 + nb2):-1:1
        blocks[j][i] = bytes[ind]
        ind -= 1
    end
    ind != 0 && throw(DecodeError("deinterleave: not all data is recorded"))
    return blocks, ecblocks
end

"""
    correct_message(msgblock::AbstractVector, ecblock::AbstractVector)

Correct the data bits.
"""
function correct_message(msgblock::AbstractVector, ecblock::AbstractVector)
    nsym = length(ecblock)
    # received polynomial = ecpoly + (msgpoly << nsym)
    recievedpoly = Poly(reverse!(vcat(msgblock, ecblock)))
    # error correction using Euclidean algorithm
    msgpoly = euclidean_decoder(recievedpoly, nsym)
    return reverse!(msgpoly.coeff)[1:length(msgblock)]
end

"""
    block2bits(blocks::AbstractVector)

Convert a block to bits.
"""
block2bits(block::AbstractVector) = vcat(int2bitarray.(block; pad=8)...)

"""
    decodemode(bits::AbstractVector)
"""
function decodemode(bits::AbstractVector)
    length(bits) != 4 && throw(DecodeError("decodemode: bits should be 4-bit long"))
    num = bitarray2int(bits)
    num == 1 && return Numeric()
    num == 2 && return Alphanumeric()
    num == 4 && return Byte() # or UTF8()
    num == 8 && return Kanji()
    throw(DecodeError("decodemode: unknown mode"))
end

"""
    decodemessage(bits::AbstractVector, msglen::Int, ::Mode)

Decode the message from bit-data.
"""
function decodemessage(bits::AbstractVector, msglen::Int, ::Numeric)
    nthrees, rem = msglen ÷ 3, msglen % 3
    ## main part of the integers
    bitlen1 = nthrees * 10
    ints1 = [bitarray2int(@view(bits[i:i+9])) for i in 1:10:bitlen1]
    str1 = join(string.(ints1; pad=3))
    ## remain part
    rem == 0 && return str1
    return str1 * string(bitarray2int(@view(bits[bitlen1 + 1:bitlen1 + 1 + rem * 3]));pad=rem)
end

function decodemessage(bits::AbstractVector, msglen::Int, ::Alphanumeric)
    ntwos, rem = msglen >> 1, msglen & 1
    ## main part of the integers
    bitlen1 = ntwos * 11
    ints1 = [bitarray2int(@view(bits[i:i+10])) for i in 1:11:bitlen1]
    function int2str(k::Int)
        head, tail = k ÷ 45, k % 45
        return antialphanumeric[head] * antialphanumeric[tail]
    end
    str1 = join(int2str.(ints1))
    ## remain part
    rem == 0 && return str1
    remnum = bitarray2int(@view(bits[bitlen1 + 1:bitlen1 + 6]))
    println(remnum)
    return str1 * antialphanumeric[remnum]
end

function decodemessage(bits::AbstractVector, msglen::Int, ::Kanji)
    ints = [bitarray2int(@view(bits[i:i+12])) for i in 1:13:msglen * 13]
    return join(getindex.(Ref(antikanji), ints))
end

function qrdecode(mat::AbstractMatrix; noerror=false)
    ## extract data bits from the QR-Matrix and 
    ## check whether the matrix is valid or not
    version, eclevel, mask, msgbits = qrdecompose(mat;noerror=noerror)

    ## bits => bytes
    nbits = length(msgbits)
    nrb = remainderbits[version] # number of remainder bits
    (nbits - nrb) & 7 == 0 || throw(
        DecodeError("number of data bits(exclude remainders) should be a multiple of 8"))
    nbytes = nbits >> 3
    bytes = bitarray2int.([@view(msgbits[(i - 1) * 8 + 1:i * 8]) for i in 1:nbytes])

    ## bits ==de-interleave=> blocks & ecblocks
    msgblocks, ecblocks = deinterleave(bytes, eclevel, version)
    ncodewords = length(first(ecblocks))

    ## msgblocks & ecblocks ==error correct=> corrected message blocks
    dmsgblocks = correct_message.(msgblocks, ecblocks)
    ## require the code to contain no errors
    if noerror
        dmsgblocks != msgblocks && throw(DecodeError("Data-bits of the QR-Code contain errors"))
        decblocks = getecblock.(dmsgblocks, ncodewords)
        decblocks != ecblocks && throw(DecodeError("Data-bits of the QR-Code contain errors"))
    end
    
    ## encoded bits
    encoded = vcat(block2bits.(dmsgblocks)...)

    ## extract format info of the encoded bits
    ## message = mode indicator + character count indicator + encoded data
    modeind = @view(encoded[1:4])
    mode = decodemode(modeind) # mode: Numeric, Alphanumeric, Byte, Kanji
    i = (version ≥ 1) + (version ≥ 10) + (version ≥ 27)
    cclength = charactercountlength[mode][i] # length of the char-count-indicator
    ccindicator = encoded[5:5 + cclength - 1] # character count indicator
    ## get message length from the c-c-indicator
    msglen = bitarray2int(ccindicator)
    ## start position of the encoded bits
    startpos = 4 + cclength + 1

    ## decode message
    msg = decodemessage(@view(encoded[startpos:end]), msglen, mode)

    return QRInfo(version, eclevel, mask, mode, msg)
end