# decode message from the QRCode
## outline
### 1. de-interleave the message
### 2. correct and decode message
### 3. pack up the informations

using QRCoders: ecblockinfo, int2bitarray, getecblock, remainderbits, 
                charactercountlength, antialphanumeric, antikanji, bits2bytes
using .Syndrome: euclidean_decoder

## de-interleave the message
"""
    deinterleave(bytes::AbstractVector, ec::ErrCorrLevel, version::Int)

De-interleave the message, i.e. sperate msgblocks and ecblocks from data bits.
"""
function deinterleave(bytes::AbstractVector, ec::ErrCorrLevel, version::Int)
    ncodewords, nb1, nc1, nb2, nc2 = ecblockinfo[ec][version, :]
    return deinterleave(bytes, ncodewords, nb1, nc1, nb2, nc2)
end

function deinterleave(bytes::AbstractVector, ncodewords::Int, 
                      nb1::Int, nc1::Int, nb2::Int, nc2::Int)
    ## blocks of group1 and group2
    blocks = vcat([Vector{UInt8}(undef, nc1) for _ in 1:nb1],
                  [Vector{UInt8}(undef, nc2) for _ in 1:nb2])
    ## error correction blocks
    ecblocks = [Vector{UInt8}(undef, ncodewords) for _ in 1:(nb1 + nb2)]
    
    ind = length(bytes) # index start from the end of the message
    ## Error correction bytes
    for i in ncodewords:-1:1, j in (nb1 + nb2):-1:1
        ecblocks[j][i] = bytes[ind]
        ind -= 1
    end
    
    ## Encoded bytes
    for i in nc2:-1:(1 + nc1), j in (nb1 + nb2):-1:(nb1 + 1)
        blocks[j][i] = bytes[ind]
        ind -= 1
    end
    for i in nc1:-1:1, j in (nb1 + nb2):-1:1
        blocks[j][i] = bytes[ind]
        ind -= 1
    end
    # ind != 0 && throw(DecodeError("deinterleave: not all data is recorded"))
    return blocks, ecblocks
end

## correct and decode message
"""
    correct_message(msgblock::AbstractVector, ecblock::AbstractVector)

Correct the msgblock using Euclidean decoder.
"""
function correct_message(msgblock::AbstractVector, ecblock::AbstractVector)
    nsym = length(ecblock)
    # received polynomial = ecpoly + (msgpoly << nsym)
    receivedpoly = Poly(reverse!(vcat(msgblock, ecblock)))
    # error correction using Euclidean algorithm
    msgpoly = euclidean_decoder(receivedpoly, nsym)
    return reverse!(msgpoly.coeff)[1:length(msgblock)]
end

"""
    block2bits(blocks::AbstractVector)

Convert a block to bits.
"""
block2bits(block::AbstractVector) = vcat(int2bitarray.(block; pad=8)...)

"""
    decodemode(bits::AbstractVector)

Decode mode from the bits of length 4.
Note: the `Byte` mode and the `UTF8` mode use the same mode indicator(0100).
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
    decodedata(bits::AbstractVector, msglen::Int, ::Mode)

Decode message from `bits` without checking the `pad_bits`.
"""
function decodedata(bits::AbstractVector, msglen::Int, ::Numeric)
    nthrees, rem = msglen ÷ 3, msglen % 3
    ## main part of the integers
    bitlen1 = nthrees * 10
    ints1 = [bitarray2int(@view(bits[i:i+9])) for i in 1:10:bitlen1]
    str1 = join(string.(ints1; pad=3))
    ## remain part
    rem == 0 && return str1
    remnum = bitarray2int(@view(bits[bitlen1 + 1:bitlen1 + 1 + rem * 3]))
    return str1 * string(remnum;pad=rem)
end

function decodedata(bits::AbstractVector, msglen::Int, ::Alphanumeric)
    ntwos, rem = msglen >> 1, msglen & 1
    ## main part of the integers
    bitlen1 = ntwos * 11
    ints1 = [bitarray2int(@view(bits[i:i+10])) for i in 1:11:bitlen1]
    function int2str(k::Int)
        head, tail = k ÷ 45, k % 45
        return antialphanumeric[head] * antialphanumeric[tail]
    end
    str1 = join(int2str.(ints1))
    ## remain part of the integers
    rem == 0 && return str1
    remnum = bitarray2int(@view(bits[bitlen1 + 1:bitlen1 + 6]))
    return str1 * antialphanumeric[remnum]
end

function decodedata(bits::AbstractVector, msglen::Int, ::Kanji)
    ints = [bitarray2int(@view(bits[i:i+12])) for i in 1:13:msglen * 13]
    return join(getindex.(Ref(antikanji), ints))
end

function decodedata(bits::AbstractVector, msglen::Int, ::Byte)
    ints = bits2bytes(@view(bits[1:msglen << 3]))
    return join(Char.(ints))
end

function decodedata(bits::AbstractVector, msglen::Int, ::UTF8)
    bytes = bits2bytes(@view(bits[1:msglen << 3]))
    return String(bytes)
end

"""
    tryutf8(bits::AbstractVector, msglen::Int)

Try decoding message by utf-8 mode.
"""
function tryutf8(bits::AbstractVector, msglen::Int)
    # trybyte(bits, msglen) || return false
    bytes = bits2bytes(@view(bits[1:msglen << 3]))
    ## 1000 0000, 0100 0000, 0010 0000, 0001 0000, 0000 1000
    b0, b1, b2, b3, b4 = 0x80, 0x40, 0x20, 0x10, 0x08
    ind = 1
    ## try to read utf-8 characters
    while ind ≤ msglen
        num = bytes[ind]
        if num & b0 == 0 # ascii character(leading 0)
            ind += 1
        elseif num & b1 == 0 # first byte should not be 10xx xxxx
            return false
        elseif num & b2 == 0 # 110x xxxx
            ind + 1 > msglen && return false
            bytes[ind + 1] & 0xc0 == 0x80 || return false # should be 10xx xxxx
            ind += 2
        elseif num & b3 == 0 # 1110 xxxx
            ind + 2 > msglen && return false
            bytes[ind + 1] & 0xc0 == 0x80 || return false # should be 10xx xxxx
            bytes[ind + 2] & 0xc0 == 0x80 || return false # should be 10xx xxxx
            ind += 3
        elseif num & b4 == 0 # 1111 0xxx
            ind + 3 > msglen && return false
            bytes[ind + 1] & 0xc0 == 0x80 || return false # should be 10xx xxxx
            bytes[ind + 2] & 0xc0 == 0x80 || return false # should be 10xx xxxx
            bytes[ind + 3] & 0xc0 == 0x80 || return false # should be 10xx xxxx
            ind += 4
        else # first byte should not be 1111 10xx
            return false
        end
    end
    return true
end

"""
    trybyte(bits::AbstractVector, msglen::Int)

Try decoding message by `Byte` mode, where `bits` = `message_bits` + `pad_bits`.
Return `true` if the `pad_bits` is correct.
"""
function trybyte(bits::AbstractVector, msglen::Int)
    nbits = length(bits)
    nbits & 7 == 4 || throw(DecodeError("trybyte: incorrect bits length"))
    (nbits >> 3) < msglen && throw(DecodeError("trybyte: bits data is too short"))
    ## check the pad bits
    remainbytes = bits2bytes(@view(bits[msglen * 8 + 5:end]))
    nrem = length(remainbytes)
    return remainbytes == repeat([0xec, 0x11], ceil(Int, nrem / 2))[1:nrem]
end


## pack up the information
"""
    qrdecode(mat::AbstractMatrix; noerror::Bool=false, preferutf8=true)::QRInfo

QR code decoder.

If `noerror` is `true`, the decoder will raise an Exception(ReedSolomonError/InfoError)
when the QR code `mat` needs error correction.

If `preferutf8` is `true`, the decoder will try to decode the message by `UTF8` mode first
when dealing with `Byte` mode.
"""
function qrdecode(mat::AbstractMatrix; noerror::Bool=false, preferutf8=true)::QRInfo
    # --- decompose --- #
    ## 1. extract data bits from the QR-Matrix and 
    ## 2. check whether the matrix is valid or not
    version, eclevel, mask, msgbits = qrdecompose(mat;noerror=noerror)
    
    # --- get encoded bits --- #
    ## bits => bytes
    nbits = length(msgbits)
    nrb = remainderbits[version] # number of remainder bits
    (nbits - nrb) & 7 == 0 || throw(
        DecodeError("number of data bits(exclude remainders) should be a multiple of 8"))
    bytes = bits2bytes(msgbits)

    ## de-interleave: bytes ==> blocks & ecblocks
    msgblocks, ecblocks = deinterleave(bytes, eclevel, version)
    ncodewords = length(first(ecblocks))

    ## error correct: msgblocks and ecblocks ==> corrected message blocks
    crr_msgblocks = correct_message.(msgblocks, ecblocks)

    ## Requires that received bits do not contain errors
    if noerror
        crr_msgblocks != msgblocks && throw(DecodeError("Data-bits of the QR-Code contain errors"))
        decblocks = getecblock.(crr_msgblocks, ncodewords)
        decblocks != ecblocks && throw(DecodeError("Data-bits of the QR-Code contain errors"))
    end

    ## encoded = mode indicator + character count indicator + message bits
    ##         + pad bits(terminate + pad to multiple of 8 + repeated pattern)
    encoded = vcat(block2bits.(crr_msgblocks)...)

    ## extract format info of the encoded bits
    modeind = @view(encoded[1:4]) # mode indicator
    mode = decodemode(modeind) # Numeric, Alphanumeric, Byte, Kanji
    i = (version ≥ 1) + (version ≥ 10) + (version ≥ 27)
    cclength = charactercountlength[mode][i] # length of the char-count-indicator
    ccindicator = encoded[5:5 + cclength - 1] # character count indicator

    ## get message length from the c-c-indicator
    msglen = bitarray2int(ccindicator)

    # --- decode message --- #
    ## start position of the message bits
    startpos = 4 + cclength + 1
    messagebits = @view(encoded[startpos:end])

    ## decode message
    ### for Byte mode
    if mode == Byte()
        trybyte(messagebits, msglen) || throw(DecodeError("Failed to decode by Byte mode or UTF8 mode"))
        ## prefer to read UTF-8 mode
        if preferutf8 && tryutf8(messagebits, msglen)
            mode = UTF8()
        end
    end
    msg = decodedata(messagebits, msglen, mode)

    ## pack up the result
    return QRInfo(version, eclevel, mask, mode, msg)
end