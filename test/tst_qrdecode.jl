
@testset "Coding & Decoding -- Numeric" begin
    # --- enoding --- #
    # message
    alphabet = join('0':'9')
    msg = join(rand(alphabet, rand(1:5596)))

    # options
    eclevel = Medium()
    mode = getmode(msg)
    version = getversion(msg, mode, eclevel)

    # Mode indicator
    modeindicator = modeindicators[mode]

    # Character count: part of the encoded message
    msglen = length(msg)
    ccindicator = getcharactercountindicator(msglen, version, mode)

    # Encoded data: main part of the encoded message
    encodeddata = encodedata(msg, mode)

    # Getting parameters for the error correction
    ncodewords, nb1, nc1, nb2, nc2 = ecblockinfo[eclevel][version, :]
    requiredbits = 8 * (nb1 * nc1 + nb2 * nc2)

    # Pad encoded message before error correction
    encoded = vcat(modeindicator, ccindicator, encodeddata)
    encoded = padencodedmessage(encoded, requiredbits)

    # Getting error correction codes
    blocks = makeblocks(encoded, nb1, nc1, nb2, nc2)
    ecblocks = getecblock.(blocks, ncodewords)

    # Interleave code blocks
    msgbits = interleave(blocks, ecblocks, ncodewords, nb1, nc1, nb2, nc2, version)

    # Pick the best mask
    matrix = emptymatrix(version)
    masks = makemasks(matrix)
    matrix = placedata!(matrix, msgbits)
    candidates = map(enumerate(masks)) do (i, m)
        i - 1, xor.(matrix, m)
    end
    mask, matrix = first(sort(candidates, by = penalty ∘ last))
    matrix = addformat(matrix, mask, version, eclevel)

    mat = qrcode(msg; compact=true)
    @test mat == matrix

    # --- decoding --- #
    d_version, d_eclevel, d_mask, d_msgbits = qrdecompose(mat; noerror=true)

    # get msg bits
    @test d_msgbits == msgbits && d_eclevel == eclevel && d_mask == mask && d_version == version

    nbits = length(d_msgbits)
    nrb = remainderbits[d_version]
    @test (nbits - nrb) & 7 == 0 # number of bits is a multiple of 8
    nbytes = nbits >> 3
    bytes = bitarray2int.([@view(d_msgbits[(i - 1) * 8 + 1:i * 8]) for i in 1:nbytes])

    # de-interleave blocks
    d_blocks, d_ecblocks = deinterleave(bytes, d_eclevel, d_version)
    @test d_blocks == blocks && d_ecblocks == ecblocks
    d_ncodewords = length(first(ecblocks))
    @test d_ncodewords == ncodewords


    # correct message using Euclidean algorithm (no error in this case)
    d_ec_blocks = correct_message.(d_blocks, d_ecblocks)
    @test d_ec_blocks == blocks
    d_ec_ecblocks = getecblock.(d_ec_blocks, d_ncodewords)
    @test d_ec_ecblocks == ecblocks
    
    d_encoded = vcat(block2bits.(d_ec_blocks)...)
    @test d_encoded == encoded

    # split up message
    ## mode + count + data + zeros
    ## decode mode indicator
    d_modeindicator = d_encoded[1:4]
    @test d_modeindicator == modeindicator
    d_mode = decodemode(d_modeindicator)
    @test d_mode == mode
    ## decode character count indicator
    i = (version ≥ 1) + (version ≥ 10) + (version ≥ 27)
    cclength = charactercountlength[d_mode][i]
    d_ccindicator = d_encoded[5:5 + cclength - 1]
    @test d_ccindicator == ccindicator
    ## message length
    d_msglen = bitarray2int(d_ccindicator)
    @test d_msglen == msglen
    ## encoded data
    startpos = 4 + cclength + 1 # start position
    if d_mode == Numeric()
        bitslen = (msglen ÷ 3) * 10 + (msglen % 3) * 3 + (msglen % 3 != 0)
        d_encodeddata = d_encoded[startpos:startpos + bitslen - 1]
    end
    @test d_encodeddata == encodeddata
    ## get data
    d_msg = decodemessage(d_encodeddata, d_msglen, d_mode)
    @test d_msg == msg
end

@testset "QR-Code Decoding -- Numeric" begin
    alphabet = join('0':'9')
    
    msg = join(rand(alphabet, rand(1:5596)))
    eclevel = Medium()
    mode = getmode(msg)
    version = getversion(msg, mode, eclevel)

    mat = qrcode(msg, eclevel=eclevel, compact=true)
    info = qrdecode(mat; noerror=false)
    @test info.message == msg && info.version == version && info.mode == mode && info.eclevel == eclevel
end


@testset "Message Decoding" begin
    ## test cases from https://www.thonky.com/qr-code-tutorial/
    ## Numeric
    msg = "8675309"
    bits = parse.(Bool, collect("110110001110000100101001"))
    @test decodemessage(bits, length(msg), Numeric()) == msg

    ## Alphanumeric
    msg = "HELLO WORLD"
    txt = "0110000101101111000110100010111001011011100010011010100001101"
    bits = parse.(Bool, collect(txt))
    @test decodemessage(bits, length(msg), Alphanumeric()) == msg

    ## Byte

    ## Kanji

    # --- random test --- #
    ## Numeric
    msg = join(rand(0:9, rand(1:5596)))
    bits = encodedata(msg, Numeric())
    @test decodemessage(bits, length(msg), Numeric()) == msg

    ## Alphanumeric
    msg = join(rand(keys(alphanumeric), rand(1:3391)))
    bits = encodedata(msg, Alphanumeric())
    @test decodemessage(bits, length(msg), Alphanumeric()) == msg

    ## Byte

    ## Kanji
end