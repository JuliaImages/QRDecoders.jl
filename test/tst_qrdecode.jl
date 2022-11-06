@testset "QR-Code Decoding" begin
    alphabets = join.(['0':'9', keys(alphanumeric), keys(kanji), Char.(0:255), Char.(0:127)])
    for (mode, alphabet) in zip(modes, alphabets), eclevel in eclevels
        cap = last(characterscapacity[(eclevel, mode)])
        msg = join(rand(alphabet, rand(1:cap)))
        @test getmode(msg) ⊆ mode
        version = getversion(msg, mode, eclevel)
        mat = qrcode(msg, eclevel=eclevel)
        info = qrdecode(mat; noerror=true, preferutf8=false)
        @test info.message == msg && 
              info.version == version &&
              info.mode ⊆ mode &&
              info.eclevel == eclevel
    end
end

@testset "Decoding from pictures" begin
    # test for encoding modes
    alphabets = join.(['0':'9', keys(alphanumeric), keys(kanji), Char.(0:255)])
    for (mode, alphabet) in zip(modes, alphabets)
        cap = last(characterscapacity[(Medium(), mode)])
        msg = join(rand(alphabet, rand(1:cap)))
        mat = qrcode(msg)
        info = qrdecode(mat; noerror=true, preferutf8=false)
        exportqrcode(msg, "testimages/test.png")
        exportqrcode(msg, "testimages/test.jpg")
        info2 = qrdecode("testimages/test.png"; noerror=true, preferutf8=false)
        info3 = qrdecode("testimages/test.jpg"; noerror=true, preferutf8=false)
        @test info == info2 == info3
    end

    # test for versions
    msg = "Hello, world!"
    for v in 1:40
        mat = qrcode(msg, version=v)
        info = qrdecode(mat; preferutf8=false)
        exportqrcode(msg, "testimages/test.png", version=v)
        exportqrcode(msg, "testimages/test.jpg", version=v)
        mat2 = getqrmatrix("testimages/test.png")
        mat3 = getqrmatrix("testimages/test.jpg")
        info2 = qrdecode("testimages/test.png"; noerror=true, preferutf8=false)
        info3 = qrdecode("testimages/test.jpg"; noerror=true, preferutf8=false)
        @test mat == mat2 == mat2
        @test info == info2 == info3
    end
end

@testset "Coding & Decoding" begin
    # --- enoding --- #
    # message
    alphabets = ['0':'9', keys(alphanumeric), keys(kanji), Char.(0:255), Char.(0:127)]
    for (mode, alphabet) in zip(modes, alphabets), eclevel in eclevels
        cap = last(characterscapacity[(eclevel, mode)])
        msg = join(rand(alphabet, rand(1:cap)))
        @test getmode(msg) ⊆ mode
        version = getversion(msg, mode, eclevel)
        @test_throws EncodeError getversion(join(rand(alphabet, cap + 1)), mode, eclevel)

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
        addversion!(matrix, version)

        # Apply mask and add format information
        maskedmats = [addformat!(xor.(matrix, mat), i-1, eclevel) 
                    for (i, mat) in enumerate(masks)]
        scores = penalty.(maskedmats)
        mask = argmin(scores) - 1
        matrix = maskedmats[mask + 1]

        mat = qrcode(msg; mode=mode, eclevel=eclevel)
        @test mat == matrix

        # --- decoding --- #
        d_version, d_eclevel, d_mask, d_msgbits = qrdecompose(mat; noerror=true)

        # get msg bits
        @test d_msgbits == msgbits && d_eclevel == eclevel && d_mask == mask && d_version == version

        nbits = length(d_msgbits)
        nrb = remainderbits[d_version]
        @test (nbits - nrb) & 7 == 0 # number of bits is a multiple of 8
        bytes = bits2bytes(d_msgbits)

        # de-interleave blocks
        d_blocks, d_ecblocks = deinterleave(bytes, d_eclevel, d_version)
        @test d_blocks == blocks && d_ecblocks == ecblocks
        d_ncodewords = length(first(ecblocks))
        @test d_ncodewords == ncodewords


        # correct message using Euclidean algorithm (no error in this case)
        d_ec_blocks = correct_message.(d_blocks, d_ecblocks, Ref(Euclidean()))
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
        @test d_mode == mode || all(∈((UTF8(), Byte())), [d_mode, mode])
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
        d_encodeddata = d_encoded[startpos:end]
        if d_mode == Numeric()
            bitslen = (msglen ÷ 3) * 10 + (msglen % 3) * 3 + (msglen % 3 != 0)
        elseif d_mode == Alphanumeric()
            bitslen = (msglen ÷ 2) * 11 + (msglen % 2) * 6
        elseif d_mode == Kanji()
            bitslen = msglen * 13
        elseif d_mode == Byte()
            bitslen = msglen * 8
            # try byte
            @test trybyte(d_encodeddata, d_msglen)
        else # UTF8() for ascii characters
            bitslen = msglen * 8
            @test trybyte(d_encodeddata, d_msglen) && tryutf8(d_encodeddata, d_msglen)
        end
        @test @view(d_encodeddata[1:bitslen]) == encodeddata
        ## get data
        d_msg = decodedata(d_encodeddata, d_msglen, d_mode)
        @test d_msg == msg
    end
end

@testset "Message Decoding" begin
    ## test cases from https://www.thonky.com/qr-code-tutorial/
    ## Numeric
    msg = "8675309"
    bits = parse.(Bool, collect("110110001110000100101001"))
    @test decodedata(bits, length(msg), Numeric()) == msg

    ## Alphanumeric
    msg = "HELLO WORLD"
    txt = "0110000101101111000110100010111001011011100010011010100001101"
    bits = parse.(Bool, collect(txt))
    @test decodedata(bits, length(msg), Alphanumeric()) == msg

    ## Byte
    msg = "Hello, world!"
    txt = "010010000110010101101100011011000110111100101100001000000" *
    "11101110110111101110010011011000110010000100001"
    bits = parse.(Bool, collect(txt))
    @test decodedata(bits, length(msg), Byte()) == msg
    
    ## Kanji
    msg = "茗荷"
    bits = parse.(Bool, collect("11010101010100011010010111"))
    @test decodedata(bits, length(msg), Kanji()) == msg

    ## test case from qrcode standard ISO/IEC 18004 (2000)
    ## Numeric
    msg = "01234567"
    txt = "000000110001010110011000011"
    bits = parse.(Bool, collect(txt))
    @test decodedata(bits, length(msg), Numeric()) == msg
    msg = "0123456789012345"
    txt = "000000110001010110011010100110111000010100111010100101"
    bits = parse.(Bool, collect(txt))
    @test decodedata(bits, length(msg), Numeric()) == msg

    ## Alphanumeric
    msg = "AC-42"
    txt = "0011100111011100111001000010"
    bits = parse.(Bool, collect(txt))
    @test decodedata(bits, length(msg), Alphanumeric()) == msg

    # --- random test --- #
    ## Numeric
    msg = join(rand(0:9, rand(1:5596)))
    bits = encodedata(msg, Numeric())
    @test decodedata(bits, length(msg), Numeric()) == msg

    ## Alphanumeric
    msg = join(rand(keys(alphanumeric), rand(1:3391)))
    bits = encodedata(msg, Alphanumeric())
    @test decodedata(bits, length(msg), Alphanumeric()) == msg

    ## Byte
    msg = "Hello, world!"
    txt = "010010000110010101101100011011000110111100101100001000000" *
    "11101110110111101110010011011000110010000100001"
    bits = parse.(Bool, collect(txt))
    @test decodedata(bits, length(msg), Byte()) == msg
    @test decodedata(bits, length(msg), UTF8()) == msg

    ## Kanji
    msg = join(rand(keys(kanji), rand(1:1435)))
    bits = encodedata(msg, Kanji())
    @test decodedata(bits, length(msg), Kanji()) == msg
end

@testset "Byte Vs UTF8" begin
    msg = "Hello, world!"
    mat = qrcode(msg, eclevel=Low())
    info = qrdecode(mat; noerror=true, preferutf8=false)
    @test info.mode == Byte() && info.eclevel == Low() && info.message == msg
    info = qrdecode(mat; noerror=true)
    @test info.mode == UTF8() && info.eclevel == Low() && info.message == msg

    msg = "123αβ"
    mat = qrcode(msg, eclevel=Low())
    info = qrdecode(mat; noerror=true)
    @test info.mode == UTF8() && info.eclevel == Low() && info.message == msg

    msg = "123ø" # ø(0b1111 1000) is not a valid UTF8 byte
    bits = vcat(int2bitarray.(UInt8.(collect(msg)))...)
    @test !tryutf8(bits, 4)

    mat = qrcode(msg, eclevel=Low())
    info = qrdecode(mat; noerror=true)
    @test info.mode == Byte() && info.eclevel == Low() && info.message == msg
    
    msg = "®123" # ®(0b10101110) can not be the first byte in UTF8
    mat = qrcode(msg, eclevel=Low())
    info = qrdecode(mat; noerror=true)
    @test info.mode == Byte() && info.eclevel == Low() && info.message == msg
end

@testset "QRCode -- UTF8 mode" begin
    mode = UTF8()
    for v in 1:40, eclevel in eclevels
        cap = last(characterscapacity[(eclevel, mode)]) >> 2
        msg = join(rand(Char, rand(1:cap)))
        version = getversion(msg, mode, eclevel)
        mat = qrcode(msg, mode=mode, eclevel=eclevel)
        info = qrdecode(mat; noerror=true, preferutf8=true)
        @test info.message == msg && 
              info.version == version &&
              info.mode ⊆ mode &&
              info.eclevel == eclevel
    end
end

@testset "test throws" begin
    ## unsupport mode
    modeindicator = BitArray([0, 1, 1, 1]) ## ECI mode
    @test_throws DecodeError decodemode(modeindicator)

    ## invalid try
    @test_throws DecodeError trybyte(BitArray(undef, 16), 2) ## missing terminate bits
    @test_throws DecodeError trybyte(BitArray(undef, 16), 3) ## bits data is too short

    ## modified Data
    mat = qrcode("123")
    mat[1, 9] = !mat[1, 9] ## format info
    @test_throws InfoError qrdecode(mat; noerror=true)
    mat[1, 9] = !mat[1, 9]
    
    mat[1, 10] = !mat[1, 10]
    @test_throws DecodeError qrdecode(mat; noerror=true)
end