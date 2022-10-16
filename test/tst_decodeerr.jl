# Decoding of QRCode with errors

@testset "Error correction of polynomials" begin
    # random test
    rawmsg = randpoly(100)
    nsym = 100
    msgpoly = rawmsg << nsym + geterrcode(rawmsg, nsym)
    received = randerr!(copy(msgpoly), 50)
    decodebyeu = RSdecoder(received, nsym, Euclidean())
    decodebybm = RSdecoder(received, nsym, BerlekampMassey())
    @test decodebyeu == msgpoly == decodebybm

    # error exceeds limitation
    rawmsg = Poly([1])
    nsym = 3
    # msgpoly = Poly([8, 14, 7, 1])
    msgpoly = rawmsg << nsym + geterrcode(rawmsg, nsym)
    received = Poly([100, 14, 7, 70])
    # It is true that the received polynomial is closer to the message polynomial
    # encode by Poly([1]) though the number of errors exceeds the limitation
    decodebyeu = RSdecoder(received, nsym, Euclidean())
    @test decodebyeu == msgpoly
    @test_throws ReedSolomonError RSdecoder(received, nsym, BerlekampMassey())
end

@testset "Error correction of QR matrix" begin
    msg = "Hello, world!"
    ## qrcode with random error
    mat = qrcode_with_randerr(msg)
    info = qrdecode(mat;  preferutf8=false)
    infobybm = qrdecode(mat; alg=BerlekampMassey(),  preferutf8=false)
    # read qr matrix from image
    exportfrommatrix(mat, "testimages/randerr.png", width=0)
    @test getqrmatrix("testimages/randerr.png") == mat

    @test info == infobybm
    @test_throws DecodeError qrdecode(mat; noerror=true)
    @test_throws DecodeError qrdecode(mat; alg=BerlekampMassey(), noerror=true)
    @test info.message == msg && info.mode == Byte() && info.eclevel == Medium()

    # QRCode with one error for each block
    mat = qrcode_with_randerr(msg ^ 176, 1) # max version
    info = qrdecode(mat;  preferutf8=false)
    infobybm = qrdecode(mat; alg=BerlekampMassey(), preferutf8=false)
    exportfrommatrix(mat, "testimages/randerr.png")
    @test getqrmatrix("testimages/randerr.png") == mat
    @test info == infobybm
    @test_throws DecodeError qrdecode(mat; noerror=true)
    @test_throws DecodeError qrdecode(mat; alg=BerlekampMassey(), noerror=true)
    @test info.message == msg ^ 176 && info.mode == Byte() && info.eclevel == Medium()

    ## qrcode with max error
    for i in 1:5:176
        mat = qrcode_with_maxerr(msg^i)
        info = qrdecode(mat; preferutf8=false)
        infobybm = qrdecode(mat; alg=BerlekampMassey(), preferutf8=false)
        exportfrommatrix(mat, "testimages/randerr$i.png")
        @test getqrmatrix("testimages/randerr$i.png") == mat
        @test info == infobybm
        @test_throws DecodeError qrdecode(mat; noerror=true)
        @test_throws DecodeError qrdecode(mat; alg=BerlekampMassey(), noerror=true)
        @test info.message == msg ^ i && info.mode == Byte() && info.eclevel == Medium()
    end

    ## qrcode with too much errors
    ## It is possible that Euclidean decoder corrects message with number of errors
    ## exceeds the limitation of RS code.
    received = Poly([100, 14, 7, 70])
    block = [received.coeff[end]]
    ecblock = reverse!(received.coeff[1:end-1])
    @test [1] == @test_logs (:warn, ReedSolomonError()
                 ) correct_message(block, ecblock, Euclidean())
    msg = "Hello, world!"
    ncodewords, _ = getecinfo(msg)
    mat = qrcode_with_randerr(msg, ncodewords >> 1 + 1)
    try
        # Euclidean decoder may correct the message
        qrdecode(mat; alg=Euclidean())
    catch
        @test_throws ReedSolomonError qrdecode(mat; alg=Euclidean())
    end
    @test_throws ReedSolomonError qrdecode(mat; alg=BerlekampMassey())

    msg = "你好αβ" ^ 233
    ncodewords, _ = getecinfo(msg)
    mat = qrcode_with_randerr(msg, ncodewords >> 1 + 1)
    try
        # Euclidean decoder may correct the message
        qrdecode(mat; alg=Euclidean())
    catch
        @test_throws ReedSolomonError qrdecode(mat; alg=Euclidean())
    end
    @test_throws ReedSolomonError qrdecode(mat; alg=BerlekampMassey())
end