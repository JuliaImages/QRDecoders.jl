# test image dectection

@testset "finder pattern" begin
    # tests needed for later works
    @test getalterpos([1,1,0,0,0,1,1,0]) == [2, 5, 7]
    @test getalterpos([1, 1, 1]) == Int[]
end

@testset "Decoding from pictures -- static image" begin
    alphabets = join.(['0':'9', keys(alphanumeric), keys(kanji), Char.(0:255)])
    eclevel = Medium()
    for (i, mode, alphabet) in zip(1:4, modes, alphabets)
        cap = last(characterscapacity[(eclevel, mode)])
        msg = join(rand(alphabet, rand(1:cap)))
        mat = qrcode(msg)
        info = qrdecode(mat; noerror=true, preferutf8=false)
        exportqrcode(msg, "testimages/test$i.png")
        exportqrcode(msg, "testimages/test$i.jpg")
        exportqrcode(msg, "testimages/test$i.gif")
        mat2 = getqrmatrix("testimages/test$i.png")
        mat3 = getqrmatrix("testimages/test$i.jpg")
        mat4 = getqrmatrix("testimages/test$i.gif")
        @test mat == mat2 == mat3 == mat4
    end
end

@testset "Decoding from pictures -- animated imags" begin
    # 5 frames
    alphabets = join.(['0':'9', keys(alphanumeric), keys(kanji), Char.(0:255)])
    eclevel = Medium()
    for (i, mode, alphabet) in zip(1:4, modes, alphabets)
        i, mode, alphabet = 1, Numeric(), '0':'9'
        cap = last(characterscapacity[(eclevel, mode)])
        msgs = [join(rand(alphabet, rand(1:cap))) for _ in 1:5]
        # get max version of codes
        codes = QRCode.(msgs, width=0)
        maxv = maximum(getproperty.(codes, :version))
        setproperty!.(codes, :version, maxv)
        exportqrcode(msgs, "testimages/test$i.gif")
        mats = getqrmatrices("testimages/test$i.gif")
        @test mats == qrcode.(codes)
        @test codes == qrdecode_animate(mats)
    end
end

@testset "Decoding from pictures -- test errors" begin
    ## unsupported types for `getqrmatrix` and `qrdecode`
    @test_throws ArgumentError getqrmatrix("test.webp")
    @test_throws ArgumentError qrdecode("test.webp")
    ## animated images for `getqrmatrix` and `qrdecode`
    exportqrcode(["hello", "julia"], "testimages/test.gif")
    @test_throws ArgumentError qrdecode("testimages/test.gif")
    @test_throws ArgumentError getqrmatrix("testimages/test.gif")

    ## unsupported types for `getqrmatrices` or `qrdecode_animate`
    @test_throws ArgumentError getqrmatrix("testimages/test.webp")
    @test_throws ArgumentError getqrmatrix("testimages/test.apng") # etc...

    ## static images for `getqrmatrices` or `qrdecode_animate`
    exportqrcode(["test"], "testimages/test-static.gif")
    @test_throws ArgumentError getqrmatrices("testimages/test-static.gif")
    @test_throws ArgumentError getqrmatrices("testimages/test.png") 
    @test_throws ArgumentError getqrmatrices("testimages/test.jpg")
    @test_throws ArgumentError qrdecode_animate("testimages/test.png")
    @test_throws ArgumentError qrdecode_animate("testimages/test.jpg")
    @test_throws ArgumentError qrdecode_animate("testimages/test-static.gif")
end