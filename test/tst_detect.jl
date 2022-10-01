# test image dectection

@testset "finder pattern" begin
    # tests needed for later works
    @test getalterpos([1,1,0,0,0,1,1,0]) == [2, 5, 7]
    @test getalterpos([1, 1, 1]) == Int[]
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
        mat2 = getqrmatrix("testimages/test.png")
        mat3 = getqrmatrix("testimages/test.jpg")
        @test mat == mat2 == mat3
    end
end