testzbar = "testimages/zbar"
mkpath(testzbar)

@testset "Decode mode" begin
    # Numeric mode
    txt = join(0:9)
    exportqrcode(txt, "$testzbar/numeric.png")
    @test decodesingle("$testzbar/numeric.png") == 
          qrdecode("$testzbar/numeric.png").message ==
          decodeimg("$testzbar/numeric.png")[1] == txt

    # Alphanumeric mode
    txt = join(rand(keys(alphanumeric), 10))
    exportqrcode(txt, "$testzbar/alphanum.png")
    @test decodesingle("$testzbar/alphanum.png") ==
            qrdecode("$testzbar/alphanum.png").message ==
            decodeimg("$testzbar/alphanum.png")[1] == txt

    # Byte mode -- ASCII
    txt = join(Char.(0:127))
    exportqrcode(txt, "$testzbar/byte-ascii.png")
    @test decodesingle("$testzbar/byte-ascii.png") ==
            qrdecode("$testzbar/byte-ascii.png").message == txt
    ## decodeimg misdecode some message
    @test_broken decodeimg("$testzbar/byte-ascii.png")[1] == txt

    # UTF8 -- ZBar do not support UTF8
    txt = "你好"
    exportqrcode(txt, "$testzbar/utf8.png")
    @test qrdecode("$testzbar/utf8.png").message == txt
    @test_broken decodesingle("$testzbar/utf8.png") == txt
    @test_broken decodeimg("$testzbar/utf8.png")[1] == txt

    # Kanji mode
    txt = "茗荷"
    exportqrcode(txt, "$testzbar/kanji.png")
    @test decodesingle("$testzbar/kanji.png") ==
            qrdecode("$testzbar/kanji.png").message == txt
    ## decodeimg misdecode some message
    @test_broken decodeimg("$testzbar/kanji.png")[1] == txt
end