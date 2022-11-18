# Tests for syndrome decoding

## Functions for polynomials
@testset "Operators for polynomials" begin
    ## polynomial_eval(p::Poly, x::Int)
    ### take values at 0 and 1
    p = randpoly(Int, 1:255)
    @test polynomial_eval(p, 1) == reduce(⊻, p.coeff)
    @test polynomial_eval(p, 0) == first(p.coeff)
    ### polynomial of degree ≤ 3
    a, b, c, d, x = rand(0:255, 5)
    @test polynomial_eval(Poly([a, b, c, d]), x) == a ⊻ mult(b ⊻ mult(c ⊻ mult(d, x), x), x)

    ## derivative_polynomial(p::Poly)
    p = Poly([1, 2, 3])
    @test derivative_polynomial(p) == Poly([2, 0])
    @test derivative_polynomial(Poly([1])) == Poly([0])

    ## example from wikiversity: https://en.wikiversity.org/wiki/Reed%E2%80%93Solomon_codes_for_coders
    n = 10
    f = Poly(reverse!([0x40, 0xd2, 0x75, 0x47, 0x76, 0x17, 0x32, 0x06,
                0x27, 0x26, 0x96, 0xc6, 0xc6, 0x96, 0x70, 0xec]))
    msg = (f << n) + geterrcode(f, n)
    msg.coeff[end] = 0 # deliberately damage the message
    syd = syndrome_polynomial(msg, n) # syndrome polynomial
    @test syd.coeff == [64, 192, 93, 231, 52, 92, 228, 49, 83, 245]
    @test haserrors(msg, n)
    errlocpoly = erratalocator_polynomial(UInt8, [0]) # error locator polynomial
    evlpoly = evaluator_polynomial(syd, errlocpoly, n) # evaluator polynomial
    xn = Poly{UInt8}(push!(zeros(Int, n), 1)) ## xn = x^n
    @test iszeropoly(evlpoly + syd * errlocpoly % xn)

    ## erratalocator_polynomial
    nerr = rand(1:255) # number of errors
    errpos = rand(0:255, nerr)
    polyroots = Int.(gfpow2.(255 .- errpos)) # roots of error locator polynomial are 2 .^ -(errpos)
    locpoly = erratalocator_polynomial(Int, errpos)
    @test all(iszero, polynomial_eval.(Ref(locpoly), polyroots)) && length(locpoly) == nerr + 1
    @test erratalocator_polynomial(Int, Int[]) == Poly([1])

    ## reducebyHorner(p::Poly, a::Int)
    p, a = randpoly(Int, 2:255), rand(0:255)
    qx = reducebyHorner(p, a)
    qeval = popfirst!(qx.coeff)
    @test qeval == polynomial_eval(p, a)
    @test iszeropoly(qx * Poly([a, 1]) + p + Poly([qeval]))
    @test qx == p ÷ Poly([a, 1])
    @test reducebyHorner(Poly([a]), rand(0:255)) == Poly([a]) ## constant

    ## findroots(p::Poly)
    roots = sort!(unique!(rand(0:255, rand(1:255))))
    p = reduce(*, Poly([r, 1]) for r in roots)
    @test sort!(findroots(p)) == roots
    push!(roots, roots[1]) # has duplicate root
    pdup =  reduce(*, Poly([r, 1]) for r in roots)
    @test isempty(findroots(pdup))
    @test isempty(findroots(Poly([1,0,1]))) ## (x - 1)²

    ## getpositions(Λx::Poly)
    positions = sort!(unique!(rand(0:254, rand(1:255))))
    Λx = erratalocator_polynomial(Int, positions)
    @test sort!(getpositions(Λx)) == positions
    positions = collect(0:254)
    Λx = erratalocator_polynomial(Int, positions)
    @test sort!(getpositions(Λx)) == positions
end

## Detect errors
@testset "Errors dectecting" begin
    ## haserrors(msg::Poly, n::Int)
    ### RS-Code can detect up to n errors
    ### message without errors
    f, n = randpoly(Int, 1:127), rand(1:127)
    msg = f << n + geterrcode(f, n)
    @test !haserrors(msg, n)

    ### message with one error
    msg.coeff[rand(eachindex(msg.coeff))] ⊻= rand(1:255)
    @test haserrors(msg, n)

    ### message with ≤n errors
    n = rand(10:127)
    msg = f << n + geterrcode(f, n)
    errors = unique!(rand(eachindex(msg.coeff), n))
    msg.coeff[errors] = (⊻).(msg.coeff[errors], rand(1:255, length(errors)))
    @test haserrors(msg, n)
end

## fillerasures
@testset "Fill erasures" begin
    ## fillerasures(received::Poly, errpos::AbstractVector, n::Int)
    ### raw message -- HELLO WORLD
    rawmsg = Poly(reverse!([0x40, 0xd2, 0x75, 0x47, 0x76, 0x17, 0x32, 0x06,
    0x27, 0x26, 0x96, 0xc6, 0xc6, 0x96, 0x70, 0xec]))
    n = 10
    msg = rawmsg << n + geterrcode(rawmsg, n)
    ### received message
    received = copy(msg)
    errpos = UInt8[0, 10, 20]
    received.coeff[1 .+ errpos] = [6, 7, 8]
    @test fillerasures(received, errpos, n) == msg
    @test fillerasures(received, vcat(errpos, UInt8[1, 2, 3]), n) == msg

    ### random test
    fdeg, n = rand(1:200), 55
    rawmsg = randpoly(Int, fdeg)
    msg = rawmsg << n + geterrcode(rawmsg, n)
    ### received message
    errpos = sample(0:(n + fdeg - 1), 55; replace=false)
    received = copy(msg)
    received.coeff[1 .+ errpos] .⊻= rand(1:255, length(errpos))
    @test fillerasures(received, errpos, n) == msg
    
    ### error exceeds limitation
    n = 10
    errpos = collect(1:11)
    rawmsg = randpoly(Int, 1:20)
    msg = rawmsg << n + geterrcode(rawmsg, n)
    @test_throws ReedSolomonError fillerasures(msg, errpos, n)

    ### no errors
    n = 55
    rawmsg = randpoly(Int, 200)
    msg = rawmsg << n + geterrcode(rawmsg, n)
    errpos = unique!(rand(0:254, 55))
    @test fillerasures(msg, errpos, n) == msg
    @test fillerasures(msg, Int[], n) == msg
end

@testset "berlekamp_massey_decoder -- without erasures" begin
    ## erratalocator_polynomial(received::Poly, nsym::Int)
    ## berlekamp_massey_decoder(received::Poly, nsym::Int)

    ### number of errors within the capacity of RS-Code
    rawmsg = randpoly(Int, 155)
    nsym = 100
    errpos = unique!(rand(0:254, 50)) # error positions
    msg = rawmsg << nsym + geterrcode(rawmsg, nsym)
    received = copy(msg)
    received.coeff[1 .+ errpos] .⊻= rand(1:255, length(errpos))
    sydpoly = syndrome_polynomial(received, nsym)
    Λx = erratalocator_polynomial(sydpoly, nsym)
    errlocpoly = erratalocator_polynomial(Int, errpos)
    @test Λx == errlocpoly
    @test berlekamp_massey_decoder(received, nsym) == msg

    ### no errors
    sydpoly = syndrome_polynomial(msg, nsym)
    Λx = erratalocator_polynomial(sydpoly, nsym)
    @test Λx == unit(Poly{Int})
    @test berlekamp_massey_decoder(msg, nsym) == msg

    ### odd number of syndromes
    rawmsg = randpoly(Int, 100)
    nsym = 155
    errpos = unique!(rand(0:254, 77)) # error positions
    msg = rawmsg << nsym + geterrcode(rawmsg, nsym)
    received = copy(msg)
    received.coeff[1 .+ errpos] .⊻= rand(1:255, length(errpos))
    sydpoly = syndrome_polynomial(received, nsym)
    Λx = erratalocator_polynomial(sydpoly, nsym)
    errlocpoly = erratalocator_polynomial(Int, errpos)
    @test Λx == errlocpoly
    @test berlekamp_massey_decoder(received, nsym) == msg

    ### too much errors(detected)
    rawmsg = randpoly(Int, 100)
    nsym = 155
    errpos = sample(0:254, 78; replace=false) # error positions
    msg = rawmsg << nsym + geterrcode(rawmsg, nsym)
    received = copy(msg)
    received.coeff[1 .+ errpos] .⊻= rand(1:255, length(errpos))
    sydpoly = syndrome_polynomial(received, nsym)
    @test_throws ReedSolomonError erratalocator_polynomial(sydpoly, nsym; check=true)
    @test_throws ReedSolomonError berlekamp_massey_decoder(received, nsym)
    
    ### too much errors(undetected!)
    ### [0] -encode> [0, 0, 0] -transfer> [2, 3, 0] -correct> [2, 3, 1] -decode> [1] 
    rawmsg = Poly([0])
    nsym = 2
    msg = rawmsg << nsym + geterrcode(rawmsg, nsym)
    errpos = [0, 1]
    received = copy(msg)
    received.coeff[1 .+ errpos] .⊻= [2, 3]
    errlocpoly = erratalocator_polynomial(Int, errpos)
    sydpoly = syndrome_polynomial(received, nsym)
    Λx = erratalocator_polynomial(sydpoly, nsym; check=true)
    @test !iszeropoly(errlocpoly + Λx)
    @test Λx == erratalocator_polynomial(Int, [2])
    @test !iszeropoly(berlekamp_massey_decoder(received, nsym) + msg)

    # test for UInt8
    rawmsg = randpoly(UInt8, 155)
    nsym = 100
    errpos = unique!(rand(0:254, 50)) # error positions
    msg = rawmsg << nsym + geterrcode(rawmsg, nsym)
    received = copy(msg)
    received.coeff[1 .+ errpos] .⊻= rand(0x1:0xff, length(errpos))
    sydpoly = syndrome_polynomial(received, nsym)
    Λx = erratalocator_polynomial(sydpoly, nsym)
    errlocpoly = erratalocator_polynomial(UInt8, errpos)
    @test Λx == errlocpoly
    @test berlekamp_massey_decoder(received, nsym) == msg

    ### no errors
    sydpoly = syndrome_polynomial(msg, nsym)
    Λx = erratalocator_polynomial(sydpoly, nsym)
    @test Λx == unit(Poly{UInt8})
    @test berlekamp_massey_decoder(msg, nsym) == msg

    ### odd number of syndromes
    rawmsg = randpoly(UInt8, 100)
    nsym = 155
    errpos = unique!(rand(0:254, 77)) # error positions
    msg = rawmsg << nsym + geterrcode(rawmsg, nsym)
    received = copy(msg)
    received.coeff[1 .+ errpos] .⊻= rand(0x1:0xff, length(errpos))
    sydpoly = syndrome_polynomial(received, nsym)
    Λx = erratalocator_polynomial(sydpoly, nsym)
    errlocpoly = erratalocator_polynomial(UInt8, errpos)
    @test Λx == errlocpoly
    @test berlekamp_massey_decoder(received, nsym) == msg

    ### too much errors(detected)
    rawmsg = randpoly(UInt8, 100)
    nsym = 155
    errpos = sample(0x0:0xfe, 78; replace=false) # error positions
    msg = rawmsg << nsym + geterrcode(rawmsg, nsym)
    received = copy(msg)
    received.coeff[1 .+ errpos] .⊻= rand(0x1:0xff, length(errpos))
    sydpoly = syndrome_polynomial(received, nsym)
    @test_throws ReedSolomonError erratalocator_polynomial(sydpoly, nsym; check=true)
    @test_throws ReedSolomonError berlekamp_massey_decoder(received, nsym)
    
    ### too much errors(undetected!)
    ### [0] -encode> [0, 0, 0] -transfer> [2, 3, 0] -correct> [2, 3, 1] -decode> [1] 
    rawmsg = Poly([0x0])
    nsym = 2
    msg = rawmsg << nsym + geterrcode(rawmsg, nsym)
    errpos = [0, 1]
    received = copy(msg)
    received.coeff[1 .+ errpos] .⊻= [2, 3]
    errlocpoly = erratalocator_polynomial(UInt8, errpos)
    sydpoly = syndrome_polynomial(received, nsym)
    Λx = erratalocator_polynomial(sydpoly, nsym; check=true)
    @test !iszeropoly(errlocpoly + Λx)
    @test Λx == erratalocator_polynomial(UInt8, [2])
    @test !iszeropoly(berlekamp_massey_decoder(received, nsym) + msg)
end

@testset "berlekamp_massey_decoder -- with erasures" begin
    ## -- The modified Forney syndromes is still debugging -- ##

    ### length of the received message is too long
    rawmsg = randpoly(Int, 100)
    nsym = 156
    msg = rawmsg << nsym + geterrcode(rawmsg, nsym)
    @test_throws DomainError berlekamp_massey_decoder(msg, nsym)
    
    ### number of erasures exceeds the capacity of RS-Code
    rawmsg = randpoly(Int, 150)
    nsym = 50
    msg = rawmsg << nsym + geterrcode(rawmsg, nsym)
    erasures = sample(0:199, 51; replace=false)
    @test_throws ReedSolomonError berlekamp_massey_decoder(msg, erasures, nsym)
end