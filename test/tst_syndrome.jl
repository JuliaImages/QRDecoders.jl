# Tests for syndrome decoding

## Functions about syndrome polynomial
@testset "Syndrome polynomials" begin
    ## polynomial_eval(p::Poly, x::Int)
    ### take values at 0 and 1
    p = randpoly(rand(1:255))
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
    msg = (f << n) + geterrorcorrection(f, n)
    msg.coeff[end] = 0 # deliberately damage the message
    syd = syndrome_polynomial(msg, n) # syndrome polynomial
    @test syd.coeff == [64, 192, 93, 231, 52, 92, 228, 49, 83, 245]
    @test haserrors(msg, n)
    errloc = erratalocator_polynomial([0]) # error locator polynomial
    evlpoly = evaluator_polynomial(syd, errloc, n) # evaluator polynomial
    xn = Poly(push!(zeros(Int, n), 1)) ## xn = x^n
    @test iszeropoly(evlpoly + syd * errloc % xn)

    ## erratalocator_polynomial(errpos::AbstractVector)
    nerr = rand(1:255) # number of errors
    errpos = rand(0:255, nerr)
    polyroots = gfpow2.(255 .- errpos) # roots of error locator polynomial are 2 .^ -(errpos)
    locpoly = erratalocator_polynomial(errpos)
    @test all(iszero, polynomial_eval.(Ref(locpoly), polyroots)) && length(locpoly) == nerr + 1
    @test erratalocator_polynomial(Int[]) == Poly([1])
end

## Detect errors
@testset "Errors dectecting" begin
    ## haserrors(msg::Poly, n::Int)
    ### RS-Code can detect up to n errors
    ### message without errors
    f, n = randpoly(rand(1:127)), rand(1:127)
    msg = f << n + geterrorcorrection(f, n)
    @test !haserrors(msg, n)

    ### message with one error
    msg.coeff[rand(eachindex(msg.coeff))] ⊻= rand(1:255)
    @test haserrors(msg, n)

    ### message with ≤n errors
    n = rand(10:127)
    msg = f << n + geterrorcorrection(f, n)
    errors = unique!(rand(eachindex(msg.coeff), n))
    msg.coeff[errors] = (⊻).(msg.coeff[errors], rand(1:255, length(errors)))
    @test haserrors(msg, n)
end

@testset "Syndrome decoding" begin
    ## fillearsed(recieved::Poly, errpos::AbstractVector, n::Int)
    ### raw message -- HELLO WORLD
    rawmsg = Poly(reverse!([0x40, 0xd2, 0x75, 0x47, 0x76, 0x17, 0x32, 0x06,
    0x27, 0x26, 0x96, 0xc6, 0xc6, 0x96, 0x70, 0xec]))
    n = 10
    msg = rawmsg << n + geterrorcorrection(rawmsg, n)
    ### recieved message
    recieved = copy(msg)
    errpos = [0, 10, 20]
    recieved.coeff[1 .+ errpos] = [6, 7, 8]
    @test fillearsed(recieved, errpos, n) == msg

    ### original message -- random
    fdeg, n = rand(1:200), 55
    rawmsg = randpoly(fdeg)
    msg = rawmsg << n + geterrorcorrection(rawmsg, n)
    ### recieved message
    errpos = unique!(rand(1:fdeg, 55))
    recieved = copy(msg)
    recieved.coeff[1 .+ errpos] .⊻= rand(1:255, length(errpos))
    @test fillearsed(recieved, errpos, n) == msg
    
    ### error exceeds limitation
    fdeg, n = rand(1:200), 10
    errpos = collect(1:11)
    rawmsg = randpoly(fdeg)
    msg = rawmsg << n + geterrorcorrection(rawmsg, n)
    @test_throws ReedSolomonError fillearsed(msg, errpos, n)

    ## reducebyHorner(p::Poly, a::Int)
    p, a = randpoly(rand(1:255)), rand(0:255)
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

    ## getposition(Λx::Poly)
    positions = sort!(unique!(rand(0:254, rand(1:255))))
    Λx = erratalocator_polynomial(positions)
    @test sort!(getposition(Λx)) == positions
    positions = collect(0:254)
    Λx = erratalocator_polynomial(positions)
    @test sort!(getposition(Λx)) == positions
end

@testset "Berlekamp-Massey-algorithm" begin
    ## erratalocator_polynomial(recieved::Poly, nsym::Int)
    ### number of errors within the capacity of RS-Code
    rawmsg = randpoly(155)
    nsym = 100
    errpos = unique!(rand(0:254, 50)) # error positions
    msg = rawmsg << nsym + geterrorcorrection(rawmsg, nsym)
    recieved = copy(msg)
    recieved.coeff[1 .+ errpos] .⊻= rand(1:255, length(errpos))
    Λx = erratalocator_polynomial(recieved, nsym)
    errloc = erratalocator_polynomial(errpos)
    @test Λx == errloc

    ### no errors
    Λx = erratalocator_polynomial(msg, nsym)
    @test Λx == unit(Poly)

    ### odd number of syndromes
    rawmsg = randpoly(100)
    nsym = 155
    errpos = unique!(rand(0:254, 77)) # error positions
    msg = rawmsg << nsym + geterrorcorrection(rawmsg, nsym)
    recieved = copy(msg)
    recieved.coeff[1 .+ errpos] .⊻= rand(1:255, length(errpos))
    Λx = erratalocator_polynomial(recieved, nsym)
    errloc = erratalocator_polynomial(errpos)
    @test Λx == errloc

    ### too much errors(detected)
    rawmsg = randpoly(100)
    nsym = 155
    errpos = rand(0:176) .+ collect(0:78) # error positions
    msg = rawmsg << nsym + geterrorcorrection(rawmsg, nsym)
    recieved = copy(msg)
    recieved.coeff[1 .+ errpos] .⊻= rand(1:255, length(errpos))
    @test_throws ReedSolomonError erratalocator_polynomial(recieved, nsym)
    
    ### too much errors(undetected!)
    ### [0] -encode> [0, 0, 0] -transfer> [2, 3, 0] -correct> [2, 3, 1] -decode> [1] 
    rawmsg = Poly([0])
    nsym = 2
    msg = rawmsg << nsym + geterrorcorrection(rawmsg, nsym)
    errpos = [0, 1]
    recieved = copy(msg)
    recieved.coeff[1 .+ errpos] .⊻= [2, 3]
    errloc = erratalocator_polynomial(errpos)
    Λx = erratalocator_polynomial(recieved, nsym)
    @test !iszeropoly(errloc + Λx)
    @test Λx == erratalocator_polynomial([2])
end