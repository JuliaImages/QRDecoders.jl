# Tests for syndrome decoding

## Functions about syndrome polynomial
@testset "Syndrome polynomial" begin
    ## polynomial_eval(p::Poly, x::Int)
    ### take values at 0 and 1
    p = randpoly(rand(1:255))
    @test polynomial_eval(p, 1) == reduce(⊻, p.coeff)
    @test polynomial_eval(p, 0) == first(p.coeff)

    ### polynomial of degree ≤ 3
    a, b, c, d, x = rand(0:255, 5)
    @test polynomial_eval(Poly([a, b, c, d]), x) == a ⊻ mult(b ⊻ mult(c ⊻ mult(d, x), x), x)

    ## ---Dividing line--- ##
    ## syndrome_polynomial(message::Poly, n::Int)
    ### example from wikiversity: https://en.wikiversity.org/wiki/Reed%E2%80%93Solomon_codes_for_coders
    n = 10
    f = Poly(reverse!([0x40, 0xd2, 0x75, 0x47, 0x76, 0x17, 0x32, 0x06,
                0x27, 0x26, 0x96, 0xc6, 0xc6, 0x96, 0x70, 0xec]))
    msg = (f << n) + geterrorcorrection(f, n)
    msg.coeff[end] = 0 # deliberately damage the message
    syd = syndrome_polynomial(msg, n)
    @test syd.coeff == [64, 192, 93, 231, 52, 92, 228, 49, 83, 245]

    ## ---Dividing line--- ##
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

## Restoring erased information
@testset "erasure" begin
    ## evaluator_polynomial(errpos::AbstractVector)
    nerr = rand(1:255) # number of errors
    errpos = rand(0:255, nerr)
    errlocs = getindex.(Ref(logtable), 255 .- errpos) # error locator 2 .^ -(errpos)
    errpoly = erratalocator_polynomial(errpos)
    @test all(iszero, polynomial_eval.(Ref(errpoly), errlocs)) && length(errpoly) == nerr + 1

    ## evaluator_polynomial(syd::Poly, errloc::Poly, n::Int)
    
end