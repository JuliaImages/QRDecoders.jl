# Tests for Euclidean decoder

@testset "Euclidean division" begin
    ## extended_euclidean_divide(r₁::Poly, r₂::Poly)
    ## Sugiyama_euclidean_divide(r₁::Poly, r₂::Poly, upperdeg::Int)
    r = randpoly(55) ## common factor of f(x) and g(x)
    f = r * randpoly(1:200)
    g = r * randpoly(1:200)
    a, b, common = extended_euclidean_divide(f, g)
    @test iszeropoly(a * f + b * g + common)
    @test iszeropoly(f % common) && iszeropoly(g % common)
    a, b, common = Sugiyama_euclidean_divide(f, g, 100)
    @test iszeropoly(a * f + b * g + common)

    f, g = randpoly(rand(1:255)), Poly([0, 0])
    a, b, common = extended_euclidean_divide(f, g)
    @test iszeropoly(a * f + b * g + common)
    @test iszeropoly(f % common) && iszeropoly(g % common)

    f, g = Poly([0, 0]), randpoly(rand(1:255))
    a, b, common = extended_euclidean_divide(f, g)
    @test iszeropoly(a * f + b * g + common)
    @test iszeropoly(f % common) && iszeropoly(g % common)

    f = randpoly(20:200)
    g = randpoly(20:200)
    a, b, common = Sugiyama_euclidean_divide(f, g, 20)
    @test iszeropoly(a * f + b * g + common)
end

@testset "Sugiyama's adaptation of ED algorithm" begin
    nothing ## TODO
end

@testset "Euclidean decoder -- without erasures" begin
    rawmsg = randpoly(155)
    nsym = 100
    errpos = unique!(rand(0:254, 50)) # error positions
    msg = rawmsg << nsym + geterrorcorrection(rawmsg, nsym)
    received = copy(msg)
    received.coeff[1 .+ errpos] .⊻= rand(1:255, length(errpos))
    @test euclidean_decoder(received, nsym) == msg

    ### no errors
    @test euclidean_decoder(msg, nsym) == msg

    ### odd number of syndromes
    rawmsg = randpoly(100)
    nsym = 155
    errpos = unique!(rand(0:254, 77)) # error positions
    msg = rawmsg << nsym + geterrorcorrection(rawmsg, nsym)
    received = copy(msg)
    received.coeff[1 .+ errpos] .⊻= rand(1:255, length(errpos))
    @test euclidean_decoder(received, nsym) == msg

    ### too much errors(detected)
    rawmsg = randpoly(100)
    nsym = 155
    errpos = sample(0:254, 78; replace=false) # error positions
    msg = rawmsg << nsym + geterrorcorrection(rawmsg, nsym)
    received = copy(msg)
    received.coeff[1 .+ errpos] .⊻= rand(1:255, length(errpos))
    # might be undetected (special property of Euclidean decoder)
    # @test_throws ReedSolomonError euclidean_decoder(received, nsym)
    
    ### too much errors(undetected!)
    ### [0] -encode> [0, 0, 0] -transfer> [2, 3, 0] -correct> [2, 3, 1] -decode> [1] 
    rawmsg = Poly([0])
    nsym = 2
    msg = rawmsg << nsym + geterrorcorrection(rawmsg, nsym)
    errpos = [0, 1]
    received = copy(msg)
    received.coeff[1 .+ errpos] .⊻= [2, 3]
    @test !iszeropoly(euclidean_decoder(received, nsym) + msg)
end

@testset "Euclidean decoder -- with erasures" begin
    ## euclidean_decoder(received::Poly, erasures::AbstractVector, nsym::Int)

    ### length of the received message is too long
    rawmsg = randpoly(100)
    nsym = 156
    msg = rawmsg << nsym + geterrorcorrection(rawmsg, nsym)
    @test_throws DomainError euclidean_decoder(msg, nsym)
    
    ### number of erasures exceeds the capacity of RS-Code
    rawmsg = randpoly(150)
    nsym = 50
    msg = rawmsg << nsym + geterrorcorrection(rawmsg, nsym)
    erasures = sample(0:199, 51; replace=false)
    @test_throws ReedSolomonError euclidean_decoder(msg, erasures, nsym)

    ### number of errors within the capacity of RS-Code
    ### 2 * v + ρ = 60 + 40 == nsym
    rawmsg = randpoly(155)
    nsym = 100
    errpos = sample(0:254, 70; replace=false) # error positions
    erasures = errpos[31:70]
    msg = rawmsg << nsym + geterrorcorrection(rawmsg, nsym)
    received = copy(msg)
    received.coeff[1 .+ errpos] .⊻= rand(1:255, length(errpos))
    @test euclidean_decoder(received, erasures, nsym) == msg
    ### Γx and Ωx by euclidean_decoder
    sydpoly = syndrome_polynomial(received, nsym)
    Γx = erratalocator_polynomial(erasures)
    xn = Poly(push!(zeros(Int, nsym), 1))
    upperdeg = (nsym + length(erasures)) ÷ 2 - 1
    Λx, _, Ωx = Sugiyama_euclidean_divide(sydpoly * Γx, xn, upperdeg)
    Λx *= Γx
    ### Λx and Ωx by direct computation
    errloc = erratalocator_polynomial(errpos)
    evlpoly = evaluator_polynomial(sydpoly, errloc, nsym)
    c0 = errloc ÷ Λx
    @test iszeropoly(c0 * Λx + errloc) && iszeropoly(c0 * Ωx + evlpoly)

    ### no errors
    received = copy(msg)
    @test euclidean_decoder(received, erasures, nsym) == msg
    ### Γx and Ωx by euclidean_decoder
    sydpoly = syndrome_polynomial(received, nsym)
    Γx = erratalocator_polynomial(erasures)
    xn = Poly(push!(zeros(Int, nsym), 1))
    upperdeg = (nsym + length(erasures)) ÷ 2 - 1
    Λx, _, Ωx = Sugiyama_euclidean_divide(sydpoly * Γx, xn, upperdeg)
    Λx *= Γx
    @test iszeropoly(Λx)

    ### errors within the capacity of the RS-Code
    ### 55 + 50 * 2 ≤ 155
    rawmsg = randpoly(100)
    nsym = 155
    errpos = sample(0:254, 105; replace=false) # error positions
    erasures = errpos[51:105]
    msg = rawmsg << nsym + geterrorcorrection(rawmsg, nsym)
    received = copy(msg)
    received.coeff[1 .+ errpos] .⊻= rand(1:255, length(errpos))
    @test euclidean_decoder(received, erasures, nsym) == msg
    ### Γx and Ωx by euclidean_decoder
    sydpoly = syndrome_polynomial(received, nsym)
    Γx = erratalocator_polynomial(erasures)
    xn = Poly(push!(zeros(Int, nsym), 1))
    upperdeg = (nsym + length(erasures)) ÷ 2 - 1
    Λx, _, Ωx = Sugiyama_euclidean_divide(sydpoly * Γx, xn, upperdeg)
    Λx *= Γx
    ### Λx and Ωx by direct computation
    errloc = erratalocator_polynomial(errpos)
    evlpoly = evaluator_polynomial(sydpoly, errloc, nsym)
    c0 = errloc ÷ Λx
    @test iszeropoly(c0 * Λx + errloc) && iszeropoly(c0 * Ωx + evlpoly)

    ### too much errors(detected)
    ### 2 * v + ρ = d + 1
    rawmsg = randpoly(100)
    v, ρ = rand(20:50), rand(20:55)
    nsym = 2 * v + ρ - 1
    errpos = sample(0:(nsym + 100 - 1), nsym + 1; replace=false) # error positions
    erasures = errpos[1:ρ]
    msg = rawmsg << nsym + geterrorcorrection(rawmsg, nsym)
    received = copy(msg)
    received.coeff[1 .+ errpos] .⊻= rand(1:255, length(errpos))
    # might be undetected (special property of Euclidean decoder)
    # @test_throws ReedSolomonError BMdecoder(received, erasures, nsym)

    ### samll case
    rawmsg = Poly([0])
    nsym = 3
    msg = rawmsg << nsym + geterrorcorrection(rawmsg, nsym)
    errpos = [0, 1]
    erasures = [0]
    received = copy(msg)
    received.coeff[1 .+ errpos] .⊻= rand(1:255, 2)
    @test euclidean_decoder(received, erasures, nsym) == msg
    ### Γx and Ωx by euclidean_decoder
    sydpoly = syndrome_polynomial(received, nsym)
    Γx = erratalocator_polynomial(erasures)
    xn = Poly(push!(zeros(Int, nsym), 1))
    upperdeg = (nsym + length(erasures)) ÷ 2 - 1
    Λx, _, Ωx = Sugiyama_euclidean_divide(sydpoly * Γx, xn, upperdeg)
    Λx *= Γx
    ### Λx and Ωx by direct computation
    errloc = erratalocator_polynomial(errpos)
    evlpoly = evaluator_polynomial(sydpoly, errloc, nsym)
    c0 = errloc ÷ Λx
    @test iszeropoly(c0 * Λx + errloc) && iszeropoly(c0 * Ωx + evlpoly)

    ### too much errors(undetected!)
    ### [0] -encode> [0, 0, 0, 0, 0] -transfer> [*, *, 54, 15, 0] -correct> [64, 120, 54, 15, 1] -decode> [1]
    ### 2 + 2 * 2 > 4
    rawmsg = Poly([0])
    nsym = 4
    msg = rawmsg << nsym + geterrorcorrection(rawmsg, nsym)
    errpos = [0, 1, 2, 3]
    erasures = [0, 1]
    received = copy(msg)
    received.coeff[1 .+ errpos] .⊻= [rand(1:255), rand(1:255), 54, 15]
    @test !iszeropoly(euclidean_decoder(received, erasures, nsym) + msg)
    
    ### wrong erausres might decrease the limitation of the RS-Code
    ### [0] -encode> [0, 0, 0, 0, 0] -transfer> [0, 0, 54, 15, 0] -correct> [0, 0, 0, 0, 0] -decode> [0]
    received = copy(msg)
    errpos = [2, 3]
    received.coeff[1 .+ errpos] .⊻= [54, 15] ## contains only two errors
    @test euclidean_decoder(received, [0, 1], nsym) == Poly([64, 120, 54, 15, 1])
end