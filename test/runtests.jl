using Test
using StatsBase: sample
# QRCoders
using QRCoders
using QRCoders:
    # encode data
    makeblocks, getecblock, interleave, 
    emptymatrix, makemasks, penalty,
    placedata!, addformat!, addversion!,
    # tables
    alphanumeric, kanji, qrversion, qrformat, 
    ecblockinfo, remainderbits,
    # encode
    getmode, characterscapacity, modeindicators, 
    getcharactercountindicator, charactercountlength,
    padencodedmessage, encodedata, encodemessage,
    # data convert
    bitarray2int, int2bitarray, bits2bytes

using QRCoders.Polynomial:
    # operator for GF(256) integers 
    gfpow2, gflog2, gfinv, mult, divide,
    # operator for polynomials
    iszeropoly, degree, unit,
    geterrorcorrection, euclidean_divide

# QRDecoders
using QRDecoders
using QRDecoders:
    # version and format
    hamming_weight, hamming_distance, 
    qrversion, qrdecode_version, qrformat, qrdecode_format,
    # decompose QR code 
    extract_databits, deinterleave, block2bits, 
    # decode data
    correct_message, decodemode, decodedata,
    # Byte mode and UTF8 mode
    trybyte, tryutf8,
    # algorithm
    Euclidean, BerlekampMassey
                  
using QRDecoders.Syndrome:
    # polynomial tools
    findroots, reducebyHorner, getpositions,
    # syndrome decoding
    polynomial_eval, syndrome_polynomial,
    derivative_polynomial, evaluator_polynomial,
    # applications
    haserrors, fillerasures, RSdecoder,
    # BM decoder
    erratalocator_polynomial, berlekamp_massey_decoder,
    # Euclidean decoder
    extended_euclidean_divide, Sugiyama_euclidean_divide, euclidean_decoder

"""
    randpoly(n::Int)

Random polynomial of degree n.
"""
randpoly(n::Int) = Poly([rand(0:255, n-1)..., rand(1:255)])
randpoly(range::AbstractVector{Int}) = randpoly(rand(range))
eclevels = [Low(), Medium(), Quartile(), High()]

include("tst_qrinfo.jl")
include("tst_syndrome.jl")
include("tst_euclidean.jl")
include("tst_qrdecode.jl")