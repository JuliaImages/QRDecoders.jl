using Test
using ZBar
using ImageTransformations: imresize
using FileIO: load, save
using StatsBase: sample
# QRCoders
using QRCoders
using QRCoders:
    # encode data
    makeblocks, getecblock, interleave, 
    emptymatrix, makemask, makemasks, penalty,
    placedata!, addformat!, addversion!,
    # tables
    alphanumeric, kanji, qrversion, qrformat, 
    ecblockinfo, remainderbits,
    # encode
    getmode, characterscapacity, modeindicators, 
    getcharactercountindicator, charactercountlength,
    padencodedmessage, encodedata, encodemessage,
    # data convert
    bitarray2int, int2bitarray, bits2bytes, utf8len

using QRCoders.Polynomial:
    # operator for GF(256) integers 
    gfpow2, gflog2, gfinv, mult, divide,
    # operator for polynomials
    iszeropoly, degree, unit,
    geterrcode, euclidean_divide

# QRDecoders
using QRDecoders
using QRDecoders:
    # version and format
    hamming_weight, hamming_distance, 
    qrdecode_version, qrdecode_format,
    # decompose QR code 
    extract_databits, deinterleave, block2bits, 
    # decode data
    correct_message, decodemode, decodedata,
    # Byte mode and UTF8 mode
    trybyte, tryutf8,
    # image processing
    getalterpos, getqrmatrix
                  
using QRDecoders.Syndrome:
    # polynomial tools
    findroots, reducebyHorner, getpositions,
    # syndrome decoding
    polynomial_eval, syndrome_polynomial,
    derivative_polynomial, evaluator_polynomial,
    # applications
    haserrors, fillerasures, RSdecoder,
    Euclidean, BerlekampMassey,
    # BM decoder
    erratalocator_polynomial, berlekamp_massey_decoder,
    # Euclidean decoder
    extended_euclidean_divide, Sugiyama_euclidean_divide, euclidean_decoder

"""
    randpoly(n::Int)

Random polynomial of degree n.
"""
randpoly(::Type{T}, n::Int) where T = Poly{T}([rand(0:255, n-1)..., rand(1:255)])
randpoly(n::Int) = randpoly(UInt8, n)
randpoly(::Type{T}, range::AbstractVector{Int}) where T = randpoly(T, rand(range))
randpoly(range::AbstractVector{Int}) = randpoly(UInt8, range)
eclevels = [Low(), Medium(), Quartile(), High()]
modes = [Numeric(), Alphanumeric(), Kanji(), Byte(), UTF8()]

include("randerr.jl")

# interact with ZBar
include("tst_zbar.jl")

# decompose of the qr code 
include("tst_qrinfo.jl")

# Syndrome decoding
include("tst_syndrome.jl")

# euclidean decoder
include("tst_euclidean.jl")

# decode message from QR matrix
include("tst_qrdecode.jl")

# decode with random errors
include("tst_decodeerr.jl")

# image detecting
include("tst_detect.jl")