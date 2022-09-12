using QRDecoders
using Test
using StatsBase: sample
using QRCoders
using QRCoders: formatinfo, versioninfo, encodemessage, qrcode, emptymatrix, makemasks, getcharactercountindicator,
                encodedata, padencodedmessage, makeblocks, getecblock, interleave, placedata!, modeindicators, 
                ecblockinfo, penalty, addformat!, bitarray2int, int2bitarray, charactercountlength, getmode, 
                remainderbits, alphanumeric, kanji, characterscapacity, issubset, bits2bytes, utf8len
using QRCoders.Polynomial: mult, geterrorcorrection, gfpow2, iszeropoly, gflog2, unit, euclidean_divide, divide
using QRDecoders: hamming_weight, hamming_distance, qrversion, qrdecode_version, qrformat, qrdecode_format, mode2bin, 
                  ReedSolomonError, InfoError, extract_databits, qrdecompose, deinterleave, correct_message, 
                  block2bits, decodemode, decodedata, qrdecode, trybyte, tryutf8, DecodeError
                  
using QRDecoders.Syndrome: polynomial_eval, syndrome_polynomial, haserrors, fillerasures,
                            derivative_polynomial, erratalocator_polynomial, evaluator_polynomial,
                            findroots, reducebyHorner, getpositions, BMdecoder,
                            extended_euclidean_divide, Sugiyama_euclidean_divide, euclidean_decoder

"""
    randpoly(n::Int)

Random polynomial of degree n.
"""
randpoly(n::Int) = Poly([rand(0:255, n-1)..., rand(1:255)])
randpoly(range::AbstractVector{Int}) = randpoly(rand(range))

include("tst_qrinfo.jl")
include("tst_syndrome.jl")
include("tst_euclidean.jl")
include("tst_qrdecode.jl")