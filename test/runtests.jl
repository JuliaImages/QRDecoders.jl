using QRDecoders
using Test
using StatsBase: sample
using QRCoders
using QRCoders: formatinfo
using QRCoders.Polynomial: mult, geterrorcorrection, gfpow2, iszeropoly, gflog2, unit, euclidean_divide, divide
using QRDecoders: hamming_weight, hamming_distance, qrformat, qrdecode_format, quality2binary, ReedSolomonError,
                  extended_euclidean_divide, Sugiyama_euclidean_divide, euclidean_decoder
using QRDecoders.Syndrome: polynomial_eval, syndrome_polynomial, haserrors, fillerasures,
                            derivative_polynomial, erratalocator_polynomial, evaluator_polynomial,
                            findroots, reducebyHorner, getpositions, BMdecoder

"""
    randpoly(n::Int)

Random polynomial of degree n.
"""
randpoly(n::Int) = Poly([rand(0:255, n-1)..., rand(1:255)])
randpoly(range::AbstractVector{Int}) = randpoly(rand(range))

include("tst_qrinfo.jl")
include("tst_syndrome.jl")
include("tst_euclidean.jl")