using QRDecoders
using Test
using QRCoders
import QRCoders: formatinfo
import QRCoders.Polynomial: mult, geterrorcorrection, gfpow2, iszeropoly, gflog2, unit
import QRDecoders: hamming_weight, hamming_distance, qrformat, qrdecode_format, quality2binary, ReedSolomonError
import QRDecoders.Syndrome: polynomial_eval, syndrome_polynomial, haserrors, fillearsed,
                            derivative_polynomial, erratalocator_polynomial, evaluator_polynomial,
                            findroots, reducebyHorner, getposition


"""
    randpoly(n::Int)

Random polynomial of degree n.
"""
randpoly(n::Int) = Poly([rand(0:255, n-1)..., rand(1:255)])
randpoly(range::AbstractVector{Int}) = randpoly(rand(range))

include("tst_qrinfo.jl")
include("tst_syndrome.jl")