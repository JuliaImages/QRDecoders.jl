using QRDecoders
using Test
using QRCoders
import QRCoders: formatinfo
import QRCoders.Polynomial: mult, geterrorcorrection
import QRDecoders: hamming_weight, hamming_distance, qrformat, qrdecode_format, quality2binary
                #    polynomial_eval, syndrome_polynomial, haserrors,
                #    erratalocator_polynomial, error_evaluator


"""
    randpoly(n::Int)

random polynomial of degree n.
"""
randpoly(n::Int) = Poly([rand(0:255, n-1)..., rand(1:255)])

include("tst_qrinfo.jl")
# include("tst_syndrome.jl")