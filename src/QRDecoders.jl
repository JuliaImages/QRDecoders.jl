module QRDecoders

export euclidean_decoder, ReedSolomonError
export Euclidean, BerlekampMassey, RSdecoder
export qrdecompose, InfoError
export qrdecode, DecodeError
export QRInfo

using QRCoders

"""
    ReedSolomonError <: Exception

An error occurs during error-correction.
"""
struct ReedSolomonError <: Exception
    st::AbstractString
    ReedSolomonError() = new("Number of errors exceeds the limitation of the RS-code")
end

"""
    ReedSolomonAlgorithm

An abstract type for error correction algorithm of Reed Solomon code.
"""
abstract type ReedSolomonAlgorithm end
"""
    Euclidean <: ReedSolomonAlgorithm

Euclidean algorithm for error correction.
"""
struct Euclidean <: ReedSolomonAlgorithm end
"""
    BerlekampMassey <: ReedSolomonAlgorithm

Berlekamp-Massey algorithm for error correction.
"""
struct BerlekampMassey <: ReedSolomonAlgorithm end

"""
    InfoError <: Exception

The non-data part of QR-matrix contains error.

For example, Finder pattern, Alignment pattern, Timing pattern, 
Format information, Version information, matrix size and etc.
"""
struct InfoError <: Exception
    st::AbstractString
end

"""
    DecodeError <: Exception

The data part of QR-matrix contains error.
"""
struct DecodeError <: Exception
    st::AbstractString
end

mutable struct QRInfo
    version::Int # version info(1 ≤ v ≤ 40)
    eclevel::ErrCorrLevel # error correction level(High, Low, Medium, Quartile)
    mask::Int # mask pattern(0-7)
    mode::Mode # encoding mode: Numeric, Alphanumeric, Byte, Kanji, UTF8
    message::AbstractString # decoded data
end

function Base.:(==)(info1::QRInfo, info2::QRInfo)
    all(getfield(info1, f) == getfield(info2, f) for f in fieldnames(QRInfo))
end


include("qrinfo.jl")
include("syndrome.jl")
include("qrdecode.jl")

using .Syndrome: euclidean_decoder

end