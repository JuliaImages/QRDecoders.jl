module QRDecoders
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
    InfoError <: Exception

The non-data part of QR-matrix contains error.

For example, Finder pattern, Alignment pattern, Timing pattern, 
Format information, Version information, matrix size and etc.
"""
struct InfoError <: Exception
    st::AbstractString
    InfoError() = new("Invalid QR-Patterns")
    InfoError(st) = new(st)
end

"""
    DecodeError <: Exception

The data part of QR-matrix contains error.
"""
struct DecodeError <: Exception
    st::AbstractString
    DecodeError(st) = new(st)
end

mutable struct QRInfo
    version::Int # version info(1 ≤ v ≤ 40)
    eclevel::ErrCorrLevel # error correction level(High, Low, Medium, Quartile)
    mask::Int # mask pattern(0-7)
    mode::Mode # encoding mode: Numeric, Alphanumeric, Byte, Kanji, UTF8
    message::AbstractString # decoded data
end

include("qrinfo.jl")
include("syndrome.jl")
include("qrdecode.jl")

end