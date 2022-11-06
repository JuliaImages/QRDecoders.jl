module QRDecoders

export Euclidean, BerlekampMassey
export euclidean_decoder, berlekamp_massey_decoder, RSdecoder
export InfoError, DecodeError, ReedSolomonError
export qrdecompose, qrdecode, qrdecodegif, getqrmatrix, getqrmatrices
export QRInfo

using QRCoders
using ImageTransformations: imresize
using FileIO: load
using ColorTypes

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

"""
    QRInfo

A struct to store the information of a QR code.

## Fields
- `version::Int`: version of the QR code
- `eclevel::Int`: error correction level of the QR code
- `mask::Int`: mask of the QR code
- `mode::Int`: mode of the QR code
- `message::String`: decoded message
"""
mutable struct QRInfo
    version::Int # version info(1 ≤ v ≤ 40)
    eclevel::ErrCorrLevel # error correction level(High, Low, Medium, Quartile)
    mask::Int # mask pattern(0-7)
    mode::Mode # encoding mode: Numeric, Alphanumeric, Byte, Kanji, UTF8
    message::AbstractString # decoded data
end

"""
    ==(info1::QRInfo, info2::QRInfo)

Compare two `QRInfo` objects.
"""
function Base.:(==)(info1::QRInfo, info2::QRInfo)
    all(getfield(info1, f) == getfield(info2, f) for f in fieldnames(QRInfo))
end

"""
    ==(info, code)

Compare a `QRInfo` object and a `QRCode` object.
"""
function Base.:(==)(info::QRInfo, code::QRCode)
    all(field -> getfield(info, field) == getfield(code, field), fieldnames(QRInfo))
end
Base.:(==)(code::QRCode, info::QRInfo) = info == code

include("qrinfo.jl")
include("syndrome.jl")
include("qrdecode.jl")
include("detect.jl")

using .Syndrome: euclidean_decoder, berlekamp_massey_decoder, RSdecoder

# The other method of qrdecode can be found in qrdecode.jl
"""
    qrdecode(path::AbstractString; keywords...)

QR code decoder.

For more information of the keywords, see `qrdecode(mat::AbstractMatrix; keywords...)`.
"""
qrdecode(path::AbstractString; keywords...) = qrdecode(getqrmatrix(path); keywords...)

"""
    qrdecodes(path::AbstractString; keywords...)

QR code decoder for animated QR code.
"""
function qrdecodegif(path::AbstractString; keywords...)::Vector{QRInfo}
    mats = getqrmatrices(path)
    return qrdecodegif(mats; keywords...)
end

end