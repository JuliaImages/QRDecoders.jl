module QRDecoders
using QRCoders

"""
An error occurs during error-correction.
"""
struct ReedSolomonError <: Exception
    st::AbstractString
    ReedSolomonError() = new("Number of errors exceeds the limitation of the RS-code")
end

"""
The QR-matrix is invalid.
"""
struct InfoError <: Exception
    st::AbstractString
    InfoError() = new("Invalid QR-Patterns")
    InfoError(st) = new(st)
end

include("qrinfo.jl")
include("syndrome.jl")

end