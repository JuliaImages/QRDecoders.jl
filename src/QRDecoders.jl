module QRDecoders
using QRCoders
using QRCoders.Polynomial: Poly

export Poly

struct ReedSolomonError <: Exception
    st::AbstractString
    ReedSolomonError() = new("Number of errors exceeds the limitation of the RS-code")
end

include("qrinfo.jl")
include("syndrome.jl")



end