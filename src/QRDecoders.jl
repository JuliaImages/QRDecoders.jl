module QRDecoders

## import functions from QRCoder.jl
using QRCode
using QRCode.Polynomial
import QRCode: formatinfo
import QRCode.Polynomial: logtable, antilogtable, mult, generator, lead, tail!, init!

## add functions to QRCode
include("qrcoder.jl")

end