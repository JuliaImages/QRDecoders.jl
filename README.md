# QRDecoders.jl

QR Codes decoder with support of `Numeric` mode, `Alphanumeric` mode, `Kanj` mode, `Byte` mode and `UTF8` mode.

The decoding rules of QRDecoders.jl follow the ISO/IEC 18004:2000 standard and are compatible with [QRCoders.jl](http://github.com/JuliaImages/QRCoders.jl).

# Installation
The package is not yet registered. One can install it with the following command:
```julia
using Pkg
# Pkg.add("QRDecoders") # not yet registered
Pkg.add("http://github.com/JuliaImages/QRDecoders.jl)
```

# General Usage
## Example 1
Decode a QR Code from a compact QR code matrix:
```julia
julia> # using Pkg; Pkg.add("QRCoders")
julia> using QRCoders, QRDecoders
julia> mat = qrcode("Hello World!"); # generate a QR code matrix
julia> info = qrdecode(mat)
QRInfo(1, Medium(), 0, UTF8(), "Hello World!")
```

The data type `QRInfo` contains information of the QR Code. 

| field | description | type
|---|---| ---|
| version | version of QR Code | `Int64` |
| eclevel | error correction level | `ErrCorrLevel` |
| mask | mask pattern | `Int64` |
| mode | mode of QR Code | `Mode` |
| message | decoded message | `AbstractString` |

The result `QRInfo(1, Medium(), 0, UTF8(), "Hello World!")` means that the QR Code is of version 1, with medium error correction level, mask pattern 1, and is encoded by UTF8 mode with input message `"Hello World!"`.

## Example 2
Decode a QR Code from an image without non-QRCode information:
```julia
julia> exportqrcode("Hello World!", "qrcode.png")
julia> qrdecode("qrcode.png")
```

For more complicate cases, it is not supported yet, but will be added in the future work.

## Options
There are some options for the decoder that can be set by the keyword arguments.

First and foremost, we implement two algorithms for error correction:
 - Berlekamp-Massey algorithm
 - Sugiyama's adaption of Euclidean algorithm

The default algorithm is `Euclidean algorithm`, and one can specify the algorithm by the keyword argument `alg`. For example:
```julia
julia> qrdecode(mat; alg = BerlekampMassey())
```

The encoding mode `Byte` and `UTF8` shared the same mode indicator `0100`. Therefore, when dealing with cases about `Byte` mode or `UTF8` mode, the decoder will try to decode the message with `UTF8` mode first. If the decoding fails, then it will use the `Byte` mode. However, one can set the keyword `preferutf8` to `false` to tell the decoder to skip the `UTF8` check.

Examples:
```julia
julia> mat = qrcode("Hello world!");
julia> qrdecode(mat) # use UTF8 mode by defauls
QRInfo(1, Medium(), 5, UTF8(), "Hello world!")

julia> qrdecode(mat; preferutf8=false) # skip UTF8 check
QRInfo(1, Medium(), 5, Byte(), "Hello world!")

julia> mat = qrcode("©®"); # fail to use UTF8 mode
julia> qrdecode(mat)
QRInfo(1, Medium(), 5, Byte(), "©®")
```

Note: One should know that the `UTF8` mode and `Byte` mode are the same for ascii characters. It only matters when the input message contains non-ascii characters.

Furthermore, there is one more option for the decoder, which is `noerror`. If `noerror` is set to be `true`, the QR matrix should be error-free. Otherwise, the decoder will raise `DecodeError` when it detects any error. This can be useful for tests of the Encode-Decode process interacting with QRCoders.jl.

# Acknowledgement
The QRDecoders.jl is created as part of the [OSPP'2022 project](https://summer-ospp.ac.cn/) guided by [Johnny Chen](http://github.com/johnnychen94).