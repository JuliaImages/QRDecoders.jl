```@meta
CurrentModule = QRDecoders
```

# QRDecoders

QR Codes decoder with support of `Numeric` mode, `Alphanumeric` mode, `Kanj` mode, `Byte` mode and `UTF8` mode.

The decoding rules of QRDecoders.jl are compatible with [QRCoders.jl](http://github.com/JuliaImages/QRCoders.jl).

## Decoding message from QR codes

```@docs
QRInfo
qrdecode
qrdecompose
getqrmatrix
```

## Decoding procedures

Check the version and format information
```@docs
QRDecoders.qrdecode_version
QRDecoders.qrdecode_format
```

Extract message bits
```@docs
QRDecoders.extract_databits
QRDecoders.deinterleave
QRDecoders.block2bits
```

Error correction
```@docs
QRDecoders.correct_message
```

Decode message
```@docs
QRDecoders.decodemode
QRDecoders.decodedata
```

## Syndrome Decoding

Algorithm for decoding error correction codes.

```@docs
ReedSolomonAlgorithm
Euclidean
BerlekampMassey
RSdecoder
berlekamp_massey_decoder
euclidean_decoder
```

Tools for Syndrome Decoding.

```@docs
QRDecoders.Syndrome.syndrome_polynomial
QRDecoders.Syndrome.evaluator_polynomial
QRDecoders.Syndrome.erratalocator_polynomial
QRDecoders.Syndrome.haserrors
QRDecoders.Syndrome.fillerasures
QRDecoders.Syndrome.Sugiyama_euclidean_divide
```

Tools for polynomials.

```@docs
QRDecoders.Syndrome.polynomial_eval
QRDecoders.Syndrome.derivative_polynomial
QRDecoders.Syndrome.findroots
QRDecoders.Syndrome.reducebyHorner
QRDecoders.Syndrome.getpositions
QRDecoders.Syndrome.extended_euclidean_divide
```
## Error types

```@docs
InfoError
DecodeError
ReedSolomonError
```