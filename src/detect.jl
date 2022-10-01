# detect QR code matrix(matrices) from an image

"""
    getalterpos(array::AbstractVector)

Get the switch-value indexes of the input `array`.
It can be used for detecting the finder patterns.

# Examples

```jl
julis> using QRDecoders: getalterpos
julia> array = Bool[1,1,0,0,0,1,1,0]
julia> getalterpos(array)
3-element Vector{Int64}:
 2
 5
 7
```
"""
getalterpos(array::AbstractVector) = findall(@views array[2:end] .!= array[1:end-1])

"""
    getqrmatrix(img)

Get the QR-code matrix from an image.

Note that the input must be a **standard** QR-code image
with or without white border, i.e. the image should not 
contain any other non-QR-Code information.
"""
function getqrmatrix(img::AbstractString)
    endswith(img, ".png") && return getqrmatrix(load(img))
    # binarize the image
    return getqrmatrix(round.(load(img)))
end

function getqrmatrix(img::AbstractMatrix)
    # find the left-top corner of the QR-code
    lt = findfirst(iszero, img)
    x1, y1 = lt.I
    # find the right down corner of the QR-code
    x2 = findlast(iszero, @view(img[:, y1]))
    y2 = y1 + x2 - x1
    rd = CartesianIndex((x2, y2))
    compactmat = @view img[lt:rd] # image crop
    # middle position of the left-top finder pattern
    mid = findfirst(!iszero, compactmat)[1] >> 1
    scale = findfirst(!iszero, @view(compactmat[mid, :])) - 1
    
    # rescale the QR-code matrix
    return .! Bool.(imresize(compactmat; ratio=1/scale))
end