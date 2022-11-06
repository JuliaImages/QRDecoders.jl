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
    getqrmatrix(imgpath::AbstractString)

Get the QR-code matrix from an image.

Note that the input must be a **standard** QR-code image
with or without white border, i.e. the image should not 
contain any non-QR-Code information.
"""
function getqrmatrix(imgpath::AbstractString)
    ext, mat = last(split(imgpath, '.')), load(imgpath)
    # png image
    ext == "png" && return getqrmatrix(mat)
    # binarize for jpg
    ext == "jpg" && return getqrmatrix(round.(mat))
    # gif image
    if ext == "gif"
        # static image
        ndims(mat) == 2 || throw(ArgumentError("The input image $imgpath is not a static image. Try `getqrmatrices` instead."))
        return getqrmatrix(Gray.(mat))
    end
    # unsupported image format
    throw(ArgumentError("Unsupported image format for $ext"))
end

"""
    getqrmatrix(mat::AbstractMatrix)

Get the standard QR-code matrix from an image-matrix.
"""
function getqrmatrix(mat::AbstractMatrix)
    # find the left-top corner of the QR-code
    lt = findfirst(iszero, mat)
    x1, y1 = lt.I
    # find the right down corner of the QR-code
    x2 = findlast(iszero, @view(mat[:, y1]))
    y2 = y1 + x2 - x1
    rd = CartesianIndex((x2, y2))
    compactmat = @view mat[lt:rd] # image crop
    # middle position of the left-top finder pattern
    mid = findfirst(!iszero, compactmat)[1] >> 1
    scale = findfirst(!iszero, @view(compactmat[mid, :])) - 1
    
    # rescale the QR-code matrix
    return .! Bool.(imresize(compactmat; ratio=1/scale))
end

"""
    getqrmatrices(imgpath::AbstractString)

Get the standard QR-code matrices from a gif file.
"""
function getqrmatrices(imgpath::AbstractString)
    ext = last(split(imgpath, '.'))
    ext == "gif" || throw(ArgumentError("The input image $imgpath should be a gif image."))
    mat = Gray.(load(imgpath))
    ndims(mat) == 3 || throw(ArgumentError("The input image $imgpath is not an animated image."))
    return @views [getqrmatrix(mat[:, :, i]) for i in 1:size(mat, 3)]
end