var documenterSearchIndex = {"docs":
[{"location":"","page":"QRDecoders","title":"QRDecoders","text":"CurrentModule = QRDecoders","category":"page"},{"location":"#QRDecoders","page":"QRDecoders","title":"QRDecoders","text":"","category":"section"},{"location":"","page":"QRDecoders","title":"QRDecoders","text":"QR Codes decoder with support of Numeric mode, Alphanumeric mode, Kanj mode, Byte mode and UTF8 mode.","category":"page"},{"location":"","page":"QRDecoders","title":"QRDecoders","text":"The decoding rules of QRDecoders.jl are compatible with QRCoders.jl.","category":"page"},{"location":"#Decoding-message-from-QR-codes","page":"QRDecoders","title":"Decoding message from QR codes","text":"","category":"section"},{"location":"","page":"QRDecoders","title":"QRDecoders","text":"QRInfo\nqrdecode\nqrdecodegif\nqrdecompose\ngetqrmatrix\ngetqrmatrices","category":"page"},{"location":"#QRDecoders.QRInfo","page":"QRDecoders","title":"QRDecoders.QRInfo","text":"QRInfo\n\nA struct to store the information of a QR code.\n\nFields\n\nversion::Int: version of the QR code\neclevel::Int: error correction level of the QR code\nmask::Int: mask of the QR code\nmode::Int: mode of the QR code\nmessage::String: decoded message\n\n\n\n\n\n","category":"type"},{"location":"#QRDecoders.qrdecode","page":"QRDecoders","title":"QRDecoders.qrdecode","text":"qrdecode(mat::AbstractMatrix\n        ; noerror::Bool=false\n        , preferutf8::Bool=true\n        , alg::ReedSolomonAlgorithm=Euclidean()\n        )::QRInfo\n\nQR code decoder.\n\nIf noerror is true, the decoder will raise an Exception(ReedSolomonError/InfoError) when the QR code mat needs error correction.\n\nIf preferutf8 is true, the decoder will try to decode the message by UTF8 mode when dealing with Byte mode.\n\nThe error correction algorithm is specified by the variable alg(default: Euclidean).\n\n\n\n\n\nqrdecode(path::AbstractString; keywords...)\n\nQR code decoder.\n\nFor more information of the keywords, see qrdecode(mat::AbstractMatrix; keywords...).\n\n\n\n\n\n","category":"function"},{"location":"#QRDecoders.qrdecodegif","page":"QRDecoders","title":"QRDecoders.qrdecodegif","text":"qrdecodes(path::AbstractString; keywords...)\n\nQR code decoder for animated QR code.\n\n\n\n\n\n","category":"function"},{"location":"#QRDecoders.qrdecompose","page":"QRDecoders","title":"QRDecoders.qrdecompose","text":"qrdecompose(mat::AbstractMatrix, noerror=false)\n\nDecompose the QR-Code into its constituent parts.\n\n\n\n\n\n","category":"function"},{"location":"#QRDecoders.getqrmatrix","page":"QRDecoders","title":"QRDecoders.getqrmatrix","text":"getqrmatrix(imgpath::AbstractString)\n\nGet the QR-code matrix from an image.\n\nNote that the input must be a standard QR-code image with or without white border, i.e. the image should not  contain any non-QR-Code information.\n\n\n\n\n\ngetqrmatrix(mat::AbstractMatrix)\n\nGet the standard QR-code matrix from an image-matrix.\n\n\n\n\n\n","category":"function"},{"location":"#QRDecoders.getqrmatrices","page":"QRDecoders","title":"QRDecoders.getqrmatrices","text":"getqrmatrices(imgpath::AbstractString)\n\nGet the standard QR-code matrices from a gif file.\n\n\n\n\n\n","category":"function"},{"location":"#Decoding-procedures","page":"QRDecoders","title":"Decoding procedures","text":"","category":"section"},{"location":"","page":"QRDecoders","title":"QRDecoders","text":"Check the version and format information","category":"page"},{"location":"","page":"QRDecoders","title":"QRDecoders","text":"QRDecoders.qrdecode_version\nQRDecoders.qrdecode_format","category":"page"},{"location":"#QRDecoders.qrdecode_version","page":"QRDecoders","title":"QRDecoders.qrdecode_version","text":"qrdecode_version(version::Int)\n\nDecode version information.(Interger to Integer)\n\n\n\n\n\nqrdecode_version(mat::AbstractMatrix; noerror=false)\n\nReturn the version of the QR-Code.(Matrix to Integer)\n\n\n\n\n\n","category":"function"},{"location":"#QRDecoders.qrdecode_format","page":"QRDecoders","title":"QRDecoders.qrdecode_format","text":"qrdecode_format(fmt::Int)::Int\n\nDecode format information.\n\n\n\n\n\nqrdecode_format(mat::AbstractMatrix; noerror=false)\n\nReturn the format of the QR-Code(ErrCorrLevel + mask).\n\n\n\n\n\n","category":"function"},{"location":"","page":"QRDecoders","title":"QRDecoders","text":"Extract message bits","category":"page"},{"location":"","page":"QRDecoders","title":"QRDecoders","text":"QRDecoders.extract_databits\nQRDecoders.deinterleave\nQRDecoders.block2bits","category":"page"},{"location":"#QRDecoders.extract_databits","page":"QRDecoders","title":"QRDecoders.extract_databits","text":"extract_databits(mat::AbstractMatrix, datapos::AbstractMatrix)\n\nExtract data bits from the QR-Code.(Inverse procedure of placedata!)\n\n\n\n\n\n","category":"function"},{"location":"#QRDecoders.deinterleave","page":"QRDecoders","title":"QRDecoders.deinterleave","text":"deinterleave(bytes::AbstractVector, ec::ErrCorrLevel, version::Int)\n\nDe-interleave the message, i.e. sperate msgblocks and ecblocks from data bits.\n\n\n\n\n\n","category":"function"},{"location":"#QRDecoders.block2bits","page":"QRDecoders","title":"QRDecoders.block2bits","text":"block2bits(blocks::AbstractVector)\n\nConvert a block to bits.\n\n\n\n\n\n","category":"function"},{"location":"","page":"QRDecoders","title":"QRDecoders","text":"Error correction","category":"page"},{"location":"","page":"QRDecoders","title":"QRDecoders","text":"QRDecoders.correct_message","category":"page"},{"location":"#QRDecoders.correct_message","page":"QRDecoders","title":"QRDecoders.correct_message","text":"correct_message(msgblock::AbstractVector\n               , ecblock::AbstractVector\n               , alg::ReedSolomonAlgorithm\n               ; noerror::Bool=false)\n\nError correction of the message block using the given algorithm.\n\nThrow DecodeError if the value noerror is true and the message need error correction.\n\n\n\n\n\n","category":"function"},{"location":"","page":"QRDecoders","title":"QRDecoders","text":"Decode message","category":"page"},{"location":"","page":"QRDecoders","title":"QRDecoders","text":"QRDecoders.decodemode\nQRDecoders.decodedata","category":"page"},{"location":"#QRDecoders.decodemode","page":"QRDecoders","title":"QRDecoders.decodemode","text":"decodemode(bits::AbstractVector)\n\nDecode mode from the bits of length 4. Note: the Byte mode and the UTF8 mode use the same mode indicator(0100).\n\n\n\n\n\n","category":"function"},{"location":"#QRDecoders.decodedata","page":"QRDecoders","title":"QRDecoders.decodedata","text":"decodedata(bits::AbstractVector, msglen::Int, ::Mode)\n\nDecode message from bits without checking the pad_bits.\n\n\n\n\n\n","category":"function"},{"location":"#Syndrome-Decoding","page":"QRDecoders","title":"Syndrome Decoding","text":"","category":"section"},{"location":"","page":"QRDecoders","title":"QRDecoders","text":"Algorithm for decoding error correction codes.","category":"page"},{"location":"","page":"QRDecoders","title":"QRDecoders","text":"ReedSolomonAlgorithm\nEuclidean\nBerlekampMassey\nRSdecoder\nberlekamp_massey_decoder\neuclidean_decoder","category":"page"},{"location":"#QRDecoders.ReedSolomonAlgorithm","page":"QRDecoders","title":"QRDecoders.ReedSolomonAlgorithm","text":"ReedSolomonAlgorithm\n\nAn abstract type for error correction algorithm of Reed Solomon code.\n\n\n\n\n\n","category":"type"},{"location":"#QRDecoders.Euclidean","page":"QRDecoders","title":"QRDecoders.Euclidean","text":"Euclidean <: ReedSolomonAlgorithm\n\nEuclidean algorithm for error correction.\n\n\n\n\n\n","category":"type"},{"location":"#QRDecoders.BerlekampMassey","page":"QRDecoders","title":"QRDecoders.BerlekampMassey","text":"BerlekampMassey <: ReedSolomonAlgorithm\n\nBerlekamp-Massey algorithm for error correction.\n\n\n\n\n\n","category":"type"},{"location":"#QRDecoders.Syndrome.RSdecoder","page":"QRDecoders","title":"QRDecoders.Syndrome.RSdecoder","text":"RSdecoder(received::Poly, nsym::Int, ::ReedSolomonAlgorithm)\n\nDecode the message polynomial using the given Reed-Solomon algorithm.\n\n\n\n\n\n","category":"function"},{"location":"#QRDecoders.Syndrome.berlekamp_massey_decoder","page":"QRDecoders","title":"QRDecoders.Syndrome.berlekamp_massey_decoder","text":"berlekamp_massey_decoder(received::Poly, erasures::AbstractVector, nsym::Int)\n\nBerlekamp-Massey algorithm, decode message polynomial from received polynomial(given erasures).\n\n\n\n\n\nberlekamp_massey_decoder(received::Poly, nsym::Int)\n\nBerlekamp-Massey algorithm, decode message polynomial from received polynomial(without erasures).\n\n\n\n\n\n","category":"function"},{"location":"#QRDecoders.Syndrome.euclidean_decoder","page":"QRDecoders","title":"QRDecoders.Syndrome.euclidean_decoder","text":"euclidean_decoder(received::Poly, erasures::AbstractVector, nsym::Int)\n\nDecode the received polynomial using the Euclidean algorithm(with erasures).\n\n\n\n\n\neuclidean_decoder(received::Poly, erasures::AbstractVector, nsym::Int)\n\nDecode the received polynomial using the Euclidean algorithm(without erasures).\n\n\n\n\n\n","category":"function"},{"location":"","page":"QRDecoders","title":"QRDecoders","text":"Tools for Syndrome Decoding.","category":"page"},{"location":"","page":"QRDecoders","title":"QRDecoders","text":"QRDecoders.Syndrome.syndrome_polynomial\nQRDecoders.Syndrome.evaluator_polynomial\nQRDecoders.Syndrome.erratalocator_polynomial\nQRDecoders.Syndrome.forney_algorithm\nQRDecoders.Syndrome.haserrors\nQRDecoders.Syndrome.fillerasures\nQRDecoders.Syndrome.Sugiyama_euclidean_divide","category":"page"},{"location":"#QRDecoders.Syndrome.syndrome_polynomial","page":"QRDecoders","title":"QRDecoders.Syndrome.syndrome_polynomial","text":"syndrome_polynomial(received::Poly, nsym::Integer)\n\nComputes the syndrome polynomial S(x) for the received polynomial where nsym is the number of syndromes.\n\n\n\n\n\n","category":"function"},{"location":"#QRDecoders.Syndrome.evaluator_polynomial","page":"QRDecoders","title":"QRDecoders.Syndrome.evaluator_polynomial","text":"evaluator_polynomial(sydpoly::Poly, errlocpoly::Poly, nsym::Int)\n\nReturn the evaluator polynomial Ω(x) where Ω(x)≡S(x)Λ(x) mod xⁿ.\n\n\n\n\n\n","category":"function"},{"location":"#QRDecoders.Syndrome.erratalocator_polynomial","page":"QRDecoders","title":"QRDecoders.Syndrome.erratalocator_polynomial","text":"erratalocator_polynomial(errpos::AbstractVector)\n\nCompute the erasures/error locator polynomial Λ(x) from the erasures/errors positions.\n\n\n\n\n\nerratalocator_polynomial(sydpoly::Poly, nsym::Int; check=false)\n\nCompute the error locator polynomial Λ(x)(without erasures). The check tag ensures that Λx can be decomposed into products of one degree polynomials.\n\n\n\n\n\nerratalocator_polynomial(sydpoly::Poly, erasures::AbstractVector, n::Int)\n\nBerlekamp-Massey algorithm, compute the error locator polynomial Λ(x)(given the erased positions). The check tag ensures that Λx can be decomposed into products of one degree polynomials.\n\n\n\n\n\n","category":"function"},{"location":"#QRDecoders.Syndrome.forney_algorithm","page":"QRDecoders","title":"QRDecoders.Syndrome.forney_algorithm","text":"forney_algorithm(Λx::Poly, Ωx::Poly, errpos::AbstractVector)\n\nForney algorithm, returns the error-corrected values. eₖ = 2^{iₖ}⋅Ω(2^{-iₖ}) / Λ'(2^{-iₖ})\n\n\n\n\n\n","category":"function"},{"location":"#QRDecoders.Syndrome.haserrors","page":"QRDecoders","title":"QRDecoders.Syndrome.haserrors","text":"haserrors(received::Poly, nsym::Int)\n\nReturns true if the received polynomial has errors. (may go undetected when the number of errors exceeds n)\n\n\n\n\n\n","category":"function"},{"location":"#QRDecoders.Syndrome.fillerasures","page":"QRDecoders","title":"QRDecoders.Syndrome.fillerasures","text":"fillerasures(received::Poly, errpos::AbstractVector, nsym::Int)\n\nForney algorithm, computes the values (error magnitude) to correct the input message.\n\nWarnning: The output polynomial might be incorrect if errpos is incomplete.\n\n\n\n\n\n","category":"function"},{"location":"#QRDecoders.Syndrome.Sugiyama_euclidean_divide","page":"QRDecoders","title":"QRDecoders.Syndrome.Sugiyama_euclidean_divide","text":"Sugiyama_euclidean_divide(r₁::Poly, r₂::Poly, upperdeg::Int)\n\nYasuo Sugiyama's adaptation of the Extended Euclidean algorithm. Find u(x), v(x) and r(x) s.t. r(x) = u(x)r₁(x) + v(x)r₂(x) where r(x) = gcd(r₁(x), r₂(x)) or deg(r(x)) ≤ upperdeg.\n\n\n\n\n\n","category":"function"},{"location":"","page":"QRDecoders","title":"QRDecoders","text":"Tools for polynomials.","category":"page"},{"location":"","page":"QRDecoders","title":"QRDecoders","text":"QRDecoders.Syndrome.polynomial_eval\nQRDecoders.Syndrome.derivative_polynomial\nQRDecoders.Syndrome.findroots\nQRDecoders.Syndrome.reducebyHorner\nQRDecoders.Syndrome.getpositions\nQRDecoders.Syndrome.extended_euclidean_divide","category":"page"},{"location":"#QRDecoders.Syndrome.polynomial_eval","page":"QRDecoders","title":"QRDecoders.Syndrome.polynomial_eval","text":"polynomial_eval(p::Poly, x::Integer)\n\nEvaluates the polynomial p(x) at x in GF256.\n\n\n\n\n\n","category":"function"},{"location":"#QRDecoders.Syndrome.derivative_polynomial","page":"QRDecoders","title":"QRDecoders.Syndrome.derivative_polynomial","text":"derivative_polynomial(p::Poly)\n\nComputes the derivative of the polynomial p.\n\n\n\n\n\n","category":"function"},{"location":"#QRDecoders.Syndrome.findroots","page":"QRDecoders","title":"QRDecoders.Syndrome.findroots","text":"findroots(p::Poly)\n\nComputes the roots of the polynomial p using Horner's method.\n\nThe output will be an empty list if p(x) contains duplicate roots or roots not in GF(256).\n\n\n\n\n\n","category":"function"},{"location":"#QRDecoders.Syndrome.reducebyHorner","page":"QRDecoders","title":"QRDecoders.Syndrome.reducebyHorner","text":"reducebyHorner(p::Poly, a::Int)\n\nHorner's rule, find the polynomial q(x) such that p(x)-(x-a)q(x) is a constant polynomial.\n\n\n\n\n\n","category":"function"},{"location":"#QRDecoders.Syndrome.getpositions","page":"QRDecoders","title":"QRDecoders.Syndrome.getpositions","text":"getpositions(Λx::Poly)\n\nCaculate positions of errors from the error locator polynomial Λx. Note that the indexs start from 0.\n\n\n\n\n\n","category":"function"},{"location":"#QRDecoders.Syndrome.extended_euclidean_divide","page":"QRDecoders","title":"QRDecoders.Syndrome.extended_euclidean_divide","text":"extended_euclidean_divide(r₁::Poly, r₂::Poly)\n\nReturn polynomials u(x) and v(x) such that u(x)r₁(x) + v(x)r₂(x) = gcd(r₁(x), r₂(x)).\n\nillustration\n\nLet\n\nr_k = u_kr_1 + v_kr_2quad k geq 2\n\nThen\n\nbeginaligned\n    r_0 = q_0r_1 + r_2 quadRightarrow r_2 = r_0 - q_0r_1 where r_0=r_2 q_0=0\n    r_1 = q_1r_2 + r_3 quadRightarrow r_3 = r_1 - q_1r_2\n    r_2 = q_2r_3 + r_4 quadRightarrow r_4 = r_2 - q_2r_3\n    vdots\n    r_k = q_kr_k+1 + r_k+2Rightarrow r_k+2 = r_k - q_kr_k+1  \n    phantom= q_kr_k+1 + r_k+2Rightarrow r_k+2 = u_kr_1 + v_kr_2 - q_k(u_k+1r_1 + v_k+1r_2)\n    phantom= q_kr_k+1 + r_k+2Rightarrow r_k+2 = (u_k - q_ku_k+1)r_1 + (v_k - q_kv_k+1)r_2\nendaligned\n\nLoop until r_t = q_tr_t+1 + 0,  then r_t+1 is the greatest common factor of  r_1 and r_2.\n\nHere we obtain the recursive formula of u_k and v_k.\n\nbeginaligned\n    u_2 v_2 = 0 1quad \n    u_3 v_3 = 1 -q_1\n    u_k+1 = u_k - q_ku_k+1quad\n    v_k+1 = v_k - q_kv_k+1\nendaligned\n\n\n\n\n\n","category":"function"},{"location":"#Error-types","page":"QRDecoders","title":"Error types","text":"","category":"section"},{"location":"","page":"QRDecoders","title":"QRDecoders","text":"InfoError\nDecodeError\nReedSolomonError","category":"page"},{"location":"#QRDecoders.InfoError","page":"QRDecoders","title":"QRDecoders.InfoError","text":"InfoError <: Exception\n\nThe non-data part of QR-matrix contains error.\n\nFor example, Finder pattern, Alignment pattern, Timing pattern,  Format information, Version information, matrix size and etc.\n\n\n\n\n\n","category":"type"},{"location":"#QRDecoders.DecodeError","page":"QRDecoders","title":"QRDecoders.DecodeError","text":"DecodeError <: Exception\n\nThe data part of QR-matrix contains error.\n\n\n\n\n\n","category":"type"},{"location":"#QRDecoders.ReedSolomonError","page":"QRDecoders","title":"QRDecoders.ReedSolomonError","text":"ReedSolomonError <: Exception\n\nAn error occurs during error-correction.\n\n\n\n\n\n","category":"type"}]
}
