immutable SplineSparse{T,I} <: Base.SparseMatrix.AbstractSparseMatrix{T,I}
    n_chunks::Int
    chunk_len::Int
    ncol::Int
    vals::Vector{T}
    cols::Vector{I}

    function SplineSparse(n_chunks, chunk_len, ncol, vals, cols)
        if length(cols)*chunk_len != length(vals)
            error("vals and cols not conformable")
        end
        new(n_chunks, chunk_len, ncol, vals, cols)
    end
end

function SplineSparse{T,I}(N, len, ncol::Int, vals::AbstractVector{T},
                           cols::AbstractVector{I})
    SplineSparse{T,I}(N, len, ncol, vals, cols)
end

n_chunks(s::SplineSparse) = s.n_chunks
chunk_len(s::SplineSparse) = s.chunk_len
Base.eltype{T}(::SplineSparse{T}) = T
ind_type{T,I}(::SplineSparse{T,I}) = I

@inline val_ix{T,I}(s::SplineSparse{T,I}, row, chunk, n) =
    s.n_chunks*s.chunk_len*(row-1) + s.chunk_len*(chunk-1) + n

@inline col_ix{N}(s::SplineSparse{N}, row, chunk) =
    s.n_chunks*(row-1) + chunk

function Base.full{T}(s::SplineSparse{T})
    nrow = _nrows(s)
    out = zeros(T, nrow, s.ncol)

    for row in 1:nrow
        for chunk in 1:s.n_chunks
            first_col = s.cols[col_ix(s, row, chunk)]
            @simd for n in 1:s.chunk_len
                out[row, first_col+n-1] = s.vals[val_ix(s, row, chunk, n)]
            end
        end
    end

    out
end

function Base.findnz{T,I}(s::SplineSparse{T,I})
    nrow = _nrows(s)
    rows = repeat(collect(I, 1:nrow), inner=[s.n_chunks*s.chunk_len])
    cols = Array(I, length(s.vals))

    for row in 1:nrow
        for chunk in 1:s.n_chunks
            first_col = s.cols[col_ix(s, row, chunk)]
            for n in 1:s.chunk_len
                cols[val_ix(s, row, chunk, n)] = first_col+n-1
            end
        end
    end

    rows, cols, s.vals
end

function Base.convert(::Type{SparseMatrixCSC}, s::SplineSparse)
    I, J, V = findnz(s)
    sparse(I, J, V, _nrows(s), s.ncol)
end

Base.sparse(s::SplineSparse) = convert(SparseMatrixCSC, s)

function Base.getindex{T}(s::SplineSparse{T}, row::Integer, cols::Integer)
    first_cols = s.cols[row]

    for chunk in 1:s.n_chunks
        first_col = s.cols[col_ix(s, row, chunk)]

        if cols < first_col || cols > (first_col + s.chunk_len-1)
            continue
        end

        n = cols-first_col+1
        return s.vals[val_ix(s, row, chunk, n)]

    end

    zero(T)
end

_nrows(s::SplineSparse) = Int(length(s.vals) / (s.n_chunks*s.chunk_len))
Base.size(s::SplineSparse) = (_nrows(s), s.ncol)
Base.size(s::SplineSparse, i::Integer) = i == 1 ? _nrows(s) :
                                         i == 2 ? s.ncol :
                                         1

function row_kron{T1,I1,T2,I2}(s1::SplineSparse{T1,I1}, s2::SplineSparse{T2,I2})

    nrow = _nrows(s1)
    _nrows(s2) == nrow || error("s1 and s2 must have same number of rows")
    N1 = s1.n_chunks
    len1 = s1.chunk_len
    N2 = s2.n_chunks
    len2 = s2.chunk_len

    # new number of chunks the length chunks times the number of chunks in
    # first matrix times number of chunks in second matrix
    N = len1*N1*N2

    # new chunk length is the chunk length from the second matrix
    len = len2

    # T and I are  easy...
    T = promote_type(T1, T2)
    I = promote_type(I1, I2)

    cols = Array(I, N*nrow)
    vals = Array(T, N*len*nrow)

    ix = 0
    c_ix = 0
    @inbounds for row in 1:nrow
        for chunk1 in 1:N1
            first_col1 = s1.cols[col_ix(s1, row, chunk1)]

            for n1 in 1:len1
                v1 = s1.vals[val_ix(s1, row, chunk1, n1)]

                for chunk2 in 1:N2
                    first_col2 = s2.cols[col_ix(s2, row, chunk2)]

                    cols[c_ix+=1] = (first_col1+n1-2)*s2.ncol+first_col2

                    @simd for n2 in 1:len2
                        vals[ix+=1] = v1*s2.vals[val_ix(s2, row, chunk2, n2)]
                    end
                end
            end
        end
    end

    SplineSparse(N, len, s1.ncol*s2.ncol, vals, cols)

end

function matvec!{T}(out::AbstractVector{T}, s::SplineSparse, v::AbstractVector)
    @inbounds for row in eachindex(out)
        val = zero(T)
        for chunk in 1:s.n_chunks
            first_col = s.cols[col_ix(s, row, chunk)]

            @simd for n in 1:s.chunk_len
                val += s.vals[val_ix(s, row, chunk, n)] * v[first_col+(n-1)]
            end
        end

        out[row] = val
    end
    out
end

function Base.(:(*)){T,I,T2}(s::SplineSparse{T,I},
                             v::AbstractVector{T2})
    size(s, 2) == size(v, 1) || throw(DimensionMismatch())

    out_T = promote_type(T, T2)
    out = Array(out_T, size(s, 1))
    matvec!(out, s, v)
end


function Base.(:(*)){T,I,T2}(s::SplineSparse{T,I},
                             m::AbstractMatrix{T2})
    size(s, 2) == size(m, 1) || throw(DimensionMismatch())

    out_T = promote_type(T, T2)
    out = Array(out_T, size(s, 1), size(m, 2))

    @inbounds for row in 1:size(s, 1)
        vals = zeros(out_T, size(m, 2))
        for chunk in 1:s.n_chunks
            first_col = s.cols[col_ix(s, row, chunk)]

            for n in 1:s.chunk_len
                s_val = s.vals[val_ix(s, row, chunk, n)]
                s_col = first_col+(n-1)

                for m_col in 1:size(m, 2)
                    vals[m_col] += s_val * m[s_col, m_col]
                end

            end
        end

        out[row, :] = vals
    end

    out

end
