immutable SplineSparse{N,len,T,I} <: Base.SparseMatrix.AbstractSparseMatrix{T,I}
    ncol::Int
    col::Vector{I}
    vals::Vector{T}

    function SplineSparse(ncol, col, vals)
        if length(col)*N*len != length(vals)
            error("vals and cols not conformable")
        end
        new(ncol, col, vals)
    end
end

ind_type(I::DataType) = isempty(I.parameters) ? I : promote_type(I.parameters...)

n_chunks{N}(::SplineSparse{N}) = N
chunk_len{N,len}(::SplineSparse{N,len}) = len
Base.eltype{N,len,T}(::SplineSparse{N,len,T}) = T
ind_type{N,len,T,I}(::SplineSparse{N,len,T,I}) = I

# TODO: I can't get this to be type stable...
flat_ind_type_impl{N,len,T,I<:Tuple}(::Type{SplineSparse{N,len,T,I}}) =
    :(reduce(promote_type, I.parameters))

flat_ind_type_impl{N,len,T,I}(::Type{SplineSparse{N,len,T,I}}) = :(I)

@generated flat_ind_type{N,len,T,I}(s::SplineSparse{N,len,T,I}) =
    flat_ind_type_impl(s)

@inline val_ix{N,len,T,I}(s::SplineSparse{N,len,T,I}, row, chunk, n) =
    N*len*(row-1) + len*(chunk-1) + n

function Base.full{N,len,T,I}(s::SplineSparse{N,len,T,I})
    nrow = length(s.col)
    out = zeros(T, nrow, s.ncol)

    # ix = 0
    for row in 1:nrow
        first_cols = s.col[row]
        for chunk in 1:N
            first_col = first_cols[chunk]
            @simd for n in 1:len
                # ix += 1
                out[row, first_col+n-1] = s.vals[val_ix(s, row, chunk, n)]
            end
        end
    end

    out
end

function Base.findnz{N,len,T,I}(s::SplineSparse{N,len,T,I})
    nrow = length(s.col)
    rows = repeat(collect(flat_ind_type(s), 1:nrow), inner=[N*len])
    cols = Array(flat_ind_type(s), length(s.vals))

    for row in 1:nrow
        first_cols = s.col[row]
        for chunk in 1:N
            first_col = first_cols[chunk]
            for n in 1:len
                cols[val_ix(s, row, chunk, n)] = first_col+n-1
            end
        end
    end

    rows, cols, s.vals
end

function Base.convert(::Type{SparseMatrixCSC}, s::SplineSparse)
    I, J, V = findnz(s)
    sparse(I, J, V, length(s.col), s.ncol)
end

Base.sparse(s::SplineSparse) = convert(SparseMatrixCSC, s)

function Base.getindex{N,len,T}(s::SplineSparse{N,len,T}, row::Integer, col::Integer)
    first_cols = s.col[row]

    for chunk in 1:N
        first_col = first_cols[chunk]

        if col < first_col || col > (first_col + len-1)
            continue
        end

        n = col-first_col+1
        return s.vals[val_ix(s, row, chunk, n)]

    end

    zero(T)
end

Base.size(s::SplineSparse) = (length(s.col), s.ncol)
Base.size(s::SplineSparse, i::Integer) = i == 1 ? length(s.col) :
                                         i == 2 ? s.ncol :
                                         1

function row_kron{N1,len1,T1,I1,N2,len2,T2,I2}(s1::SplineSparse{N1,len1,T1,I1},
                                               s2::SplineSparse{N2,len2,T2,I2})

    nrow = length(s1.col)
    length(s2.col) == nrow || error("s1 and s2 must have same number of rows")

    # new number of chunks the length chunks times the number of chunks in
    # first matrix times number of chunks in second matrix
    N = len1*N1*N2

    # new chunk length is the chunk length from the second matrix
    len = len2

    # T is easy...
    T = promote_type(T1, T2)

    # I is almost as easy
    single_I = promote_type(flat_ind_type(s1), flat_ind_type(s2))
    I = NTuple{N,single_I}

    col = Array(I, nrow)
    col_buf = Array(single_I, N)
    vals = Array(T, N*len*nrow)

    ix = 0
    for row in 1:nrow
        # extract first columns for this row from each matrix
        first_cols1 = s1.col[row]
        first_cols2 = s2.col[row]

        col_ix = 0

        for chunk1 in 1:N1
            first_col1 = first_cols1[chunk1]

            for n1 in 1:len1
                v1 = s1.vals[val_ix(s1, row, chunk1, n1)]

                for chunk2 in 1:N2
                    first_col2 = first_cols2[chunk2]

                    col_ix += 1
                    col_buf[col_ix] = (first_col1+n1-2)*s2.ncol+first_col2

                    for n2 in 1:len2
                        v2 = s2.vals[val_ix(s2, row, chunk2, n2)]
                        v = v1*v2

                        vals[ix+=1] = v
                    end
                end
            end
        end
        col[row] = tuple(col_buf...)
    end

    SplineSparse{N,len,T,I}(s1.ncol*s2.ncol, col, vals)

end

function Base.(:(*)){N,len,T,I,T2}(s::SplineSparse{N,len,T,I},
                                   v::AbstractVector{T2})
    size(s, 2) == size(v, 1) || throw(DimensionMismatch())
    out_T = promote_type(T, T2)
    out = Array(out_T, size(s, 1))

    @inbounds for row in eachindex(out)
        val = zero(out_T)
        first_cols = s.col[row]

        for chunk in 1:N
            first_col = first_cols[chunk]

            @simd for n in 1:len
                val += s.vals[val_ix(s, row, chunk, n)] * v[first_col+(n-1)]
            end
        end

        out[row] = val
    end

    out

end
