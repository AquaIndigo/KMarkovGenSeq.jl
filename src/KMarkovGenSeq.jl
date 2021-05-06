module KMarkovGenSeq

using FASTX, CUDA, LinearAlgebra

export count_freq!, db_freq, query_freq, query_freq_cup, to

function count_freq!(cnt::AbstractVector, seq::AbstractVector{UInt8}, len)
    mask = typemax(UInt) >> (64 - 2 * len)
    len += 1
    idx = UInt(0)
    for i in 1:len-1
        @inbounds idx += seq[i]
        idx <<= 2
    end
    for num in @view(seq[len:end])
        idx += num
        @inbounds cnt[idx + 1] += 1
        idx &= mask
        idx <<= 2
    end
end

function db_freq(dbseqs, len)
    map = zeros(UInt8, 127)
    map[['A', 'G', 'C', 'T'] .|> UInt8] .= [0x0, 0x1, 0x2, 0x3] 
    res = zeros(Float32, 4 ^ (len + 1), length(dbseqs))

    descriptions = String[] #empty array

    @inbounds for (idx, dbseq) in dbseqs |> enumerate
        reader = open(FASTA.Reader, dbseq) #read one db seq
        record = FASTA.Record()
        read!(reader, record)
        push!(descriptions, description(record)) #store the descriptions
        @simd for i in record.sequence
            record.data[i] = map[record.data[i]]
        end
        
        count_freq!(@view(res[:, idx]), @view(record.data[record.sequence]), len)
        
        for i in 1:4:4^(len + 1) # calculate the P(Oi | O_{i+1}) = F(O)
            res[i:i+3, idx] ./= sum(@view(res[i:i+3, idx])) # if all 0, gives NaN
        end

        close(reader)
    end
    res .= .-log.(res) # if 0, gives inf; if NaN, gives NaN
    res[isnan.(res)] .= 10 # treat inf and NaN as 10
    res[isinf.(res)] .= 10 
    descriptions, res
end

"""
calculate C = A' * B
dims: m×n = (p×m)'p×n
"""

function cu_mat_mul!(C::CuMatrix, A::CuMatrix, B::CuMatrix)
    handle = CUBLAS.cublasCreate()
    row_A, col_A = size(A)
    row_B, col_B = size(B)
    CUBLAS.cublasSgemm_v2(
        handle, 
        CUBLAS.CUBLAS_OP_T,
        CUBLAS.CUBLAS_OP_N,
        col_A, col_B, row_B,
        1.0,
        A, row_A,
        B, row_B,
        0.0,
        C, col_A
    )
end

function query_freq(query_seqs, dbseq_dir, len)
    batch_size = 100
    map = zeros(UInt8, 127)
    map[['A', 'G', 'C', 'T'] .|> UInt8] .= [0x0, 0x1, 0x2, 0x3] 

    query_freqs_h = zeros(Float32, 4 ^ (len + 1), batch_size)
    query_freqs_d = CUDA.zeros(Float32, 4 ^ (len + 1), batch_size)
    

    dbseqs = dbseq_dir * "/" .* readdir(dbseq_dir)
    descriptions, db_freq_h = db_freq(dbseqs, len)
    db_freq_d = cu(db_freq_h)

    product_h = zeros(Float32, batch_size, length(dbseqs))
    product_d = CUDA.zeros(Float32, batch_size, length(dbseqs))

    reader = open(FASTA.Reader, query_seqs)
    total_query = 0
    io = open("output/res.txt", "w")
    while !eof(reader)
        cnt = 0
        for idx in 1:batch_size
            eof(reader) && break

            record = FASTA.Record()
            read!(reader, record)
            @inbounds @simd for i in record.sequence
                record.data[i] = map[record.data[i]]
            end

            count_freq!(@view(query_freqs_h[:, idx]), @view(record.data[record.sequence]), len)
            cnt += 1
        end
    
        copyto!(query_freqs_d, query_freqs_h)
        cu_mat_mul!(product_d, query_freqs_d, db_freq_d);
        
        copyto!(product_h, product_d)

        argmin_score = argmin(product_h, dims=2)

        for i in 1:cnt
            println(io, i + total_query, "\t", descriptions[argmin_score[i][2]])
        end

        total_query += cnt
    end
    close(io)
    CUDA.unsafe_free!(product_d)
    CUDA.unsafe_free!(query_freqs_d)
    CUDA.unsafe_free!(db_freq_d)
end

function query_freq_cup(query_seqs, dbseq_dir, len)
    map = zeros(UInt8, 127)
    map[['A', 'G', 'C', 'T'] .|> UInt8] .= [0x0, 0x1, 0x2, 0x3] 

    query_freqs_h = zeros(Float32, 4 ^ (len + 1), 100)
   
    dbseqs = dbseq_dir * "/" .* readdir(dbseq_dir)
    descriptions, db_freq_h = db_freq(dbseqs, len)
    
    product_h = zeros(Float32, 100, length(dbseqs))

    reader = open(FASTA.Reader, query_seqs)
    total_query = 0
    io = open("output/res.txt", "w")
    while !eof(reader)
        cnt = 0
        for idx in 1:100
            eof(reader) && break

            record = FASTA.Record()
            read!(reader, record)
            @inbounds @simd for i in record.sequence
                record.data[i] = map[record.data[i]]
            end

            count_freq!(@view(query_freqs_h[:, idx]), @view(record.data[record.sequence]), len)
            cnt += 1
        end

        # @time BLAS.gemm!('T', 'N', 1.0f0, query_freqs_h, db_freq_h, 0.0f0, product_h)
        mul!(product_h, query_freqs_h', db_freq_h)

        argmin_score = argmin(product_h, dims=2)

        for i in 1:cnt
            println(io, i + total_query, "\t", descriptions[argmin_score[i][2]])
        end
        total_query += cnt
    end
    close(io)
end

end