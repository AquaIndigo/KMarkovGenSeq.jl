module KMarkovGenSeq

using FASTX, LinearAlgebra
import CUDA: copyto!, CuMatrix, CUBLAS
import CUDA

export count_freq!, db_freq, query_freq

function count_freq!(cnt::AbstractVector, seq::AbstractVector{UInt8}, len)
    mask = typemax(UInt) >> (64 - 2 * len)
    len += 1
    idx = UInt(0)
    flg = 0
    for i in 1:len-1
        @inbounds idx += seq[i]
        (seq[i] == 4) && (flg = len)
        idx <<= 2
        flg -= 1
    end
    
    for num in @view(seq[len:end])
        (num != 0x4) ? (idx += num) : (flg = len)
        (flg <= 0) && (@inbounds cnt[idx + 1] += 1)
        flg -= 1
        idx &= mask
        idx <<= 2
    end
end

function db_freq(dbseqs, len)
    map = fill(0x4, 127)
    map[['A', 'G', 'C', 'T'] .|> UInt8] .= [0x0, 0x1, 0x2, 0x3] 
    map[['a', 'g', 'c', 't'] .|> UInt8] .= [0x0, 0x1, 0x2, 0x3]
    res = zeros(Float32, 4 ^ (len + 1), length(dbseqs), 2)

    descriptions = String[] #empty array
    @sync @inbounds for (idx, dbseq) in dbseqs |> enumerate
        reader = open(FASTA.Reader, dbseq) #read one db seq
        record = FASTA.Record()
        read!(reader, record)
        (description(record) !== nothing) && push!(descriptions, description(record)) #store the descriptions
        @simd for i in record.sequence
            record.data[i] = map[record.data[i]]
        end

        Threads.@spawn begin
            count_freq!(@view(res[:, idx, 1]), @view(record.data[record.sequence]), len)
            @simd for i in record.sequence
                record.data[i] == 0x4 || (record.data[i] = 0x3 - record.data[i]) 
            end
            count_freq!(@view(res[:, idx, 2]), @view(record.data[record.sequence |> reverse]), len)        
        end

        close(reader)
    end
    @inbounds Threads.@threads for idx in 1:length(dbseqs)
        for i in 1:4:4^(len + 1) # calculate the P(Oi | O_{i+1}) = F(O)
            res[i:i+3, idx, 1] ./= sum(@view(res[i:i+3, idx, 1])) # if all 0, gives NaN
            res[i:i+3, idx, 2] ./= sum(@view(res[i:i+3, idx, 2]))
        end
    end

    res .= .-log.(res) # if 0, gives -inf; if NaN, gives NaN
    res[isnan.(res)] .= 10 # treat inf and NaN as 10
    res[isinf.(res)] .= 10 
    descriptions, res
end

function query_freq(query_seqs, dbseq_dir, len, out_file=nothing, use_gpu=false)
    batch_size = 100
    map = fill(0x4, 127)
    map[['A', 'G', 'C', 'T'] .|> UInt8] .= [0x0, 0x1, 0x2, 0x3] 
    map[['a', 'g', 'c', 't'] .|> UInt8] .= [0x0, 0x1, 0x2, 0x3] 

    query_freqs_h = zeros(Float32, 4 ^ (len + 1), batch_size)
    query_freqs_d = use_gpu ? CUDA.cu(query_freqs_h) : nothing
    
    dbseqs = dbseq_dir * "/" .* readdir(dbseq_dir)
    descriptions, db_freq_h = db_freq(dbseqs, len)
    db_freq_d = use_gpu ? CUDA.cu(db_freq_h) : nothing

    product_h = zeros(Float32, batch_size, length(dbseqs), 2)
    product_d = use_gpu ? CUDA.cu(product_h) : nothing

    reader = open(FASTA.Reader, query_seqs)

    total_query = 0
    cls_res = Int[]

    io = !isnothing(out_file) ? open(joinpath("output/", out_file), "w") : nothing
    while !eof(reader)
        cnt = 0
        # VERY IMPORTANT INITIALIZE!!!
        fill!(query_freqs_h, 0.0)
        @sync for idx in 1:batch_size
            eof(reader) && break

            record = FASTA.Record()
            read!(reader, record)
            @inbounds @simd for i in record.sequence
                record.data[i] = map[record.data[i]]
            end
            Threads.@spawn count_freq!(@view(query_freqs_h[:, idx]), @view(record.data[record.sequence]), len)
            cnt += 1
        end
        if db_freq_d !== nothing
            copyto!(query_freqs_d, query_freqs_h)
            cu_mat_mul!(@view(product_d[:, :, 1]), query_freqs_d, @view(db_freq_d[:, :, 1]))
            cu_mat_mul!(@view(product_d[:, :, 2]), query_freqs_d, @view(db_freq_d[:, :, 2]))
            copyto!(product_h, product_d)
        else
            mul!(@view(product_h[:, :, 1]), query_freqs_h', @view(db_freq_h[:, :, 1]))
            mul!(@view(product_h[:, :, 2]), query_freqs_h', @view(db_freq_h[:, :, 2]))
        end

        argmin_score = argmin(product_h, dims=(2, 3))

        append!(cls_res, argmin_score[i][2] for i in 1:cnt)
        if !isnothing(io)
            for i in 1:cnt
                println(io, i + total_query, "\t", descriptions[argmin_score[i][2]])
            end
        end
        total_query += cnt
    end
    !isnothing(io) && close(io)

    if use_gpu
        CUDA.unsafe_free!(product_d)
        CUDA.unsafe_free!(query_freqs_d)
        CUDA.unsafe_free!(db_freq_d)
    end

    total_query, cls_res, descriptions
end

"""
calculate C = A' * B
dims: m×n = (p×m)'p×n
"""
function cu_mat_mul!(C::CuMatrix, A::CuMatrix, B::CuMatrix)
    handle = CUBLAS.handle()
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

end