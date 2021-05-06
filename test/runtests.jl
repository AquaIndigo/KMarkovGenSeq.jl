using KMarkovGenSeq
using Test

@testset "counting test" begin
    ispath("data/dummy_genome") || mkdir("data/dummy_genome")
    open("data/dummy_genome/test_genome.fna", "w") do io
        write(io, """
        >seq1| This is a dummy ref seq 
        TACTGTACTTTATCAACGGGTTTACTTGTGAATAAGGCGTAAAAATCAACTATGTTATCCACAATAGATT
        GTGTAAAAAACCGTACAGCCTCAAATTTGCAAGGTGATCAAGCGGTTTTTAAAGCGTAATTTAGAATGAA
        """)
    end
    s = "TACTGTACTTTATCAACGGGTTTACTTGTGAATAAGGCGTAAAAATCAACTATGTTATCCACAATAGATTGTGTAAAAAACCGTACAGCCTCAAATTTGCAAGGTGATCAAGCGGTTTTTAAAGCGTAATTTAGAATGAA"
    bases = "AGCT"
    d = Dict{String, Int}()
    seq2idx = Dict{String, Int}()
    idx = 0
    db_freq2 = zeros(Float32, 64, 1)
    for b1 in bases, b2 in bases, b3 in bases 
        qs = b1 * b2 * b3
        d[qs] = 0
        idx += 1
        seq2idx[qs] = idx
    end
    for b1 in bases, b2 in bases, b3 in bases 
        qs = b1 * b2 * b3
        d[qs] += findall(qs, s, overlap=true) |> length
    end
    for (k, v) in d
        db_freq2[seq2idx[k]] = v
    end
    for i in 1:4:64
        db_freq2[i:i+3, 1] ./= sum(@view(db_freq2[i:i+3, 1]))
    end
    db_freq2 .= .-log.(db_freq2)
    db_freq2[isnan.(db_freq2)] .= 10
    db_freq2[isinf.(db_freq2)] .= 10 

    desc, db_freq1 = KMarkovGenSeq.db_freq(["data/dummy_genome/test_genome.fna"], 2)
    
    @test db_freq2 == db_freq1

    rm("data/dummy_genome", recursive=true)
end
