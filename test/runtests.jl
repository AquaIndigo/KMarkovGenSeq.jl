using KMarkovGenSeq
using Test

@testset "counting test" begin
    function test_func(s::String)
        bases = "AGCT"
        d = Dict{String,Int}()
        seq2idx = Dict{String,Int}()
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
            d[qs] += length(findall(qs, s; overlap=true))
        end
        for (k, v) in d
            db_freq2[seq2idx[k]] = v
        end
        for i in 1:4:64
            db_freq2[i:(i + 3), 1] ./= sum(@view(db_freq2[i:(i + 3), 1]))
        end
        db_freq2 .= .-log.(db_freq2)
        db_freq2[isnan.(db_freq2)] .= 10
        db_freq2[isinf.(db_freq2)] .= 10

        _, db_freq1 = KMarkovGenSeq.db_freq(["test_genome.fna"], 2)

        @test db_freq2 == db_freq1[:, :, 1]

        return rm("test_genome.fna")
    end
    @testset "normal sequence" begin
        open("test_genome.fna", "w") do io
            return write(io, """
                   >seq1| This is a dummy ref seq 
                   TTTATTTATTCAGAATAATAAGTTATTCATTATTTTTATACTGAAGTTTTTATACTTAGGTATTAATATA
                   AGAGTAATTAGTTTTAGTCAGTTATATTGAGTCGCAGCGCAAGCTGCAAATTAAACTTTAAATTGAAGAG
                   TTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCTTAACACATGCAAGTCGAACGGTAACAGAGAATA
                   GCTTGCTATTTTGCTGACGAGTGGCGGACGGGTGAGTAATGCTTGGGAATTTGCCTTGTTGCGGGGGACA
                   """)
        end
        s = "TTTATTTATTCAGAATAATAAGTTATTCATTATTTTTATACTGAAGTTTTTATACTTAGGTATTAATATAAGAGTAATTAGTTTTAGTCAGTTATATTGAGTCGCAGCGCAAGCTGCAAATTAAACTTTAAATTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCTTAACACATGCAAGTCGAACGGTAACAGAGAATAGCTTGCTATTTTGCTGACGAGTGGCGGACGGGTGAGTAATGCTTGGGAATTTGCCTTGTTGCGGGGGACA"
        test_func(s)
    end

    @testset "sequence with unknown (with M)" begin
        open("test_genome.fna", "w") do io
            return write(io, """
                   >seq1| This is a dummy ref seq
                   CGCAAGCTGCAAATTAAACTTTAAATTGAAGAM
                   """)
        end
        s = "CGCAAGCTGCAAATTAAACTTTAAATTGAAGAM"
        test_func(s)
    end
end
