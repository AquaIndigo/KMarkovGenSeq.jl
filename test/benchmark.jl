import Pkg
Pkg.activate(".")

using KMarkovGenSeq
println("Using CPU")
for i in 3:9
    @time query_freq("data/test.fa", "data/genomes/", i, "res.txt", false)    
end

println("Using GPU")
for i in 3:9
    @time query_freq("data/test.fa", "data/genomes/", i, "res.txt", true)    
end