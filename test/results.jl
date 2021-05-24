using KMarkovGenSeq

function main()
    cnt, pred_cls, des = query_freq("data/reads.fa", "data/genomes", 9, nothing, false)
    len = length(des)
    println("Total query sequences: $cnt")
    counter = fill(0, len)
    for cls in pred_cls
        counter[cls] += 1
    end
    println()
    for threshold in [1, 5, 10, 50]
        println("Number of groups containing at least $threshold query sequences: ", count(>=(threshold), counter))
    end
    for i in eachindex(counter)
        if counter[i] >= 10
            println(des[i], "contains\t", counter[i], "\tquery sequences")
        end
    end
    println()
    @show(counter)
end

main()

"""
Total query sequences: 1876

Number of groups containing at least 1 query sequences: 10
Number of groups containing at least 5 query sequences: 10
Number of groups containing at least 10 query sequences: 10
Number of groups containing at least 50 query sequences: 8

Baumannia cicadellinicola str. Hc (Homalodisca coagulata), complete genomecontains      298     query sequences
Psychromonas ingrahamii 37 chromosome, complete genomecontains  22      query sequences
Sphingomonas wittichii RW1 chromosome, complete genomecontains  81      query sequences
Roseiflexus castenholzii DSM 13941 chromosome, complete genomecontains  375     query sequences
Hydrogenobaculum sp. Y04AAS1 chromosome, complete genomecontains        104     query sequences
Alteromonas macleodii str. 'Deep ecotype' chromosome, complete genomecontains   279     query sequences
Denitrovibrio acetiphilus DSM 12809 chromosome, complete genomecontains 111     query sequences
Frankia symbiont of Datisca glomerata chromosome, complete genomecontains       23      query sequences
Candidatus Midichloria mitochondrii IricVA chromosome, complete genomecontains  348     query sequences
Corynebacterium variabile DSM 44702 chromosome, complete genomecontains 235     query sequences

counter = [298, 22, 81, 375, 104, 279, 111, 23, 348, 235]
"""