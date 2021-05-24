using KMarkovGenSeq, FASTX

function ground_truth_classification()
    descriptions_dict = Dict{String, Int}()
    for (idx, dbseq) in readdir("data/genomes") |> enumerate
        db_seq_path = joinpath("data/genomes", dbseq)
        reader = open(FASTA.Reader, db_seq_path) #read one db seq
        record = FASTA.Record()
        read!(reader, record)
        descriptions_dict[description(record)] = idx
        close(reader)
    end

    gt_cls = Int[]
    for line in eachline("data/seq_id.map")
        push!(gt_cls, descriptions_dict[split(line, '\t')[end]])
    end
    gt_cls
end

function main()
    gt_cls = ground_truth_classification()
    len = gt_cls |> length
    best_k, best_acc = 0, 0
    des = String[]
    for k in 3:9
        _, pred_cls, des = query_freq("data/test.fa", "data/genomes", k, nothing, false)
        acc = sum(gt_cls[i] == pred_cls[i] for i in 1:len)
        if acc > best_acc
            best_acc = acc
            best_k = k
        end
        @info "k = $k, accuracy = $(acc/len)"
    end
    println("The best k is $best_k with accuracy ", best_acc / len)

end

main()

"""
[ Info: k = 3, accuracy = 0.60725
[ Info: k = 4, accuracy = 0.63815
[ Info: k = 5, accuracy = 0.65865
[ Info: k = 6, accuracy = 0.68475
[ Info: k = 7, accuracy = 0.7535
[ Info: k = 8, accuracy = 0.8921
[ Info: k = 9, accuracy = 0.9837
The best k is 9 with accuracy 0.9837
"""