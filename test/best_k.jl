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
    for k in 3:7
        _, pred_cls, des = query_freq("data/test.fa", "data/genomes", k, "res.txt", false)
        acc = sum(gt_cls[i] == pred_cls[i] for i in 1:len)
        if acc > best_acc
            best_acc = acc
            best_k = k
        end
    end
    println("The best k is $best_k with accuracy ", best_acc / len)
    @show des
end

main()

print()