## Method

We use integers $A=0, G=1, C=2, T=3$ to represent each possible sequence, $N_1N_2\ldots N_k$, of length $k$ as an integer $n=\sum_1^k int(N_i)^{k+1-i}$ ranging from 0 to $4^{k}-1$. A state in the Markov model is defined as an oligonucleotide of length $k$, and each state connects to 4 other states ($({k-1})\text{mer}+A/G/C/T$). The previous state shares $k-1$ bases with the next state. Therefore, there are $4^{k+1}$ transitions in total. A genomic sequence under the $k$ th-order Markov model can be viewed as a sequence of state-transitions. The transition probabilities can be calculated for each genome in the training data set according to its Markov model as following:

$$
kMM_{i, mn}=P_{i}\left(O_{m}\rightarrow O_{n}\right)=\frac{F_{i}\left(O_{m}\rightarrow O_{n}\right)}{F_{i}\left(O_{m}\right)}=\sum_{b \in \{\text{A,G,C,T}\}}\frac{F_{i}\left(O_{m}\rightarrow O_{n}\right)}{F_{i}\left(O_{m}\rightarrow O_b\right)}
$$

where $O_m$ and $O_n$ are oligonucleotides of length $k$. If $O_n$ is ended with base $b$, we can also denote it as $O_b$. $P_{i}\left(O_{m}\rightarrow O_{n}\right)$ represents the transition probability from $O_m$ to $O_n$, $F_{i}\left(O_{m}\rightarrow O_{n}\right)$ represents observed count of transitions from $O_m$ to $O_n$ in a genomic sequence $i$ and $F(O_m)$ is the observed count of $O_m$. Since $F_{i}\left(O_{m}\right) = \sum_{b \in \{\text{A,G,C,T}\}}F_{i}\left(O_{m}\rightarrow O_b\right)$, we can get the result using only the counts of $(k+1)$mers and there's no need for us to recount $k$ mers. A $4^{k+1}$ dimension vector is created to represent each genome.  Because we have to evaluate the transformation from $k$ mer to the next base, a sequence of $k+1$ length is needed to be counted and recorded. We have to initialize a count vector with the length of $4^{k+1}$. In practice, the minus logarithm value of each transition probability is saved.

We read the sequence and count the transition states along the way, using sliding windows to improve the efficiency. A short sequence of length $l$ can be considered as $l-k$ transitions and a score $S_i$ , which represents the distance between the short sequence and a genome $i$, can be computed as following:

$$
\begin{align}
S_{i} &=-\sum_{j=0}^{l-k-1} \ln \left(P_{i}\left(O_{j} \rightarrow O_{j+1}\right)\right)\\
&=-\sum_{mn\in (k+1)mer}f_i(O_m\rightarrow O_n)\ln(P_{i}(O_{m}\rightarrow O_{n}))\\
&=-\sum_{mn\in (k+1)mer}f_i(O_m\rightarrow O_n)\ln(\sum_{b \in \{\text{A,G,C,T}\}}\frac{F_{i}\left(O_{m}\rightarrow O_{n}\right)}{F_{i}\left(O_{m}\rightarrow O_b\right)})
\end{align}
$$

where $O_j$ and $O_{j+1}$ are two oligonucleotides of length $k$, and $mn$ is $O_m$ and $O_n$ combined with the length of $k+1$, and $P_{i}(O_{j} \rightarrow O_{j+1})$ is the transition probability from $O_j$ to $O_{j+1}$ observed in the $i$-th genome. $F_{i}\left(O_{m}\rightarrow O_{n}\right)$ represents observed count of transitions from $O_m$ to $O_n$ in a genomic sequence $i$ and $f_{i}\left(O_{m}\rightarrow O_{n}\right)$ is that in a query sequence. When the transition from $O_j$ to $O_{j+1}$ does not exist in the $i$-th genome, the logarithm value of the transition probability will be set to a constant (default is 10). For each sequence, a genome in the database with the minimum score is selected as the source genome. At the end, each sequence will be annotated with the taxonomy information of its source genome. $S_i$ can be calculated by the group of $4^{k+1}$ categories of transition state, which enables that the counts of transition states in query sequence (denoted by a vector of $4^{k+1}$ length for each query sequence and a matrix $A$ of $z \times 4^{k+1}$ scale for all $z$ queries) and transition probabilities in the genome database (denoted by a vector of $4^{k+1}$ length for each genome database and a matrix $B$ of $2w \times 4^{k+1}$ scale for all $w$ databases, considering the reversed sequences) can be calculated respectively. The score can be calculated by $A*B^T$ and a scoring matix with a scale of $z\times 2w$ will return. Choose the best score in $2w$ scores of each query sequences, and the corresponding database can be the best-hit result.

![best-hit](https://raw.githubusercontent.com/AquaIndigo/KMarkovGenSeq.jl/master/docs/src/assets/image-20210517103144335.png)

The algorithm complexity is determined by the number of genomes in the database and the order of Markov Models. It can be defined as follows: 
$$
T(k,N')\sim O(N'\times 4^k)
$$
where $k$ represents the length of oligonucleotides and $N'$ stands for the number of genomes.

In practice, the score $S_i$ in equation is calculated by matrix multiplication. First, the transitions generated from each genome in the reference database are converted into a $4^{k+1}$ diemension vector. Then, a matrix can be created from all vectors generated from genomes in the reference database. These can be prebuilt. For each short metagenomic sequence, the transitions generated from it are converted into a $4^{k+1}$ dimension vector as well. Then, the scores are computed by matrix multiplication, which is done by calling the SGEMM() function of CUBLAS library. At the end, the best score is picked and the associated genome is selected as its source genome. These are done by GPUs and the taxonomy information about the source genomes is printed out by CPUs.

In our project, we use numeric vector rather than string to record the sequences and bitwise operation rather than hash table to record the $k$ th-order Markov model, so the efficiency can be improved. Also, matrix multiplication can be operated in GPUs thus our algorithm is capable even in extremely high dimensional cases.

## Improvement

The most time-consuming procedure of our project is counting the frequency of transition state, whereas matrix multiplication doesn't take long in comparison. Thus, matrix multiplication in GPU doesn't make it in speeding up in the scenario of this scale since transferring data (IO procedure of video memory and RAM) takes time, however, when in high dimensional cases, GPU computation makes sense.

Though using hash table to record the frequencies of all the sub-strings with length $k+1$ is the most straightforward way, considering the inefficiency and inconvenience in following operation, we used dynamic programming with sliding windows to count the transition states with the help of array, rather than hash table. Here, we call the sub-strings with length $k+1$ as $k$-sstr.

In order to map the sub-strings to the indices of the array, we treat the sub-strings as quadratic numbers, since there are four cases for each bit of a $k$-sstr: `A`, `G`, `C` and `T`. And we can assign $0$ to `A`, $1$ to `G`, $2$ to `C` and $3$ to `T`. Therefore, all the $k$-sstr can be one-by-one-mapped to an integer from $0$ to $4^{k+1}-1$. And for example, when $k=3$, `AAAA` can be mapped to $(0000)_4=0$,  `TTTT` can be mapped to $(3333)_4=255$，while `CGAT` mapped to $(2301)_4=2\times4^3+3\times4^2+0\times4^1+1\times4^0=177$. And we can mapped any $k$-sstr to an integer index of an array.

Given a sequence with length $l$, if using the most naive algorithm to iterate through each $k$-sstr in the sequence and then compute the integer $n=\sum_1^k int(N_i)^{k+1-i}(N_i:A=0, G=1, C=2, T=3)$ corresponding to this $k$-sstr, adding one to the number of subscript digits $n$ in the counting array. It is easy to show that such a computation requires a complexity of $\mathcal O(k(l-k))$. To speed up the counting operation, we adopt the idea of dynamic programming to see the overhead of mapping: Considering the simple fact if there is a sequence `AGCTTCG`, when we want to count the $3$-sstr of it, at the beginning, we should count the integer $n_1$ corresponding to `ss1=AGCT`, and for the next $3$-sstr `ss2=GCTT` corresponding to the integer $n_2$, we can give the following equation: 
$$
4n_1 \bmod (3333)_4 = n_2 - n(T), n(T) = 3
$$
Thus, we can use the integer of the previous $3$-sstr to compute the current $3$-sstr value, and this computation requires only $\mathcal O(1)$ overhead, so counting $k$-sstr in a sequence of length $N$ requires only $\mathcal O(N)$ time complexity.

Since we are dealing with a 4-decimal number and 4 is an integer power of 2, we can take a two-digit shift to the right instead of multiplying by 4 and use the with operation instead of the remainder operation, thus having the following implementation for this algorithm: 

```julia
function count_freq!(cnt::AbstractVector, seq::AbstractVector{UInt8}, len)
    mask = typemax(UInt) >> (64 - 2 * len) # 2 * len bits of 1
    len += 1
    idx = UInt(0)
    for i in 1:len-1 # calculate the prefix
        idx += seq[i]
        idx = idx << 2
    end
    
    for num in seq[len:end]
        idx += num
        cnt[idx] += 1
        idx = idx & mask
        idx = idx << 2
    end
end
```

Here we note that the time complexity does not increase with increasing $k$, but the running time of the actual program increases significantly with $k$. This is because, on the one hand, we need an array of length $4^{k+1}$, and the operations of creating and initializing this array grow significantly as $k$ increases; on the other hand, as $k$ increases its corresponding integer range becomes larger, the gap between two adjacent integers corresponding to $k$-sstr increases, so the localization principle of memory access is broken, and thus the increase of $k$ degrades the performance of the actual program.

Moreover, multi-threaded parallelism and synchronization in counting the transition states of a sequence also help, since the frequencies of each sequence counts are independent, so different tasks can be assigned to different processors to take full advantage of multi-core processors.

In fact, it's questionable whether matrix multiplication can help in improving the time and space efficiency since the matrix of so short query sequence can be sparse. We considered using sparse matrices in order to speed up the computation of matrix multiplication and reduce the memory overhead. However, matrix multiplication only takes about 20% of the time, so while it does speed up in this part, the sparse matrix itself is not very accessible and lowers the efficiency.

## Result

```bash
## ten complete genome
$ cat /home/faculty/ccwei/courses/2021/pab/proj1/seq_id.map | awk '!a[$2]++'
0       Frankia symbiont of Datisca glomerata chromosome, complete genome
1975    Hydrogenobaculum sp. Y04AAS1 chromosome, complete genome
3985    Candidatus Midichloria mitochondrii IricVA chromosome, complete genome
5908    Corynebacterium variabile DSM 44702 chromosome, complete genome
7943    Psychromonas ingrahamii 37 chromosome, complete genome
10002   Roseiflexus castenholzii DSM 13941 chromosome, complete genome
12024   Alteromonas macleodii str. 'Deep ecotype' chromosome, complete genome
14040   Denitrovibrio acetiphilus DSM 12809 chromosome, complete genome
16045   Sphingomonas wittichii RW1 chromosome, complete genome
18003   Baumannia cicadellinicola str. Hc (Homalodisca coagulata), complete genome

## total number of short sequences
$ grep '>' /home/faculty/ccwei/courses/2021/pab/proj1/reads.fa | wc -l
1876
```

1. the total number of short sequences in the file `reads.fa`: 1876

4. the number of groups with at least j short sequences assigned, where j = 1, 5, 10 or 50. You can describe the results in a table similar to the one below. 

| Minimum number of short sequences in a group | Number of groups |
| :------------------------------------------: | :--------------: |
|                      1                       |        10        |
|                      5                       |        10        |
|                      10                      |        10        |
|                      50                      |        8         |

5. for those groups with at least 10 short sequences assigned, list the total numbers of short sequences assigned to those groups. Then you can draw a pie chart to show the relative frequency of each group. 

- **Baumannia cicadellinicola str. Hc (Homalodisca coagulata), complete genome** contains **298** query sequences

- **Psychromonas ingrahamii 37 chromosome, complete genome** contains **22** query sequences

- **Sphingomonas wittichii RW1 chromosome, complete genome** contains **81** query sequences

- **Roseiflexus castenholzii DSM 13941 chromosome, complete genome** contains **375** query sequences

- **Hydrogenobaculum sp. Y04AAS1 chromosome, complete genome contains**  **104** query sequences

- **Alteromonas macleodii str. 'Deep ecotype' chromosome, complete genome** contains  **279** query sequences

- **Denitrovibrio acetiphilus DSM 12809 chromosome, complete genome** contains **111** query sequences

- **Frankia symbiont of Datisca glomerata chromosome, complete genome** contains **23** query sequences

- **Candidatus Midichloria mitochondrii IricVA chromosome, complete genome** contains **348** query sequences

- **Corynebacterium variabile DSM 44702 chromosome, complete genome** contains **235** query sequences

  |   Roseiflexus    | 375  |
  | :--------------: | :--: |
  |    Candidatus    | 348  |
  |    Baumannia     | 289  |
  |   Alteromonas    | 279  |
  | Corynebacterium  | 235  |
  |  Denitrovibrio   | 111  |
  | Hydrogenobaculum | 104  |
  |   Sphingomonas   |  81  |
  |     Frankia      |  23  |
  |   Psychromonas   |  22  |

![frequency](https://raw.githubusercontent.com/AquaIndigo/KMarkovGenSeq.jl/master/docs/src/assets/image-20210524202132318.png)

6. For k=3-9, which k value gives the best classification result? Please give your own criteria for “good” result, and explain why this k value gives the best result. 

|  k   | accuracy |
| :--: | :------: |
|  3   | 0.60725  |
|  4   | 0.63815  |
|  5   | 0.65865  |
|  6   | 0.68475  |
|  7   |  0.7535  |
|  8   |  0.8921  |
|  9   |  0.9837  |

Bonus: 

If you report how you use GPUs to do metagenomic sequence classification, and give the speedup compared to your CPU version, you will get a bonus of 5 points. 

```bash
Using CPU

0.151729 seconds (226.68 k allocations: 76.954 MiB)
0.192870 seconds (230.42 k allocations: 77.655 MiB)
0.311679 seconds (245.79 k allocations: 80.482 MiB)
0.553789 seconds (307.25 k allocations: 91.793 MiB, 4.00% gc time)
1.901767 seconds (553.01 k allocations: 137.254 MiB)
5.736246 seconds (1.54 M allocations: 323.395 MiB, 0.32% gc time)
27.172187 seconds (5.47 M allocations: 1.076 GiB, 0.39% gc time)

Using GPU

0.186823 seconds (268.18 k allocations: 77.272 MiB)
0.202832 seconds (315.47 k allocations: 78.642 MiB)
0.368757 seconds (513.61 k allocations: 84.308 MiB, 4.24% gc time)
0.814391 seconds (612.34 k allocations: 96.138 MiB)
2.781158 seconds (1.32 M allocations: 148.650 MiB, 0.63% gc time)
10.318713 seconds (3.99 M allocations: 360.579 MiB, 0.38% gc time)
38.076606 seconds (19.35 M allocations: 1.283 GiB, 0.22% gc time)
```

Thus, matrix multiplication in GPU doesn't make it in speeding up in the scenario of this scale since transferring data (IO procedure of Memory and RAM) takes time, however, when in high dimensional cases, GPU computation makes sense.

Since the matrix of so short query sequence can be sparse, we considered using sparse matrices in order to speed up the computation of matrix multiplication and reduce the memory overhead. With sparse array:

```
  0.435236 seconds (105.41 k allocations: 67.727 MiB, 1.84% gc time)
  1.590763 seconds (109.26 k allocations: 70.391 MiB)
  8.156021 seconds (124.62 k allocations: 77.046 MiB, 0.05% gc time)
 49.118710 seconds (186.06 k allocations: 89.670 MiB, 0.01% gc time)
104.519134 seconds (431.82 k allocations: 128.380 MiB)
141.280736 seconds (1.41 M allocations: 239.521 MiB, 0.02% gc time)
148.108550 seconds (5.35 M allocations: 695.861 MiB, 0.07% gc time)
```

Matrix multiplication only takes about 20% of the time, so while it does speed up in this part, the sparse matrix itself is not very accessible and lowers the efficiency.

