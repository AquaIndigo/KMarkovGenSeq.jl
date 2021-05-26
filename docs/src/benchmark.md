## Benchmark

```julia
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

With sparse array:

```julia
  0.435236 seconds (105.41 k allocations: 67.727 MiB, 1.84% gc time)
  1.590763 seconds (109.26 k allocations: 70.391 MiB)
  8.156021 seconds (124.62 k allocations: 77.046 MiB, 0.05% gc time)
 49.118710 seconds (186.06 k allocations: 89.670 MiB, 0.01% gc time)
104.519134 seconds (431.82 k allocations: 128.380 MiB)
141.280736 seconds (1.41 M allocations: 239.521 MiB, 0.02% gc time)
148.108550 seconds (5.35 M allocations: 695.861 MiB, 0.07% gc time)
```

## Questions

```
1. 

Total query sequences: 1876

4.
Number of groups containing at least 1 query sequences: 10
Number of groups containing at least 5 query sequences: 10
Number of groups containing at least 10 query sequences: 10
Number of groups containing at least 50 query sequences: 8

5.
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

6.
[ Info: k = 3, accuracy = 0.60725
[ Info: k = 4, accuracy = 0.63815
[ Info: k = 5, accuracy = 0.65865
[ Info: k = 6, accuracy = 0.68475
[ Info: k = 7, accuracy = 0.7535
[ Info: k = 8, accuracy = 0.8921
[ Info: k = 9, accuracy = 0.9837
The best k is 9 with accuracy 0.9837
```