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