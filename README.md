# KMarkovGenSeq

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://AquaIndigo.github.io/KMarkovGenSeq.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://AquaIndigo.github.io/KMarkovGenSeq.jl/dev)
[![Build Status](https://github.com/AquaIndigo/KMarkovGenSeq.jl/workflows/CI/badge.svg)](https://github.com/AquaIndigo/KMarkovGenSeq.jl/actions)
[![Coverage](https://codecov.io/gh/AquaIndigo/KMarkovGenSeq.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/AquaIndigo/KMarkovGenSeq.jl)

## Usage

You should run the package with a Julia at version 1.6

To get started, you can clone this repo and start your julia REPL by

```bash
$ git clone https://github.com/AquaIndigo/KMarkovGenSeq.jl.git

$ cd KMarkovGenSeq.jl

$ julia
```

In julia repl, you should first activate the environment:
```julia
(press ]) pkg> activate .

pkg> instantiate

(press backspace) julia> using KMarkovGenSeq

julia> query_file = ...# (the seqs you want to query)

julia> db_dir = ...# (the directory containing database seqs)

julia> out_file = ...# (the output file name)

julia> query_freq_cup(query_file, db_dir, 6, out_file)
```

Then the output file should in `output/out_file`, if `out_file` is not given,
no output will be created!

And the function returns the number of sequences queried and 
the classification results.

## 性能优化

程序最耗时的部分可能是数数据库和查询序列对应各种状态转移的频数/频率，
而计算得分反而不是非常耗时，因此性能优化主要优化的部分是频数的统计

- **动态规划（滑动窗口）求频数**
- 多线程并行+同步来数序列的频数（由于每条序列数频数是相互独立的，因而可以给不同处理机分配不同的任务来充分利用多核处理器）
- 矩阵乘法求打分（这个能否提升性能存疑，因为实际上query序列对应的向量十分稀疏）
- GPU求矩阵乘法（这个实际效果欠佳，一方面矩阵乘法未必有效，另一方面内存与显存的数据传输很费时）
