# KMarkovGenSeq

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://sjtu_suyao.github.io/KMarkovGenSeq.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://sjtu_suyao.github.io/KMarkovGenSeq.jl/dev)
[![Build Status](https://github.com/sjtu_suyao/KMarkovGenSeq.jl/workflows/CI/badge.svg)](https://github.com/sjtu_suyao/KMarkovGenSeq.jl/actions)
[![Coverage](https://codecov.io/gh/sjtu_suyao/KMarkovGenSeq.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/sjtu_suyao/KMarkovGenSeq.jl)

## Usage

You should run the package with a Julia at version 1.6

To get started, you can clone this repo and start your julia repl by

```bash
$ git clone https://github.com/AquaIndigo/JLSponge.jl.git

$ cd JLSponge

$ julia
```

```julia
query_freq_cup("data/reads.fa", "data/genomes", 6)
```