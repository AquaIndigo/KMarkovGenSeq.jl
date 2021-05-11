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