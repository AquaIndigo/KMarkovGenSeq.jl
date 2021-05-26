using KMarkovGenSeq
using Documenter

DocMeta.setdocmeta!(KMarkovGenSeq, :DocTestSetup, :(using KMarkovGenSeq); recursive=true)

makedocs(;
    modules=[KMarkovGenSeq],
    authors="suyaoIndigo <j.c.f.gauss@sjtu.edu.cn> and contributors",
    repo="https://github.com/AquaIndigo/KMarkovGenSeq.jl/blob/{commit}{path}#{line}",
    sitename="KMarkovGenSeq.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://aquaIndigo.github.io/KMarkovGenSeq.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Report" => "report.md",
        "Tests" => "tests.md"
    ],
)

deploydocs(;
    repo="github.com/AquaIndigo/KMarkovGenSeq.jl",
)
