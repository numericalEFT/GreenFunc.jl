using GreenFunc
using Documenter

DocMeta.setdocmeta!(GreenFunc, :DocTestSetup, :(using GreenFunc); recursive=true)

makedocs(;
    modules=[GreenFunc],
    authors="Tao Wang, Xiansheng Cai",
    repo="https://github.com/fsxbhyy/GreenFunc.jl/blob/{commit}{path}#{line}",
    sitename="GreenFunc.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://fsxbhyy.github.io/GreenFunc.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/numericalEFT/GreenFunc.jl",
    devbranch="master",
)
