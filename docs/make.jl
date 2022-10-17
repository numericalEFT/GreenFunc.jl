using GreenFunc
using Documenter

DocMeta.setdocmeta!(GreenFunc, :DocTestSetup, :(using GreenFunc); recursive=true)

makedocs(;
    modules=[GreenFunc],
    authors="Kun Chen, Tao Wang, Xiansheng Cai, PengCheng Hou, and Zhiyi Li",
    repo="https://github.com/numericaleft/GreenFunc.jl/blob/{commit}{path}#{line}",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://numericaleft.github.io/GreenFunc.jl",
        assets=String[]
    ),
    pages=[
        "Home" => "index.md",
        "Reference" => [
            "GreenFunc" => "lib/greenfunc.md",
            "MeshArrays" => "lib/mesharrays.md",
            "MeshGrids" => "lib/meshgrids.md",
            "Deprecated" => "lib/deprecated.md",
        ]
    ]
)

deploydocs(;
    repo="github.com/numericalEFT/GreenFunc.jl",
    devbranch="master"
)
