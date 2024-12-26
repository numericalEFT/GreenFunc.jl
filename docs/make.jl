using GreenFunc
using Documenter
using Documenter.Remotes: GitHub

DocMeta.setdocmeta!(GreenFunc, :DocTestSetup, :(using GreenFunc); recursive=true)

makedocs(;
    modules=[GreenFunc],
    authors="Kun Chen, Tao Wang, Xiansheng Cai, PengCheng Hou, and Zhiyi Li",
    repo=GitHub("numericaleft/GreenFunc.jl"),
    sitename="GreenFunc.jl",
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
        ]
    ]
)

deploydocs(;
    repo="github.com/numericalEFT/GreenFunc.jl",
    devbranch="master"
)
