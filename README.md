# GreenFunc

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://numericalEFT.github.io/GreenFunc.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://numericalEFT.github.io/GreenFunc.jl/dev)
[![Build Status](https://github.com/numericalEFT/GreenFunc.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/numericalEFT/GreenFunc.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/numericalEFT/GreenFunc.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/numericalEFT/GreenFunc.jl)

This library provides structures of different types of Green's function. Here we give a introduction of these Green's function.
## Features

We provide the following containers to save different Green's functions:

 -One body Green's function that has a built-in discrete Lehamnn representation (DLR),  which is a generic and  compact representation of Green's functions proposed in the Ref. [1]. 

For all Green's functions we provide the following manipulations:

- Fast and accurate Fourier transform between the imaginary-time domain and the Matsubara-frequency domain.

- Fast and accurate Green's function interpolation.

## Installation
This package has been registered. So, simply type `import Pkg; Pkg.add("GreenFunc")` in the Julia REPL to install.

## Basic Usage

```julia
using AbstractTrees, StaticArrays, ExpressionTree

# Define the type of the weight of a 4-vertex diagram
const Weight = SVector{2,Float64} 

# Define the possibilities to build a 4-vertex diagram from the left and right 4-vertex subdiagrams.
chan = [Parquet.T, Parquet.U, Parquet.S] 

# Generate the parameter for the expression tree. The second argument lists the possible number of imaginary-time variables in the bare 4-vertex (namely, the bare interaction of your model). For example, the instaneous Coulomb interaction only has one time variable, while the retared effective interaction has two time variables.
para = Parquet.Para(chan, [1, 2]) 

# Generate an expression tree for a set of 4-vertex diagrams with loop order 1, initial imaginary-time index 1, and the parameter set para.
ver4 = Parquet.Ver4{Weight}(1, 1, para) 

# use AbstractTrees interface to print/manipulate the tree
print_tree(ver4)

# [println(node) for node in Leaves(ver4)]  #print all loopNum=0 ver4
# [println(node) for node in PostOrderDFS(ver4)]  # iterator ver4 in depth-first search (children before parents)
# [println(node) for node in PreOrderDFS(ver4)]  # iterator ver4 (parents before children)

# You can also use ete3 python3 package to visualize tree
Parquet.showTree(ver4, para, verbose = 1, depth = 3)

# You can also print tree to a newick format file, there are many visualization software for the newick format
io = open("./test.newick", "w")
write(io, Parquet.newick(ver4))
close(io)
```

