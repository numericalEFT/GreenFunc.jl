using GreenFunc
using Test, StaticArrays, LinearAlgebra, Printf, Statistics, Lehmann, CompositeGrids
using JLD2, FileIO
using CodecZlib

if isempty(ARGS)
    include("test_Green.jl")
    include("test_GreenDLR.jl")
    include("test_MeshProduct.jl")
    include("test_GreenNew.jl")
    include("test_MeshGrids.jl")
else
    include(ARGS[1])
end
