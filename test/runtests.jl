using GreenFunc
using Test, StaticArrays, LinearAlgebra, Printf, Statistics

if isempty(ARGS)
    include("test_Green.jl")
    #include("interpolate.jl")
else
    include(ARGS[1])
end
