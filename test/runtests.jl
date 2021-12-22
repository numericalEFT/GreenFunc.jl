using GreenFunc
using Test, StaticArrays, LinearAlgebra, Printf, Statistics, Lehmann

if isempty(ARGS)
    include("test_Green.jl")
    #include("interpolate.jl")
else
    include(ARGS[1])
end
