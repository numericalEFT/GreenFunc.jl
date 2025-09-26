#!/usr/bin/env julia
# Broadcast operation benchmark for Issue #78

using Pkg
Pkg.activate(dirname(@__DIR__) * "/GreenFunc.jl")

using GreenFunc
using GreenFunc.MeshArrays
using CompositeGrids

println("Issue #78 - Broadcast Operation Benchmark")
println("=" ^ 60)

# Create grids for test
k1 = CompositeGrids.SimpleG.Uniform([0.0, 1.0], 10)
k2 = CompositeGrids.SimpleG.Uniform([0.0, 1.0], 10)

# Create MeshArrays
ma1 = MeshArray(k1, k2; dtype=Float64)
ma2 = MeshArray(k1, k2; dtype=Float64)

# Fill with data
ma1.data .= rand(size(ma1.data)...)
ma2.data .= rand(size(ma2.data)...)

println("\nTesting broadcast operations on 10x10 MeshArray:")
println("-" ^ 40)

# Warmup
ma1 .+= ma2
ma1 .*= 2.0

# Reset
ma1.data .= rand(size(ma1.data)...)

# Measure broadcast allocations
println("\nIn-place broadcast operations:")
alloc_add = @allocated ma1 .+= ma2
println("ma1 .+= ma2: $(alloc_add) bytes")

alloc_mul = @allocated ma1 .*= 2.0
println("ma1 .*= 2.0: $(alloc_mul) bytes")

alloc_scalar = @allocated ma1 .+= 1.0
println("ma1 .+= 1.0: $(alloc_scalar) bytes")

# Test with larger arrays
println("\nTesting with 100x100 MeshArray:")
println("-" ^ 40)

k1_large = CompositeGrids.SimpleG.Uniform([0.0, 1.0], 100)
k2_large = CompositeGrids.SimpleG.Uniform([0.0, 1.0], 100)

ma1_large = MeshArray(k1_large, k2_large; dtype=Float64)
ma2_large = MeshArray(k1_large, k2_large; dtype=Float64)

ma1_large.data .= rand(size(ma1_large.data)...)
ma2_large.data .= rand(size(ma2_large.data)...)

# Warmup
ma1_large .+= ma2_large

# Reset and measure
ma1_large.data .= rand(size(ma1_large.data)...)

alloc_add_large = @allocated ma1_large .+= ma2_large
println("ma1 .+= ma2 (100x100): $(alloc_add_large) bytes")

alloc_mul_large = @allocated ma1_large .*= 2.0
println("ma1 .*= 2.0 (100x100): $(alloc_mul_large) bytes")

println("\nConclusion:")
if alloc_add == 0 && alloc_mul == 0
    println("✓ Broadcast operations are allocation-free for small arrays")
else
    println("⚠ Broadcast operations still allocate for small arrays")
end

if alloc_add_large < 1000
    println("✓ Large array broadcast allocations minimal: $(alloc_add_large) bytes")
else
    println("⚠ Large array broadcast allocations: $(round(alloc_add_large / 1024, digits=2)) KB")
end