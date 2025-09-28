using Random

@testset "MeshArray" begin
    function test_shape(N1, N2, innermesh)
        ############# basic test ################
        mesh1 = SimpleGrid.Uniform{Float64}([0.0, 1.0], N1)
        mesh2 = SimpleGrid.Uniform{Float64}([0.0, 1.0], N2)
        if isempty(innermesh)
            g = MeshArray(mesh1, mesh2)
            @test length(g) == N1 * N2
        else
            g = MeshArray(innermesh..., mesh1, mesh2)
            @test length(g) == N1 * N2 * reduce(*, length.(innermesh))
        end
        show(g)
        @test size(g) == (length.(innermesh)..., N1, N2)
        @test eltype(typeof(g)) == Float64
        gc = similar(g, ComplexF64)
        @test eltype(typeof(gc)) == ComplexF64

        ############ broadcast test ###################
        rand!(g.data)
        if isempty(innermesh)
            g2 = MeshArray(mesh1, mesh2; data=rand(g.dims...))
        else
            g2 = MeshArray(innermesh..., mesh1, mesh2; data=rand(g.dims...))
        end

        MeshArrays._check(g, g2) #check if the two GreenFuncs have the same shape

        # sum/minus/mul/div
        g3 = g .+ g2
        @test g3.data ≈ g.data .+ g2.data
        @time g3 = g .+ g2 #call similar(broadcaststyle...) to create one copy
        @time g3 = g .+ g2.data #call similar(broadcaststyle...) to create one copy

        g4 = g .- g2
        @test g4.data ≈ g.data .- g2.data

        g5 = g .* g2
        @test g5.data ≈ g.data .* g2.data

        g6 = g ./ g2
        @test g6.data ≈ g.data ./ g2.data

        # inplace operation
        _g = deepcopy(g) # store a copy of g first

        g = deepcopy(_g)
        g .+= g2
        @test _g.data .+ g2.data ≈ g.data

        g = deepcopy(_g)
        g .+= g2
        println(".+= test time")
        @time g .+= g2
        @time g .+= g2.data

        g = deepcopy(_g)
        g .-= g2
        @test _g.data .- g2.data ≈ g.data

        g = deepcopy(_g)
        g .*= 2.0
        @test _g.data .* 2.0 ≈ g.data

        g = deepcopy(_g)
        g .*= g2
        @test _g.data .* g2.data ≈ g.data

        g = deepcopy(_g)
        g ./= 2.0
        @test _g.data ./ 2.0 ≈ g.data

        zg = zero(g)
        @test zg.data ≈ zero(g.data)
    end

    test_shape(5, 7, ())
    test_shape(5, 7, (1:2, 1:3))
    test_shape(10, 140, (1:2, 1:3))
end


@testset "MeshArray Mesh Type" begin
    N1, N2 = 16, 8

    mesh1 = SimpleGrid.Uniform{Float64}([0.0, 1.0], N1)
    mesh2 = SimpleGrid.Uniform{Float64}([0.0, 1.0], N2)

    g = MeshArray(mesh1, mesh2)
    @test isconcretetype(typeof(g.mesh))

    g = MeshArray(mesh1, mesh1, mesh1, mesh1, mesh1, mesh1)
    @test isconcretetype(typeof(g.mesh))
end

@testset "MeshArray Type Stability Helpers (Issue #78)" begin
    N1, N2 = 8, 6
    mesh1 = SimpleGrid.Uniform{Float64}([0.0, 1.0], N1)
    mesh2 = SimpleGrid.Uniform{Float64}([0.0, 1.0], N2)

    mesh_tuple = (mesh1, mesh2)
    stabilized = MeshArrays._stabilize_mesh_type(mesh_tuple)
    @test isconcretetype(typeof(stabilized))
    @test stabilized === mesh_tuple

    # NOTE: Cannot test the non-concrete tuple branch (line 69 in dense.jl) because
    # typeof() always returns the actual runtime type, which is concrete. The type
    # annotation ::Tuple{Any, Any} does not change the actual type of the tuple.
    # non_concrete_tuple = tuple(mesh1, mesh2)::Tuple{Any, Any}
    # stabilized_non_concrete = MeshArrays._stabilize_mesh_type(non_concrete_tuple)
    # @test isconcretetype(typeof(stabilized_non_concrete))

    non_concrete_mesh = Any[mesh1, mesh2]
    stabilized_from_vec = MeshArrays._stabilize_mesh_type(non_concrete_mesh)
    @test isconcretetype(typeof(stabilized_from_vec))

    data = rand(N1, N2)
    result = MeshArrays._create_mesharray_typed(data, mesh_tuple, Float64, 2)
    @test result isa MeshArray{Float64, 2}
    @test result.data === data

    data_float = ones(Float64, N1, N2)
    result_same_type = MeshArrays._create_mesharray_typed(data_float, mesh_tuple, Float64, 2)
    @test result_same_type isa MeshArray{Float64, 2}
    @test result_same_type.data === data_float

    data_int = ones(Int, N1, N2)
    result_converted = MeshArrays._create_mesharray_typed(data_int, mesh_tuple, Float64, 2)
    @test result_converted isa MeshArray{Float64, 2}
    @test eltype(result_converted.data) == Float64

    MT = typeof(mesh_tuple)
    @test MeshArrays._infer_mesh_type(mesh_tuple) == MT
end

@testset "MeshArray Broadcast Type Stability" begin
    N1, N2 = 10, 12
    mesh1 = SimpleGrid.Uniform{Float64}([0.0, 1.0], N1)
    mesh2 = SimpleGrid.Uniform{Float64}([0.0, 1.0], N2)

    g1 = MeshArray(mesh1, mesh2; data=rand(N1, N2))
    g2 = MeshArray(mesh1, mesh2; data=rand(N1, N2))

    g3 = similar(g1)
    g3 .= g1 .+ g2
    @test g3.data ≈ g1.data .+ g2.data

    g4 = similar(g1)
    Base.copyto!(g4, Base.Broadcast.broadcasted(+, g1, g2))
    @test g4.data ≈ g1.data .+ g2.data

    indices = CartesianIndices(g1.data)
    bcf = Base.Broadcast.flatten(Base.Broadcast.broadcasted(*, g1, 2.0))
    g5 = similar(g1)
    MeshArrays._copyto_typed!(g5, bcf, indices)
    @test g5.data ≈ g1.data .* 2.0
end

@testset "MeshArray with AbstractVector mesh" begin
    N1, N2 = 8, 6
    mesh1 = SimpleGrid.Uniform{Float64}([0.0, 1.0], N1)
    mesh2 = SimpleGrid.Uniform{Float64}([0.0, 1.0], N2)

    mesh_vector = [mesh1, mesh2]
    g = MeshArray(; mesh=mesh_vector, dtype=Float64)

    @test size(g) == (N1, N2)
    @test eltype(g.data) == Float64
    @test isconcretetype(typeof(g.mesh))
end
