@testset "GreenFunc" begin
    @testset "MeshProduct" begin

        N1, N2 = 5, 7
        mesh1 = SimpleGrid.Uniform{Float64}([0.0, 1.0], N1)
        mesh2 = SimpleGrid.Uniform{Float64}([0.0, 1.0], N2)
        x = 3
        y = 4
        I = 14
        meshprod = MeshProduct(mesh1, mesh2)

        println("test meshprod: ", meshprod)

        @test meshprod isa MeshProduct{Tuple{typeof(mesh1),typeof(mesh2)}}
        @test GreenFunc.rank(meshprod) == 2
        @test size(meshprod, 1) == N1
        @test size(meshprod, 2) == N2
        @test size(meshprod) == (N1, N2)

        # @inferred meshprod[1]
        # @inferred meshprod[2]

        function test_linear_index(mp, x, y, i)
            @test GreenFunc.index_to_linear(meshprod, x, y) == i
            @test GreenFunc.linear_to_index(meshprod, i) == (3, 4)
            @test mp[x, y] == mp[i]
            println(mp[i])
        end

        test_linear_index(meshprod, x, y, I)

        # println(typeof(meshprod))
        # println("rank of meshprod is $(GreenFunc.rank(meshprod))")
        # println("size of mesh1 is $(size(meshprod,1)), size of mesh2 is $(size(meshprod,2))")
        # println("size of meshprod is $(size(meshprod)), length of meshprod is $(length(meshprod))")
        # println("test index_to_linear for index(($x),($y)):$(GreenFunc.index_to_linear(meshprod,x,y))")
        # println("test linear_to_index for I=$(I):$(GreenFunc.index_to_linear(meshprod,I))")
        # println("test getindex with index input:\n", meshprod[x, y])
        # println("test getindex with linearindex input:\n", meshprod[I])
        # println("Print test:\n", show(meshprod))
    end
end

