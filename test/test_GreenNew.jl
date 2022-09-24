@testset "GreenNew" begin
    function test_shape(N1, N2, innermesh)
        ############# basic test ################
        mesh1 = SimpleGrid.Uniform{Float64}([0.0, 1.0], N1)
        mesh2 = SimpleGrid.Uniform{Float64}([0.0, 1.0], N2)
        if isempty(innermesh)
            g = GreenNew(mesh1, mesh2)
            @test length(g) == N1 * N2
        else
            g = GreenNew(innermesh..., mesh1, mesh2)
            @test length(g) == N1 * N2 * reduce(*, length.(innermesh))
        end
        @test size(g) == (length.(innermesh)..., N1, N2)
        @test eltype(typeof(g)) == Float64
        gc = similar(g, ComplexF64)
        @test eltype(typeof(gc)) == ComplexF64

        ############ broadcast test ###################
        g.data = rand(g.dims...)
        if isempty(innermesh)
            g2 = GreenNew(mesh1, mesh2; data=rand(g.dims...))
        else
            g2 = GreenNew(innermesh..., mesh1, mesh2; data=rand(g.dims...))
        end

        GreenFunc._check(g, g2) #check if the two GreenFuncs have the same shape

        # sum/minus/mul/div
        g3 = g .+ g2
        @test g3.data ≈ g.data .+ g2.data

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

    end

    test_shape(5, 7, ())
    test_shape(5, 7, (1:2, 1:3))

end