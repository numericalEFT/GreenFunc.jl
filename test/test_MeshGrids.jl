@testset "MeshGrids" begin
    @testset "TimeGrid" begin
        β = 50.0
        isFermi = true
        Euv = 80.0
        rtol = 1e-9
        tsym = :ph

        DLR = DLRGrid(Euv, β, rtol, isFermi, tsym)

        tg1 = MeshGrids.TimeGrid(DLR)
        tg2 = MeshGrids.TimeGrid(ImTime, DLR)
        tg3 = MeshGrids.TimeGrid(ImFreq, DLR)

        @test length(tg1) == length(DLR)
        @test length(tg2) == length(DLR)
        @test length(tg3) == length(DLR)
        # @test size(tg1) == size(DLR) 
        # @test size(tg2) == size(DLR)
        # @test size(tg3) == size(DLR)
        # TODO: seems size(DLR) = N instead of (N, ), while size(Vector{N})=(N, )
        @test tg1[1] == DLR.ω[1]
        @test tg2[1] == DLR.τ[1]
        @test tg3[1] == DLR.n[1]

        # eltype
        @test eltype(tg1) <: Real
        # the following tests should be correct when CompositeGrids updated
        # @test eltype(tg2) <: Real
        # @test eltype(tg3) <: Int

        for (ti, t) in enumerate(tg2)
            @test t == DLR.τ[ti]
            MeshGrids.locate(tg2, t) == ti
        end
        volume(tg2) == sum(volume(tg2, i) for i in 1:length(tg2))

        for (i, n) in enumerate(tg3)
            @test MeshGrids.ωn(tg3, i) ≈ DLR.ωn[i]
            MeshGrids.locate(tg3, n) == i
        end
        volume(tg3) == sum(volume(tg3, i) for i in 1:length(tg3))

    end
end