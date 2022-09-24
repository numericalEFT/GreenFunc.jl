@testset "MeshGrids" begin

    β = 50.0
    isFermi = true
    Euv = 80.0
    rtol = 1e-9
    tsym = :none

    DLR = DLRGrid(Euv, β, rtol, isFermi, tsym)

    # @testset "TimeGrid" begin

    #     tg1 = MeshGrids.TimeGrid(DLR)
    #     tg2 = MeshGrids.TimeGrid(ImTime, DLR)
    #     tg3 = MeshGrids.TimeGrid(ImFreq, DLR)

    #     @test length(tg1) == length(DLR)
    #     @test length(tg2) == length(DLR)
    #     @test length(tg3) == length(DLR)
    #     # @test size(tg1) == size(DLR) 
    #     # @test size(tg2) == size(DLR)
    #     # @test size(tg3) == size(DLR)
    #     # TODO: seems size(DLR) = N instead of (N, ), while size(Vector{N})=(N, )
    #     @test tg1[1] == DLR.ω[1]
    #     @test tg2[1] == DLR.τ[1]
    #     @test tg3[1] == DLR.n[1]

    #     # eltype
    #     @test eltype(tg1) <: Real
    #     # the following tests should be correct when CompositeGrids updated
    #     # @test eltype(tg2) <: Real
    #     # @test eltype(tg3) <: Int

    #     for (ti, t) in enumerate(tg2)
    #         @test t == DLR.τ[ti]
    #         MeshGrids.locate(tg2, t) == ti
    #     end
    #     volume(tg2) == sum(volume(tg2, i) for i in 1:length(tg2))

    #     for (i, n) in enumerate(tg3)
    #         @test MeshGrids.ωn(tg3, i) ≈ DLR.ωn[i]
    #         MeshGrids.locate(tg3, n) == i
    #     end
    #     volume(tg3) == sum(volume(tg3, i) for i in 1:length(tg3))

    # end

    @testset "ImTime Grid" begin
        tg2 = MeshGrids.ImTime(β; grid=DLR.τ)

        @test length(tg2) == length(DLR)
        # @test size(tg1) == size(DLR) 
        # @test size(tg2) == size(DLR)
        # @test size(tg3) == size(DLR)
        # TODO: seems size(DLR) = N instead of (N, ), while size(Vector{N})=(N, )
        @test tg2[1] == DLR.τ[1]

        # eltype
        @test eltype(typeof(tg2)) == Float64

        for (ti, t) in enumerate(tg2)
            @test t == DLR.τ[ti]
            @test MeshGrids.locate(tg2, t) == ti
        end
        @test volume(tg2) ≈ sum(volume(tg2, i) for i in 1:length(tg2))
    end

    @testset "ImFreq Grid" begin
        tg2 = MeshGrids.ImFreq(β, MeshGrids.FERMI; grid=DLR.n)
        println(tg2)

        @test length(tg2) == length(DLR)
        # @test size(tg2) == size(DLR)
        # TODO: seems size(DLR) = N instead of (N, ), while size(Vector{N})=(N, )
        @test tg2[1] == DLR.ωn[1]

        # eltype
        @test eltype(typeof(tg2)) == Int

        for (ti, t) in enumerate(tg2)
            @test tg2.grid[ti] == DLR.n[ti]
            @test tg2[ti] ≈ DLR.ωn[ti] #DLR.ωn is read from files, cannot exactly match (2n+1)/β exactly
            @test MeshGrids.locate(tg2, t) == ti
            @test MeshGrids.locate(tg2, tg2.grid[ti]) == ti
        end
        @test volume(tg2) == sum(volume(tg2, i) for i in 1:length(tg2))
    end

    @testset "DLRFreq Grid" begin
        tg2 = MeshGrids.DLRFreq(β, MeshGrids.FERMI; dlr=DLR)
        println(tg2)

        @test length(tg2) == length(DLR)
        # @test size(tg2) == size(DLR)
        # TODO: seems size(DLR) = N instead of (N, ), while size(Vector{N})=(N, )
        @test tg2[1] == DLR.ω[1]

        # eltype
        @test eltype(typeof(tg2)) == Int

        for (ti, t) in enumerate(tg2)
            @test tg2.grid[ti] ≈ DLR.ω[ti]
            @test tg2[ti] ≈ DLR.ω[ti] #DLR.ωn is read from files, cannot exactly match (2n+1)/β exactly
            @test MeshGrids.locate(tg2, t) == ti
        end
        @test volume(tg2) ≈ sum(volume(tg2, i) for i in 1:length(tg2))
    end
end