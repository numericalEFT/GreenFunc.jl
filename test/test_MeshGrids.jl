@testset "MeshGrids" begin

    β = 50.0
    isFermi = true
    Euv = 80.0
    rtol = 1e-9
    tsym = :none

    DLR = DLRGrid(Euv, β, rtol, isFermi, tsym)

    @testset "ImTime Grid" begin

        tg2 = MeshGrids.ImTime(β, isFermi) # check default reverse
        @test length(tg2) > 0

        tg2 = MeshGrids.ImTime(β, isFermi; grid=DLR.τ)

        @test length(tg2) == length(DLR)
        @test size(tg2) == (length(DLR),)
        @test tg2[1] == DLR.τ[1]

        # eltype
        @test eltype(typeof(tg2)) == Float64

        # test locate and volume
        for (ti, t) in enumerate(tg2)
            @test t == DLR.τ[ti]
            @test MeshGrids.locate(tg2, t) == ti
        end
        @test volume(tg2) ≈ sum(volume(tg2, i) for i in 1:length(tg2))

        # test floor
        δ = 1e-6
        N = length(tg2)
        @test floor(tg2, tg2[1] - δ) == 1
        @test floor(tg2, tg2[1] + δ) == 1
        for i in 2:N-1
            @test floor(tg2, tg2[i] - δ) == i - 1
            @test floor(tg2, tg2[i] + δ) == i
        end
        @test floor(tg2, tg2[N] - δ) == N - 1
        @test floor(tg2, tg2[N] + δ) == N - 1
    end

    @testset "ImTime Grid Reversed" begin
        tg2 = MeshGrids.ImTime(β, isFermi; grid=reverse(DLR.τ))
        println(typeof(tg2))
        println(tg2)

        @test length(tg2) == length(DLR)
        @test size(tg2) == (length(DLR),)
        @test tg2[1] == DLR.τ[end]
        @test tg2[end] == DLR.τ[1]
        @test tg2.grid[1] == DLR.τ[1]
        @test tg2.grid[end] == DLR.τ[end]

        # eltype
        @test eltype(typeof(tg2)) == Float64

        for (ti, t) in enumerate(tg2)
            @test t == DLR.τ[length(tg2)-ti+1]
            # @test MeshGrids.locate(tg2, t) == ti #TODO: add test for locate
        end
        @test volume(tg2) ≈ sum(volume(tg2, i) for i in 1:length(tg2))

        # test floor
        δ = 1e-6
        N = length(tg2)
        @test floor(tg2, tg2[1] - δ) == 1
        @test floor(tg2, tg2[1] + δ) == 1
        for i in 2:N-1
            # reversed 
            @test floor(tg2, tg2[i] - δ) == i
            @test floor(tg2, tg2[i] + δ) == i - 1
        end
        @test floor(tg2, tg2[N] - δ) == N - 1
        @test floor(tg2, tg2[N] + δ) == N - 1

    end

    @testset "ImFreq Grid" begin

        tg2 = MeshGrids.ImFreq(β, isFermi) # check default reverse
        @test length(tg2) > 0

        tg2 = MeshGrids.ImFreq(β, isFermi; grid=DLR.n)
        println(tg2)

        @test length(tg2) == length(DLR)
        @test size(tg2) == (length(DLR),)
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

        MeshGrids.matfreq(tg2) ≈ DLR.ωn

        # test floor
        # notice for imfreq tg2[i]!=tg2.grid[i]
        δ = 1e-6
        N = length(tg2)
        @test floor(tg2, tg2.grid[1] - δ) == 1
        @test floor(tg2, tg2.grid[1] + δ) == 1
        for i in 2:N-1
            @test floor(tg2, tg2.grid[i] - δ) == i - 1
            @test floor(tg2, tg2.grid[i] + δ) == i
        end
        @test floor(tg2, tg2.grid[N] - δ) == N - 1
        @test floor(tg2, tg2.grid[N] + δ) == N - 1

    end

    @testset "ImFreq Grid Reversed" begin

        tg2 = MeshGrids.ImFreq(β, isFermi; grid=reverse(DLR.n))
        # println(typeof(tg2))
        # println(tg2)

        @test length(tg2) == length(DLR)
        @test size(tg2) == (length(DLR),)
        @test tg2[1] == DLR.ωn[end]
        @test tg2[end] == DLR.ωn[1]
        @test tg2.grid[1] == DLR.n[1]
        @test tg2.grid[end] == DLR.n[end]

        # eltype
        @test eltype(typeof(tg2)) == Int

        for (ti, t) in enumerate(tg2)
            @test tg2.grid[ti] == DLR.n[ti]
            @test tg2[ti] ≈ DLR.ωn[length(DLR)-ti+1] #DLR.ωn is read from files, cannot exactly match (2n+1)/β exactly
            # @test MeshGrids.locate(tg2, t) == ti #TODO: add test for locate
            @test MeshGrids.locate(tg2, tg2.grid[ti]) == ti
        end
        @test volume(tg2) == sum(volume(tg2, i) for i in 1:length(tg2))

        MeshGrids.matfreq(tg2) ≈ reverse(DLR.ωn)

        # test floor
        # notice for imfreq tg2[i]!=tg2.grid[i]
        δ = 1e-6
        N = length(tg2)
        @test floor(tg2, tg2.grid[N] - δ) == 1
        @test floor(tg2, tg2.grid[N] + δ) == 1
        for i in 2:N-1
            @test floor(tg2, tg2.grid[N+1-i] - δ) == i
            @test floor(tg2, tg2.grid[N+1-i] + δ) == i - 1
        end
        @test floor(tg2, tg2.grid[1] - δ) == N - 1
        @test floor(tg2, tg2.grid[1] + δ) == N - 1
    end

    @testset "DLRFreq Grid" begin
        tg2 = MeshGrids.DLRFreq(DLR)
        println(tg2)

        @test length(tg2) == length(DLR)
        # @test size(tg2) == size(DLR)
        # TODO: seems size(DLR) = N instead of (N, ), while size(Vector{N})=(N, )
        @test tg2[1] == DLR.ω[1]

        # eltype
        @test eltype(typeof(tg2)) == Float64

        for (ti, t) in enumerate(tg2)
            @test tg2.grid[ti] ≈ DLR.ω[ti]
            @test tg2[ti] ≈ DLR.ω[ti] #DLR.ωn is read from files, cannot exactly match (2n+1)/β exactly
            @test MeshGrids.locate(tg2, t) == ti
        end
        @test volume(tg2) ≈ sum(volume(tg2, i) for i in 1:length(tg2))
    end

    @testset "CompositeGrids and BrillouinZoneMeshes" begin
        # test availability of imported package

        # CompositeGrids
        N1 = 7
        mesh1 = MeshGrids.SimpleGrid.Uniform{Float64}([0.0, 1.0], N1)
        for (i, x) in enumerate(mesh1)
            @test MeshGrids.locate(mesh1, x) == i
        end
        @test MeshGrids.volume(mesh1) ≈ sum(MeshGrids.volume(mesh1, i) for i in 1:length(mesh1))

        # BrillouinZoneMeshes
        lattice = Matrix([2.0 0 0; 1 sqrt(3) 0; 7 11 19]')
        msize = (3, 5, 7)
        br = MeshGrids.BZMeshes.Cell(lattice=lattice)
        mesh2 = MeshGrids.BZMeshes.UniformBZMesh(cell=br, size=msize)
        for (i, x) in enumerate(mesh2)
            @test MeshGrids.locate(mesh2, x) == i
        end
        @test MeshGrids.volume(mesh2) ≈ sum(MeshGrids.volume(mesh2, i) for i in 1:length(mesh2))
    end
end