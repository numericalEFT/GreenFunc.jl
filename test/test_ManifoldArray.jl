SemiCircle(dlr, grid, type) = Sample.SemiCircle(dlr.Euv, dlr.β, dlr.isFermi, grid, type, dlr.symmetry; rtol=dlr.rtol, degree=24, regularized=true)

@testset "ManifoldArray" begin
    function test_shape(N1, N2, innermesh)
        ############# basic test ################
        mesh1 = SimpleGrid.Uniform{Float64}([0.0, 1.0], N1)
        mesh2 = SimpleGrid.Uniform{Float64}([0.0, 1.0], N2)
        if isempty(innermesh)
            g = ManifoldArray(mesh1, mesh2)
            @test length(g) == N1 * N2
        else
            g = ManifoldArray(innermesh..., mesh1, mesh2)
            @test length(g) == N1 * N2 * reduce(*, length.(innermesh))
        end
        show(g)
        @test size(g) == (length.(innermesh)..., N1, N2)
        @test eltype(typeof(g)) == Float64
        gc = similar(g, ComplexF64)
        @test eltype(typeof(gc)) == ComplexF64

        ############ broadcast test ###################
        g.data = rand(g.dims...)
        if isempty(innermesh)
            g2 = ManifoldArray(mesh1, mesh2; data=rand(g.dims...))
        else
            g2 = ManifoldArray(innermesh..., mesh1, mesh2; data=rand(g.dims...))
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
    function test_fourier(N1,beta, statistics, innermesh)
        mesh1 = SimpleGrid.Uniform{Float64}([0.0, 1.0], N1)
        mesh2 = MeshGrids.DLRFreq(beta,statistics)
        if isempty(innermesh)
            g = ManifoldArray(mesh1, mesh2)
        else
            g = ManifoldArray(innermesh..., mesh1, mesh2)
        end
        g_freq = dlr_to_imfreq(g)
        Gτ = SemiCircle(mesh2.dlr, mesh2.dlr.τ, :τ)
        Gn = SemiCircle(mesh2.dlr, mesh2.dlr.n, :n)
        for (ni,n) in enumerate(mesh2.dlr.n)
            g_freq.data[:,ni] .= Gn[ni]
        end
        g_dlr=to_dlr(g_freq)
        rtol = mesh2.dlr.rtol

        g_time = dlr_to_imtime(g_dlr)
        err = maximum(abs.(g_time.data[1,:] .- Gτ))
        printstyled("test dlr_to_imtime dlr->τ $err\n", color=:white)
        
        @test err < 50 * rtol
        g_freq1 = dlr_to_imfreq(g_dlr)
        err = maximum(abs.(g_freq1.data[1,:] .- Gn))
        printstyled("test dlr_to_imfreq $err\n", color=:white)
        @test err < 50 * rtol
        g_freq1<<g_dlr
        err = maximum(abs.(g_freq1.data[1,:] .- Gn))
        printstyled("test  imfreq<<dlr $err\n", color=:white)
        @test err < 50 * rtol
        g_time <<g_dlr
        err = maximum(abs.(g_time.data[1,:] .- Gτ))
        printstyled("test  imtime<<dlr $err\n", color=:white)
        
        @test err < 50 * rtol
        
    end
    test_fourier(5,100.0,FERMI, ())
end
