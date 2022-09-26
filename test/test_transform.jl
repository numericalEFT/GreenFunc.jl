SemiCircle(dlr, grid, type) = Sample.SemiCircle(dlr.Euv, dlr.β, dlr.isFermi, grid, type, dlr.symmetry; rtol=dlr.rtol, degree=24, regularized=true)

@testset "Transform" begin
    function test_fourier(N1, beta, statistics)
        mesh1 = SimpleGrid.Uniform{Float64}([0.0, 1.0], N1)
        mesh2 = MeshGrids.DLRFreq(beta, statistics)
        g = MeshArray(mesh1, mesh2; data=zeros(N1, length(mesh2)))

        g_freq = dlr_to_imfreq(g)
        Gτ = SemiCircle(mesh2.dlr, mesh2.dlr.τ, :τ)
        Gn = SemiCircle(mesh2.dlr, mesh2.dlr.n, :n)
        # iterate over the last dimension of g_freq.data
        for (ni, n) in enumerate(mesh2.dlr.n)
            g_freq.data[:, ni] .= Gn[ni]
        end
        g_dlr = imfreq_to_dlr(g_freq)
        rtol = mesh2.dlr.rtol

        g_time = dlr_to_imtime(g_dlr)
        err = maximum(abs.(g_time.data[1, :] .- Gτ))
        printstyled("test dlr_to_imtime dlr->τ $err\n", color=:white)

        @test err < 50 * rtol
        g_freq1 = dlr_to_imfreq(g_dlr)
        err = maximum(abs.(g_freq1.data[1, :] .- Gn))
        printstyled("test dlr_to_imfreq $err\n", color=:white)
        @test err < 50 * rtol

        ########### test pipe operation #############
        g_freq2 = g_freq |> to_dlr |> dlr_to_imfreq
        err = maximum(abs.(g_freq2.data[1, :] .- Gn))
        printstyled("test dlr_to_imfreq with pipe $err\n", color=:white)
        @test err < 50 * rtol

        g_time2 = g_freq |> to_dlr |> dlr_to_imtime
        err = maximum(abs.(g_time2.data[1, :] .- Gτ))
        printstyled("test dlr_to_imtime with pipe $err\n", color=:white)
        @test err < 50 * rtol


        g_freq1 << g_dlr
        err = maximum(abs.(g_freq1.data[1, :] .- Gn))
        printstyled("test  imfreq<<dlr $err\n", color=:white)
        @test err < 50 * rtol
        g_time << g_dlr
        err = maximum(abs.(g_time.data[1, :] .- Gτ))
        printstyled("test  imtime<<dlr $err\n", color=:white)

        @test err < 50 * rtol


    end
    test_fourier(5, 100.0, FERMI)
    test_fourier(5, 100.0, BOSE)
end