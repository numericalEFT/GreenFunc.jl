SemiCircle(dlr, grid, type) = Sample.SemiCircle(dlr.Euv, dlr.β, dlr.isFermi, grid, type, dlr.symmetry; rtol=dlr.rtol, degree=24, regularized=true)

@testset "GreenFunc" begin
    @testset "GreenDLR" begin
        mesh = [0.0, 1.0]
        β = 10.0
        isFermi = true
        Euv = 80.0
        rtol = 1e-9
        innerstate = (1,)
        isFermi = true
        tsym = :ph
        # data = [1.0im, 1.0im]
        data = zeros(ComplexF64, (1, 2, 3))
        tgrid = [1, 2, 3]
        DLR = DLRGrid(Euv, β, rtol, isFermi, tsym)
        # green_freq = GreenDLR{ComplexF64}(; domain=GreenFunc.IMFREQ, DLR=DLR, tgrid=tgrid, mesh=mesh, β=β, isFermi=isFermi, Euv=Euv, rtol=rtol, tsym=tsym, innerstate=innerstate, data=data)
        green_freq = GreenDLR{ComplexF64}(; domain=GreenFunc.IMFREQ, DLR=DLR, tgrid=tgrid, mesh=mesh, innerstate=innerstate, data=data)
        # green_freq = GreenDLR(; mesh=mesh, β=β, tsym=tsym, Euv=Euv, rtol=rtol)
        println(typeof(green_freq))
        show(green_freq)

        println("size of green_freq is $(size(green_freq))\n")
        @test size(green_freq) == (1, 2, length(green_freq.tgrid))
        println("test getindex (1, 2, 3): ", green_freq[1, 2, 3])
        @test green_freq[1, 2, 3] == 0
        println("test view:\n", green_freq[:, 1:2, 1:3])
        green_freq[1, 1, 1] = im * 2.0
        println("test setindex:", green_freq[1, 1, 1])
        @test green_freq[1, 1, 1] == 2.0im
        rank_var = GreenFunc.rank(green_freq)
        println("rank of green_freq is $rank_var")
        @test rank_var == 3
        # println("density matrix: ", GreenFunc.density(green_freq))  # Wait for toTau()

        green_freq2 = deepcopy(green_freq)
        g1, g2, g3, g4 = (-green_freq)[1, 1, 1], (green_freq+green_freq2)[1, 1, 1], (green_freq-green_freq2)[1, 1, 1], (green_freq*green_freq2)[1, 1, 1]
        println("test math: $g1, $g2, $g3, $g4")
        @test green_freq[1, 1, 1] == 2.0im && green_freq2[1, 1, 1] == 2.0im
        @test g1 == -2.0im
        @test g2 == 4.0im
        @test g3 == 0
        @test g4 == -4

        green3 = similar(green_freq)
        println("\nsimialr green:")
        show(green3)
        println("view similar green:\n", green3[:, 1:2, 1:3])
        green3[1, 2, 3] = -1 - 2.0im
        green_freq[1, 2, 3] = 1 + 2.0im
        println("getindex (1, 2, 3): similar green $(green3[1, 2, 3]), original green $(green_freq[1, 2, 3])")
        @test green3[1, 2, 3] == -1 - 2.0im && green_freq[1, 2, 3] == 1 + 2.0im

        # testing iteration
        green4 = similar(green_freq)

        for (id, d) in enumerate(green4)
            inds = GreenFunc.ind2sub_gen(size(green4), id)
            # println("$id -> $inds")
            @test d == green4[id]
            @test d == green4[inds...]
        end
        green4 << :(1 / (ωn^2 + p^4))
        for (id, d) in enumerate(green4)
            inds = GreenFunc.ind2sub_gen(size(green4), id)
            @test d == green4[id]
            @test d == green4[inds...]
        end

        #     Gτ = SemiCircle(green_freq.DLR, green_freq.DLR.τ, :τ)
        #     Gn = SemiCircle(green_freq.DLR, green_freq.DLR.n, :n)
        #     green_dum = zeros(ComplexF64, (green_freq.color, green_freq.color, green_freq.spaceGrid.size, green_freq.timeGrid.size))
        #     for (ti, t) in enumerate(green_freq.timeGrid)
        #         for (qi, q) in enumerate(green_freq.spaceGrid)
        #             for (c1i, c1) in enumerate(color_n)
        #                 for (c2i, c2) in enumerate(color_n)
        #                     green_dum[c1i, c2i, qi, ti] = Gn[ti]
        #                 end
        #             end
        #         end
        #     end
    end

    #     green_freq.dynamic = green_dum
    #     green_tau = toTau(green_freq)
    #     err = maximum(abs.(green_tau.dynamic[1, 1, 1, :] .- Gτ))
    #     printstyled("SemiCircle Fourier ωn->τ $err\n", color = :white)
    #     @test err < 50 * rtol

    #     green_freq_compare = toMatFreq(green_tau)
    #     err = maximum(abs.(green_freq_compare.dynamic[1, 1, 1, :] .- Gn))
    #     printstyled("SemiCircle Fourier τ->ωn $err\n", color = :white)
    #     @test err < 50 * rtol

    #     green_dlr = toDLR(green_freq)
    #     green_tau = toTau(green_dlr)
    #     err = maximum(abs.(green_tau.dynamic[1, 1, 1, :] .- Gτ))
    #     printstyled("SemiCircle Fourier ωn->dlr->τ $err\n", color = :white)
    #     @test err < 50 * rtol

    #     #test JLD2
    #     ############# FileIO API #################
    #     save("example.jld2", Dict("green" => green_freq), compress = true)
    #     d = load("example.jld2")
    #     green_read = d["green"]
    #     @test green_read.dynamic == green_freq.dynamic
    #     #deeptest(green_read, green_freq)
    #     rm("example.jld2")
    # end

    # @testset "find" begin
    #     sgrid = [0.0, 1.0]
    #     color_n = [0.0]
    #     β = 10.0
    #     isFermi = true
    #     Euv = 1000.0

    #     green_linear = Green2DLR{Float64}(:green, GreenFunc.IMTIME, β, isFermi, Euv, sgrid)
    #     rtol = green_linear.DLR.rtol
    #     green_dum = zeros(Float64, (green_linear.color, green_linear.color, green_linear.spaceGrid.size, green_linear.timeGrid.size))
    #     for (ti, t) in enumerate(green_linear.timeGrid)
    #         for (qi, q) in enumerate(green_linear.spaceGrid)
    #             for (c1i, c1) in enumerate(color_n)
    #                 for (c2i, c2) in enumerate(color_n)
    #                     green_dum[c1i, c2i, qi, ti] = t * q
    #                 end
    #             end
    #         end
    #     end

    #     green_linear.dynamic = green_dum
    #     green_dum_ins = zeros(Float64, (green_linear.color, green_linear.color, green_linear.spaceGrid.size))
    #     for (qi, q) in enumerate(green_linear.spaceGrid)
    #         for (c1i, c1) in enumerate(color_n)
    #             for (c2i, c2) in enumerate(color_n)
    #                 green_dum_ins[c1i, c2i, qi] = q
    #             end
    #         end
    #     end

    #     green_linear.instant = green_dum_ins

    #     τ = 0.5
    #     x = 0.3
    #     interp_dym = dynamic(green_linear, τ, x, 1, 1)
    #     @test interp_dym - τ * x < 1e-8
    #     interp_ins = instant(green_linear, x, 1, 1)
    #     @test interp_ins - x < 1e-8
    #     interp_ins = dynamic(green_linear, τ, x, 1, 1, GreenFunc.DEFAULTINTERP, GreenFunc.DEFAULTINTERP)
    #     @test interp_ins - x < 1e-8
    #     interp_ins = dynamic(green_linear, τ, x, 1, 1, GreenFunc.DLRINTERP, GreenFunc.DEFAULTINTERP)
    #     @test interp_ins - x < 1e-8
    # end
end

