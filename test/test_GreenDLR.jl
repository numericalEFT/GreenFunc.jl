SemiCircle(dlr, grid, type) = Sample.SemiCircle(dlr.Euv, dlr.β, dlr.isFermi, grid, type, dlr.symmetry; rtol = dlr.rtol, degree = 24, regularized = true)

@testset "GreenFunc" begin
    @testset "GreenDLR" begin
        mesh = [0.0, 1.0]
        β = 10.0
        isFermi = true
        Euv = 100.0
        rtol = 1e-10
        innerstate = (1,)
        #data = zeros()
        isFermi = true
        tsym = :none
        #data =  Array{ComplexF64, 2}[]
        data = [1.0im,1.0im]
        tgrid = [1,2]
        DLR = DLRGrid(Euv, β, rtol, isFermi, tsym)
        #green_freq = GreenDLR{ComplexF64}(; domain=GreenFunc.IMFREQ, DLR = DLR, tgrid = tgrid, mesh = mesh, β=β , isFermi=isFermi, Euv=Euv, rtol=rtol, tsym=tsym, innerstate=innerstate, data = data)
        green_freq = GreenDLR(; mesh=mesh)
        rtol = green_freq.DLR.rtol
        println(typeof(green_freq))
        println("size of green_freq is $(size(green_freq.data))\n")
        println("test getindex:\n",green_freq[1,2,3])
        println("test view:\n",green_freq[:,1:2,1:3])
        println("size of green is", size(green_freq))
        green_freq[1,1,1] = im*2.0
        println("test setindex:",green_freq[1,1,1])
        green_freq2 = deepcopy(green_freq)
        println("test math:$((green_freq+green_freq2)[1,1,1]), $((green_freq-green_freq2)[1,1,1]), $((green_freq*green_freq2)[1,1,1])")
        rank_var = GreenFunc.rank(green_freq)
        println("rank of green_freq is $rank_var")
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

