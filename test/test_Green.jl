SemiCircle(dlr, grid, type) = Sample.SemiCircle(dlr.Euv, dlr.β, dlr.isFermi, dlr.symmetry, grid, type, dlr.rtol, 24, true)
@testset "GreenFunc" begin
    # @testset "Green2" begin
    #     tgrid = [0.0,1.0]
    #     sgrid = [0.0,1.0]
    #     color_n = [0.0,1.0]
    #     beta = 20.0
	  #     green_simple = GreenBasic.Green2{Float64}(:freq,:mom,true,:ph,nothing,beta,color_n,tgrid,sgrid)
    #     println(green_simple.dynamic)        
    # end
    @testset "Green2DLR" begin
        sgrid = [0.0]
        color_n = [0.0]
        β = 10.0
        isFermi = true
        rtol = 1e-10
        Euv = 1000.0
        symmetry=:none
        dlr = DLRGrid(Euv, β, rtol, isFermi, symmetry)
        tgrid = dlr.n
        green_freq = GreenBasic.Green2DLR{ComplexF64}(true,Euv,rtol,:k,sgrid, β, :n, timeSymmetry=symmetry)
        Gτ = SemiCircle(dlr, dlr.τ, :τ)
        Gn = SemiCircle(dlr, dlr.n, :ωn)

        for (ti,t) in enumerate(tgrid)
            for (qi, q) in enumerate(sgrid)
                for (c1i,c1) in enumerate(color_n)
                    for (c2i,c2) in enumerate(color_n)
                        green_freq.dynamic[c1i,c2i,qi,ti] = Gn[ti]
                    end
                end
            end
        end

        green_tau = GreenBasic.toTau(green_freq)
        err =  maximum(abs.(green_tau.dynamic[1,1,1,:].-Gτ))
        printstyled("SemiCircle Fourier ωn->τ $err\n",color = :white)
        @test err<50*rtol

        green_freq_compare = GreenBasic.toMatFreq(green_tau)
        err = maximum(abs.(green_freq_compare.dynamic[1,1,1,:].-Gn))
        printstyled("SemiCircle Fourier τ->ωn $err\n",color = :white)        
        @test err<50*rtol

        green_dlr = GreenBasic.toDLR(green_freq)
        green_tau = GreenBasic.toTau(green_dlr)
        err =  maximum(abs.(green_tau.dynamic[1,1,1,:].-Gτ))
        printstyled("SemiCircle Fourier ωn->dlr->τ $err\n",color = :white)
        @test err<50*rtol
    end
end

