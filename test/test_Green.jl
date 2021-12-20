@testset "GreenFunc" begin
    @testset "Green2" begin
        tgrid = [0.0,1.0]
        sgrid = [0.0,1.0]
        color_n = [0.0,1.0]
        beta = 20.0
	      green_simple = GreenBasic.Green2{Float64}(:freq,:mom,true,:particlehole,:none,beta,color_n,tgrid,sgrid)
        println(green_simple.dynamic)        
    end
end

