@testset "GreenFunc" begin
    @testset "Green2" begin
        tgrid = [0.0,1.0]
        sgrid = [0.0,1.0]
        color_n = [0.0,1.0]
	      green_simple = GreenBasic.Green2{Float64}(:freq,:mom,:fermi,color_n,tgrid,sgrid)
        println(green_simple.value,green_simple[1,1,1])
        green_simple.value[1,1,1] = 1.0
        green_simple.timeType=:time
        println(green_simple.timeType,green_simple[1,1,1])        
    end
end

