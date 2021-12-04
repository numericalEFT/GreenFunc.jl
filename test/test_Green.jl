@testset "GreenFunc" begin
    @testset "GreenUR" begin
        init = zeros(Float64,2,2,2)
	      green_simple = Greenfunc.GreenUR{Float64}(:freq,:mom,2,[0.0,1.0],[0.0,1.0],init)
        println(green_simple.value,green_simple[1,2,3])
    end
end
n
