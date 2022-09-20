@testset "GreenFunc" begin
    @testset "MeshProduct" begin
        mesh1 = SimpleGrid.Uniform{Float64}([0.0, 1.0], 5)
        mesh2 = SimpleGrid.Uniform{Float64}([0.0, 1.0], 5)
        x = 3
        y = 4
        I = 14
        meshprod = MeshProduct(mesh1,mesh2)
        println(typeof(meshprod))
        println("rank of meshprod is $(rank(meshprod))")
        println("size of mesh1 is $(size(meshprod,1)), size of mesh2 is $(size(meshprod,2))")
        println("size of meshprod is $(size(meshprod)), length of meshprod is $(length(meshprod))")
        println("test index_to_linear for index(($x),($y)):$(index_to_linear(meshprod,x,y))")
        println("test linear_to_index for I=$(I):$(index_to_linear(meshprod,I))") 
        println("test getindex with index input:\n",meshprod[x,y])
        println("test getindex with linearindex input:\n",meshprod[I])
        println("Print test:\n",show(meshprod))
    end   
     
end

