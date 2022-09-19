@testset "GreenFunc" begin
    @testset "MeshProduct" begin
        mesh1 = [0.0, 1.0]
        mesh2 = [0.0,1.0,2.0]
        meshprod = MeshProduct(mesh1,mesh2)
        println(typeof(meshprod))
        println("size of mesh is $(size(mesh)), size of meshprod is $(size(meshprod))")
        println("length of mesh is $(size(mesh)), length of meshprod is $(size(meshprod))")
        value = [0.0,1.0]
    end   
     
end

