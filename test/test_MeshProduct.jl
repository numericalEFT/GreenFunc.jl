@testset "GreenFunc" begin
    @testset "MeshProduct" begin
        mesh = [0.0, 1.0]
        meshprod = MeshProduct(mesh,mesh)
        println(typeof(meshprod))
        println("size of mesh is $(size(mesh)), size of meshprod is $(size(meshprod))")
        println("length of mesh is $(size(mesh)), length of meshprod is $(size(meshprod))")
        value = [0.0,1.0]
    end   
     
end

