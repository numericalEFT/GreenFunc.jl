using PythonCall

@testset "Triqs interface" begin
    gf = pyimport("triqs.gf")
    np = pyimport("numpy")

    ############ test Matsubara frequency mesh #############
    miw = gf.MeshImFreq(beta=1.0, S="Fermion", n_max=100)
    l = @py len(miw)
    lj = pyconvert(Int, l)

    ############ test GreenNew constructor #########
    G_w = gf.GfImFreq(mesh=miw, data=np.random.rand(lj, 2, 3)) #target_shape = [2, 3] --> innerstate = [3, 2]
    gw = GreenNew(G_w)
    @test size(gw) == (3, 2, lj)
    i1, i2, t = 1, 2, 3
    @test gw[i1, i2, t] ≈ pyconvert(Float64, G_w.data[t-1, i2-1, i1-1])
    ########### test << ###########################
    G_w = gf.GfImFreq(mesh=miw, data=np.random.rand(lj, 2, 3)) #target_shape = [2, 3] --> innerstate = [3, 2]
    gw << G_w
    @test size(gw) == (3, 2, lj)
    i1, i2, t = 1, 2, 3
    @test gw[i1, i2, t] ≈ pyconvert(Float64, G_w.data[t-1, i2-1, i1-1])

    ############ test imaginary-time mesh #############
    mt = gf.MeshImTime(beta=1.0, S="Fermion", n_max=100)
    l = @py len(mt)
    lj = pyconvert(Int, l)

    ############ test GreenNew constructor #############
    G_t = gf.GfImTime(mesh=mt, data=np.random.rand(lj, 2, 3)) #target_shape = [2, 3] --> innerstate = [3, 2]
    gt = GreenNew(G_t)
    @test size(gt) == (3, 2, lj)
    i1, i2, t = 1, 2, 3
    @test gt[i1, i2, t] ≈ pyconvert(Float64, G_t.data[t-1, i2-1, i1-1])
    ######### test << ###############################
    G_t = gf.GfImTime(mesh=mt, data=np.random.rand(lj, 2, 3)) #target_shape = [2, 3] --> innerstate = [3, 2]
    gt << G_t
    @test size(gt) == (3, 2, lj)
    i1, i2, t = 1, 2, 3
    @test gt[i1, i2, t] ≈ pyconvert(Float64, G_t.data[t-1, i2-1, i1-1])

end