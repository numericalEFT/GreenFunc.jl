using PythonCall

@testset "Triqs interface" begin
    gf = pyimport("triqs.gf")
    np = pyimport("numpy")

    ############ test imaginary-time mesh #############
    mt = gf.MeshImTime(beta=1.0, S="Fermion", n_max=100)
    l = @py len(mt)
    lj = pyconvert(Int, l)

    ############ test MeshArray constructor #############
    G_t = gf.GfImTime(mesh=mt, data=np.random.rand(lj, 2, 3)) #target_shape = [2, 3] --> innerstate = [3, 2]
    gt = MeshArray(G_t)
    @test size(gt) == (3, 2, lj)
    i1, i2, t = 1, 2, 3
    @test gt[i1, i2, t] ≈ pyconvert(Float64, G_t.data[t-1, i2-1, i1-1])
    ######### test << ###############################
    G_t = gf.GfImTime(mesh=mt, data=np.random.rand(lj, 2, 3)) #target_shape = [2, 3] --> innerstate = [3, 2]
    gt << G_t
    @test size(gt) == (3, 2, lj)
    i1, i2, t = 1, 2, 3
    @test gt[i1, i2, t] ≈ pyconvert(Float64, G_t.data[t-1, i2-1, i1-1])

    ############ test Matsubara frequency mesh #############
    miw = gf.MeshImFreq(beta=1.0, S="Fermion", n_max=100)
    l = @py len(miw)
    lj = pyconvert(Int, l)

    ############ test MeshArray constructor #########
    G_w = gf.GfImFreq(mesh=miw, data=np.random.rand(lj, 2, 3)) #target_shape = [2, 3] --> innerstate = [3, 2]
    gw = MeshArray(G_w)
    @test size(gw) == (3, 2, lj)
    i1, i2, t = 1, 2, 3
    @test gw[i1, i2, t] ≈ pyconvert(Float64, G_w.data[t-1, i2-1, i1-1])
    ########### test << ###########################
    G_w = gf.GfImFreq(mesh=miw, data=np.random.rand(lj, 2, 3)) #target_shape = [2, 3] --> innerstate = [3, 2]
    gw << G_w
    @test size(gw) == (3, 2, lj)
    i1, i2, t = 1, 2, 3
    @test gw[i1, i2, t] ≈ pyconvert(Float64, G_w.data[t-1, i2-1, i1-1])

    ########## test MeshBrZone #######################
    lat = pyimport("triqs.lattice")
    BL = lat.BravaisLattice(units=((2, 0, 0), (1, sqrt(3), 0))) # testing with a triangular lattice so that exchanged index makes a difference
    BZ = lat.BrillouinZone(BL)
    nk = 8
    mk = gf.MeshBrillouinZone(BZ, nk)
    mprod = gf.MeshProduct(mk, miw)
    G_k_w = gf.GfImFreq(mesh=mprod, target_shape=[1, 1]) #G_k_w.data.shape will be [nk^2, lj, 1, 1]
    gkw = MeshArray(G_k_w)
    umesh = gkw.mesh[4]
    for p in mk
        ilin = pyconvert(Int, p.linear_index) + 1
        inds = pyconvert(Array, p.index)[1:2] .+ 1
        pval = pyconvert(Array, p.value)
        # note that while linear_index is kept, the index reversed
        @test pval[1:2] ≈ umesh[ilin]
        @test pval[1:2] ≈ umesh[reverse(inds)...]
    end
    @test size(gkw) == (1, 1, lj, nk^2)
    ik, iw = 12, 10
    @test gkw[1, 1, iw, ik] ≈ pyconvert(Float64, G_k_w.data[ik-1, iw-1, 0, 0])
    ########### test << ###########################
    G_k_w = gf.GfImFreq(mesh=mprod, data=np.random.rand(nk^2, lj, 1, 1)) #G_k_w.data.shape will be [nk^2, lj, 1, 1]
    gkw << G_k_w
    ik, iw = 12, 10
    @test gkw[1, 1, iw, ik] ≈ pyconvert(Float64, G_k_w.data[ik-1, iw-1, 0, 0])

end