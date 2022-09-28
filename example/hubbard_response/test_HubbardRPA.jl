using GreenFunc
using Test
using CondaPkg
using PythonCall

@testset "HubbardRPA" begin

    @testset "Test PythonCall and tprf" begin
        # run tprf official example to make sure it works
        tb = pyimport("triqs_tprf.tight_binding")
        gf = pyimport("triqs.gf")
        lat = pyimport("triqs_tprf.lattice")
        sq_lat = tb.create_square_lattice(norb=1, t=1.0)

        nk = 32
        e_k = sq_lat.get_kmesh((nk, nk, 1))

        wmesh = gf.MeshImFreq(beta=10, S="Fermion", n_max=100)
        g0_wk = lat.lattice_dyson_g0_wk(mu=0, e_k=e_k, mesh=wmesh)
    end
end
