using GreenFunc
using Test
using CondaPkg
using PythonCall

@testset "HubbardRPA" begin
    include("./HubbardRPA.jl")
    using .HubbardRPA

    @testset "tprf_rpa.py" begin
        f = open("./example/hubbard_response/tprf_rpa.py", "r")
        pyexec(read(f, String), Main)
        gamma_rpa = pyeval(Py, "gamma_rpa", Main)
        gamma, green0 = gamma_rpa(norb=1, t=1.0, nk=32, dim=2, beta=10, n_max=100, mu=0, U=1.0)
        jgamma = GreenFunc.MeshArray(gamma)
    end

    @testset "HubbardRPA.jl" begin
        para = HubbardRPA.Para()
        jgamma, jgreen0 = hubbard_rpa(para)
    end

    @testset "save and load" begin
        save_hubbard_rpa_list()
        paras, greens, gammas = load_hubbard_rpa_list()
        print(paras)
        print(greens)
        print(gammas)
    end
end
