using GreenFunc
using Test
using CondaPkg
using PythonCall

@testset "HubbardRPA" begin
    include("HubbardRPA.jl")

    @testset "tprf_rpa.py" begin
        f = open("tprf_rpa.py", "r")
        pyexec(read(f, String), Main)
        gamma_rpa = pyeval(Py, "gamma_rpa", Main)
        gamma = gamma_rpa(norb=1, t=1.0, nk=32, dim=2, beta=10, n_max=100, mu=0, U=1.0)
        jgamma = GreenFunc.MeshArray(gamma)
    end

    @testset "HubbardRPA.jl" begin
        jgamma = Gamma(norb=1, t=1.0, nk=32, dim=2, beta=10, n_max=100, mu=0, U=1.0)
    end
end
