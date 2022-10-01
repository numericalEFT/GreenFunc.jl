
module HubbardRPA

using GreenFunc
using GreenFunc.Triqs.PythonCall
using Parameters
using JLD2

export hubbard_rpa, save_hubbard_rpa_list, load_hubbard_rpa_list

#TODO: load rpa kernel of hubbard model from triqs and triqs_tprf

# load python script(import from file not allowed in PythonCall)

f = open("./example/hubbard_response/tprf_rpa.py", "r")
# exec python script
pyexec(read(f, String), Main)
gamma_rpa = pyeval(Py, "gamma_rpa", Main)

@with_kw struct Para
    norb::Int = 1
    t::Float64 = 1.0
    nk::Int = 32
    dim::Int = 2
    beta::Float64 = 10.0
    nw::Int = 100
    mu::Float64 = 0.0
    U::Float64 = 1.0
end

function hubbard_rpa(; norb=1, t=1.0, nk=32, dim=2, beta=10, nw=100, mu=0, U=1.0)
    gamma, green0 = gamma_rpa(norb, t, nk, dim, beta, nw, mu, U)
    return GreenFunc.MeshArray(gamma), GreenFunc.MeshArray(green0)
end

function hubbard_rpa(para::Para)
    @unpack norb, t, nk, dim, beta, nw, mu, U = para
    gamma, green0 = gamma_rpa(norb, t, nk, dim, beta, nw, mu, U)
    return GreenFunc.MeshArray(gamma), GreenFunc.MeshArray(green0)
end

function save_hubbard_rpa_list(;
    para::Para=Para(), betas=[para.beta,],
    fname="./run/hubbard_rpa.jld2")
    paras = [Para(para, beta=beta) for beta in betas]
    funcs = [hubbard_rpa(param) for param in paras]
    greens = [funcs[i][2] for i in 1:length(funcs)]
    gammas = [funcs[i][1] for i in 1:length(funcs)]
    save(fname, Dict(
            "paras" => paras,
            "greens" => greens,
            "gammas" => gammas
        ), compress=true)
    return paras, greens, gammas
end

function load_hubbard_rpa_list(;
    fname="./run/hubbard_rpa.jld2")
    try
        f = load(fname)
        paras = f["paras"]
        greens = f["greens"]
        gammas = f["gammas"]
        return paras, greens, gammas
    catch
        println("Warning[load_hubbard_rpa_list]:" * fname * " fail to open!")
    end
end

end
