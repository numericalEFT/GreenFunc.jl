module HubbardRPA

using GreenFunc
using GreenFunc.Triqs.PythonCall
using Parameters
using JLD2

export hubbard_rpa, save_hubbard_rpa_list

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
    paras = []
    gammas = []
    greens = []
    for beta in betas
        parab = Para(para, beta=beta)
        push!(paras, parab)
        gamma, green = hubbard_rpa(para)
        push!(gammas, gamma)
        push!(greens, green)
    end
    save(fname, Dict(
            "paras" => paras,
            "greens" => greens,
            "gammas" => gammas
        ), compress=true)
end

function load_hubbard_rpa_list(;
    fname="./run/hubbard_rpa.jld2")
    f = load(fname)
    paras = f["paras"]
    greens = f["greens"]
    gammas = f["gammas"]
    return paras, greens, gammas
end

end