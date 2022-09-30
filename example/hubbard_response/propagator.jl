using GreenFunc

include("./HubbardRPA.jl")
using .HubbardRPA


if !isempty(ARGS)
    fname = ARGS[1]
else
    fname = pwd() * "/run/hubbard_rpa.jld2"
end

paras, greens, gammas = load_hubbard_rpa_list(fname=fname)

# try
#     global paras, greens, gammas = load_hubbard_rpa_list(fname=fname)
# catch
#     println("Warning:" * fname * " does not exist!")
#     println("Generating new " * fname)
#     betas = [10, 20, 40, 80]
#     global paras, greens, gammas = save_hubbard_rpa_list(betas=betas, fname=fname)
# end

# print(paras)
# print(greens)
# print(gammas)


function green(k, t1, t2)
    return 1.0
end

function bareW0(k)
    return 1.0
end

function dynamicW0(k, t1, t2)
    return 1.0
end

function bareR(k)
    return 1.0
end

function dynamicR(k, t1, t2)
    return 1.0
end
