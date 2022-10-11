
include("./HubbardRPA.jl")
using .HubbardRPA
using .HubbardRPA: GreenFunc
using .GreenFunc
using GreenFunc: CompositeGrids

if !isempty(ARGS)
    fname = ARGS[1]
    pid = parse(Int, ARGS[2])
else
    fname = pwd() * "/run/hubbard_rpa.jld2"
    pid = 1
end

try
    global paras, greens, gammas = load_hubbard_rpa_list(fname=fname)
catch
    println("Warning:" * fname * " does not exist!")
    println("Generating new " * fname)
    betas = [10,]
    global paras, greens, gammas = save_hubbard_rpa_list(betas=betas, fname=fname)
end

gr_w = greens[pid]
ga_w = gammas[pid]
para = paras[pid]


gr_dlr = gr_w |> to_dlr
gr_imt = dlr_to_imtime(gr_dlr, CompositeGrids.SimpleG.Uniform{Float64}([0, para.beta], para.nw))

ga0 = ga_w[1, 1, 1, 1, 1, 1]
ga_w .-= ga0
ga_dlr = ga_w |> to_dlr
ga_imt = dlr_to_imtime(ga_dlr, CompositeGrids.SimpleG.Uniform{Float64}([0, para.beta], para.nw))
# use default here cause ERROR:
# when dlr.τ is converted to CompositeGrid.SimpleG.Arbitrary, the bound is not [0, β]

rdyn = MeshArray(gr_imt.mesh[3:4]...)
rdyn.data .= 0.0
r0 = MeshArray(gr_imt.mesh[3])
r0.data .= 0.0

function green(k, t1, t2)
    t = t2 - t1
    factor = 1.0
    if t < 0
        t = t + para.beta
        factor = -1.0
    end

    ik, it = locate(gr_imt.mesh[3], k), locate(gr_imt.mesh[4], t)

    return real(gr_imt[1, 1, ik, it])
end

function bareW0(k)
    return real(ga0)
end

function dynamicW0(k, t1, t2)
    t = t2 - t1
    if t < 0
        t = t + para.beta
    end

    ik, it = locate(ga_imt.mesh[5], k), locate(ga_imt.mesh[6], t)

    return real(ga_imt[1, 1, 1, 1, ik, it])
end

function bareR(k; data=r0.data)
    ik = locate(r0.mesh[1], k)
    return real(data[ik])
end

function dynamicR(k, t1, t2; data=rdyn.data)
    t = t2 - t1
    factor = 1.0
    if t < 0
        t = t + para.beta
        factor = -1.0
    end

    ik, it = locate(rdyn.mesh[1], k), locate(rdyn.mesh[2], t)

    return real(data[ik, it])
end
