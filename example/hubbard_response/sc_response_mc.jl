
using MCIntegration
using CompositeGrids
using Measurements

include("propagator.jl")

function integrand(vars, config)
    ### unpack variables ####
    Kxy, T, ExtKidx, ExtTidx = vars
    Kx, Ky = Kxy #unpack CompositeVar
    beta, extKgrid, extTgrid = config.userdata

    extK = extKgrid[ExtKidx[1]] #set external kx, ky
    k = (Kx[1], Ky[1]) #set internal kx, ky
    q = (k[1] - extK[1], k[2] - extK[2]) #set q = k - k'

    t1, t2 = 0.0, extTgrid[ExtTidx[1]] #set external t1, t2
    t3, t4 = T[3], T[4]  # set internal t3, t4

    #### instant R = instant W0 * G * G * (instant R + dynamic R) ############
    instant = green(k, t1, t3) * dynamicR(k, t3, t4) * green(k, t4, t1)
    instant += green(k, t1, t3) * bareR(k) * green(k, t3, t1) / beta
    # both t3 and t4 will be integrated over, but bareR doesn't depend on t3, t4. So one must divide by beta
    instant *= bareW0(q)

    #### dynamic R = dynamic W0 * G * G * (instant R + dynamic R) ############
    dynamic = green(k, t1, t3) * dynamicR(k, t3, t4) * green(k, t4, t2)
    dynamic += green(k, t1, t3) * bareR(k) * green(k, t3, t2) / beta
    dynamic *= dynamicW0(q, t1, t2)

    return instant, dynamic
end

# for vegas algorithm
function measure(vars, obs, weights, config)
    Kxy, T, ExtKidx, ExtTidx = vars
    obs[1][ExtKidx[1]] += weights[1]
    obs[2][ExtKidx[1], ExtTidx[1]] += weights[2]
end

function PPver(;
    neval=1e6, #number of evaluations
    print=-1,
    alpha=3.0, #learning ratio
    config=nothing,
    para::HubbardRPA.Para=para,
    kwargs...)
    HubbardRPA.@unpack t, nk, dim, beta, U = para

    Euv = 4t

    kmesh = ga_w.mesh[5]

    ##### prepare external K grid and tau grid ##############
    extKgrid = [(k[1], k[2]) for k in kmesh]
    extTgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, beta], [0.0, beta], 8, 1 / Euv, 8) #roughly ~100 points if resolution = β/128

    nt = length(extTgrid)

    latvec = 1.0 * π

    Kx = Continuous(-latvec, latvec; alpha=alpha)
    Ky = Continuous(-latvec, latvec; alpha=alpha)
    K = CompositeVar(Kx, Ky)
    T = Continuous(0.0, beta; offset=2, alpha=alpha) #the first tau is fixed to 0.0, the second tau is reserved for the external time
    T.data[1] = 0.0
    ExtKxyidx = Discrete(1, nk * nk, alpha=alpha)
    ExtTidx = Discrete(1, nt, alpha=alpha)

    dof = [[1, 2, 1, 0], [1, 2, 1, 1]] # K, T, ExtKidx, ExtTidx
    # there are only 2 time variables (T[3], T[4]), while T[1] and T[2] are fixed to 0.0 and the external time

    obs = [zeros(Float64, nk * nk), zeros(Float64, nk * nk, nt)] # instant R and dynamic R

    if isnothing(config)
        config = Configuration(;
            # var=(R, Theta, Phi, T, X, ExtKidx),
            var=(K, T, ExtKxyidx, ExtTidx),
            dof=dof,
            type=Float64, # type of the integrand
            obs=obs,
            userdata=(beta, extKgrid, extTgrid),
            kwargs...
        )
    end

    result = integrate(integrand; config=config, measure=measure, print=print, neval=neval, kwargs...)
    # # result = integrate(integrandKW; config=config, print=print, neval=neval, kwargs...)
    # # niter=niter, print=print, block=block, kwargs...)

    if isnothing(result) == false

        if print >= -1
            report(result.config)
            println(report(result, pick=o -> first(o)))
            println(result)
        end
        if print >= -2
            println(result)
        end

        datadict = Dict{Int,Any}()

        for o in 1:length(dof)
            avg, std = result.mean[o], result.stdev[o]
            r = measurement.(real(avg), real(std))
            i = measurement.(imag(avg), imag(std))
            data = Complex.(r, i)
            # datadict[partition[o]] = data
            datadict[o] = data
        end
        return datadict, result
    else
        return nothing, nothing
    end
end


# solve linear response function and compute Tc 

PPver(neval=1e6)
