"""
A SYK model solver based on a forward fixed-point iteration method.

 The self-energy of the SYK model is given by,

    Σ(τ) = J² * G(τ) * G(τ) * G(β-τ),
    
 where Green's function of the SYK model is given by the Dyson equation,

    G(iωₙ) = -1/(iωₙ -μ + Σ(iωₙ))

 We solve the Dyson equation self-consistently by a weighted fixed point iteration, 
 with weight `mix` assigned to the new iterate and weight `1-mix` assigned to the previous iterate. 

 The self-energy is evaluated in the imaginary time domain, 
 and the Dyson equation is solved in the Matsubara frequency domain.

 The SYK Green's function has particle-hole symmetry when μ=0. 
 You may enforce such symmetry by setting `symmetry = :ph` when initialize the DLR grids.
 A symmetrized solver tends to be more robust than a unsymmetrized one.
"""

using GreenFunc
using Printf

diff(a, b) = maximum(abs.(a - b)) # return the maximum deviation between a and b

conformal_tau(τ, β) = π^(1 / 4) / sqrt(2β) * 1 / sqrt(sin(π * τ / β))

function syk_sigma_dlr(d, G_x, J = 1.0; sumrule = nothing, verbose = false)

    tau_k = d.τ # DLR imaginary time nodes
    tau_k_rev = d.β .- tau_k # Reversed imaginary time nodes

    G_x_rev = tau2tau(d, G_x, tau_k_rev, sumrule = sumrule, verbose = verbose) # G at beta - tau_k

    Sigma_x = J .^ 2 .* G_x .^ 2 .* G_x_rev # SYK self-energy in imaginary time

    return Sigma_x
end

function dyson(d, sigma_q, mu)
    if d.symmetry == :ph #symmetrized G
        @assert mu ≈ 0.0 "Only the case μ=0 enjoys the particle-hole symmetry."
        return 1im * imag.(-1 ./ (d.ωn * 1im .- mu .+ sigma_q))
    elseif d.symmetry == :none
        return -1 ./ (d.ωn * 1im .- mu .+ sigma_q)
    else
        error("Not implemented!")
    end
end

function solve_syk_with_fixpoint_iter(d, mu, tol = d.rtol * 10; mix = 0.1, maxiter = 5000, G_x = zeros(ComplexF64, length(d)), sumrule = nothing, verbose = true)

    for iter in 1:maxiter

        Sigma_x = syk_sigma_dlr(d, G_x, sumrule = sumrule, verbose = verbose)

        G_q_new = dyson(d, tau2matfreq(d, Sigma_x), mu)

        G_x_new = matfreq2tau(d, G_q_new, sumrule = sumrule, verbose = verbose)


        if verbose
            if iter % (maxiter / 10) == 0
                println("round $iter: change $(diff(G_x_new, G_x))")
            end
        end
        if maximum(abs.(G_x_new .- G_x)) < tol && iter > 10
            break
        end

        G_x = mix * G_x_new + (1 - mix) * G_x # Linear mixing
    end
    return G_x
end

function printG(d, G_x)
    @printf("%15s%40s%40s%40s\n", "τ", "DLR imag", "DLR real", "asymtotically exact")
    for i in 1:d.size
        if d.τ[i] <= d.β / 2
            @printf("%15.8f%40.15f%40.15f%40.15f\n", d.τ[i], imag(G_x[i]), real(G_x[i]), conformal_tau(d.τ[i], d.β))
        end
    end
    println()
end

verbose = false

printstyled("=====    Prepare the expected Green's function of the SYK model     =======\n", color = :yellow)
dsym_correct = DLRGrid(Euv = 5.0, β = 10000.0, isFermi = true, rtol = 1e-14, symmetry = :ph) # Initialize DLR object
G_x_correct = solve_syk_with_fixpoint_iter(dsym_correct, 0.00, mix = 0.1, verbose = false)
printG(dsym_correct, G_x_correct)

printstyled("=====    Test Symmetrized and Unsymmetrized DLR solver for SYK model     =======\n", color = :yellow)

@printf("%30s%30s%30s%30s%20s\n", "Euv", "sym_solver", "unsym_solver", "unsym_solver+sum_rule", "good or bad")
for Euv in LinRange(5.0, 10.0, 50)

    rtol = 1e-10
    β = 10000.0
    # printstyled("=====     Symmetrized DLR solver for SYK model     =======\n", color = :yellow)
    mix = 0.01
    dsym = DLRGrid(Euv = Euv, β = β, isFermi = true, rtol = rtol, symmetry = :ph, rebuild = true, verbose = false) # Initialize DLR object
    G_x_ph = solve_syk_with_fixpoint_iter(dsym, 0.00, mix = mix, sumrule = nothing, verbose = verbose)

    # printstyled("=====     Unsymmetrized DLR solver for SYK model     =======\n", color = :yellow)
    mix = 0.01
    dnone = DLRGrid(Euv = Euv, β = β, isFermi = true, rtol = rtol, symmetry = :none, rebuild = true, verbose = false) # Initialize DLR object
    G_x_none = solve_syk_with_fixpoint_iter(dnone, 0.00, mix = mix, sumrule = nothing, verbose = verbose)

    # printstyled("=====     Unsymmetrized DLR solver for SYK model     =======\n", color = :yellow)
    mix = 0.01
    G_x_none_sumrule = solve_syk_with_fixpoint_iter(dnone, 0.00, mix = mix, sumrule = 1.0, verbose = verbose)
    # printG(dnone, G_x_none)

    # printstyled("=====     Unsymmetrized versus Symmetrized DLR solver    =======\n", color = :yellow)
    # @printf("%15s%40s%40s%40s\n", "τ", "sym DLR (interpolated)", "unsym DLR", "difference")
    # G_x_interp = tau2tau(dsym_correct, G_x_correct, dnone.τ)
    # for i in 1:dnone.size
    #     if dnone.τ[i] <= dnone.β / 2
    #         @printf("%15.8f%40.15f%40.15f%40.15f\n", dnone.τ[i], real(G_x_interp[i]), real(G_x_none[i]), abs(real(G_x_interp[i] - G_x_none[i])))
    #     end
    # end

    G_x_interp_ph = tau2tau(dsym_correct, G_x_correct, dsym.τ)
    G_x_interp_none = tau2tau(dsym_correct, G_x_correct, dnone.τ)
    G_x_interp_none_sumrule = tau2tau(dsym_correct, G_x_correct, dnone.τ)
    d_ph = diff(G_x_interp_ph, G_x_ph)
    d_none = diff(G_x_interp_none, G_x_none)
    d_none_sumrule = diff(G_x_interp_none_sumrule, G_x_none_sumrule)
    flag = (d_ph < 100rtol) && (d_none < 100rtol) && (d_none_sumrule < 100rtol) ? "good" : "bad"

    @printf("%30.15f%30.15e%30.15e%30.15e%20s\n", Euv, d_ph, d_none, d_none_sumrule, flag)
    # println("symmetric Euv = $Euv maximumal difference: ", diff(G_x_interp, G_x_ph))
    # println("non symmetric Euv = $Euv maximumal difference: ", diff(G_x_interp, G_x_none))

end
