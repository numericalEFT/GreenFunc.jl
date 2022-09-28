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

function syk_sigma(mesh_dlr, G_t, J=1.0)
    # tgrid = G_t.mesh[findfirst(x -> (x isa MeshGrids.ImTime), G_t.mesh)]
    tau_rev = reverse(mesh_dlr.β .- mesh_dlr.dlr.τ) # Reversed imaginary time nodes
    G_t_rev = dlr_to_imtime(to_dlr(G_t, mesh_dlr), tau_rev) # G at beta - tau

    Sigma_t = J .^ 2 .* G_t .^ 2 .* G_t_rev # SYK self-energy in imaginary time

    return Sigma_t
end

function dyson(dlr, sigma_freq, mu)
    if dlr.symmetry == :ph #symmetrized G
        @assert mu ≈ 0.0 "Only the case μ=0 enjoys the particle-hole symmetry."
        return 1im * imag.(-1 ./ (dlr.ωn * 1im .- mu .+ sigma_freq))
    elseif dlr.symmetry == :none
        return -1 ./ (dlr.ωn * 1im .- mu .+ sigma_freq)
    else
        error("Not implemented!")
    end
end

function solve_syk_with_fixpoint_iter(mesh_dlr, G_t, mu, tol=mesh_dlr.rtol * 10; mix=0.1, maxiter=5000, verbose=true)
    G_t.data = zeros(ComplexF64, G_t.dims...)
    for iter in 1:maxiter
        sigma_t = syk_sigma(mesh_dlr, G_t)
        sigma_freq = dlr_to_imfreq(to_dlr(sigma_t, mesh_dlr))
        # sigma_freq = sigma_t |> to_dlr |> dlr_to_imfreq

        G_freq_new = dyson(mesh_dlr.dlr, sigma_freq, mu)

        G_t_new = dlr_to_imtime(to_dlr(G_freq_new, mesh_dlr))

        if verbose
            if iter % (maxiter / 10) == 0
                println("round $iter: change $(diff(G_t_new, G_t))")
            end
        end
        if maximum(abs.(G_t_new .- G_t)) < tol && iter > 10
            break
        end

        G_t = mix .* G_t_new + (1 - mix) .* G_t # Linear mixing
    end
    return G_t
end

function printG(dlr, G_t)
    @printf("%15s%40s%40s%40s\n", "τ", "DLR imag", "DLR real", "asymtotically exact")
    for i in 1:dlr.size
        if dlr.τ[i] <= dlr.β / 2
            @printf("%15.8f%40.15f%40.15f%40.15f\n", dlr.τ[i], imag(G_t[i]), real(G_t[i]), conformal_tau(dlr.τ[i], dlr.β))
        end
    end
    println()
end

verbose = false
β = 1e4
isFermi = true

printstyled("=====    Prepare the expected Green's function of the SYK model     =======\n", color=:yellow)

mesh_dlr = MeshGrids.DLRFreq(β, isFermi; Euv=5.0, rtol=1e-14, sym=:ph)   # Initialize DLR grid
mesh = MeshGrids.ImTime(β, isFermi; Euv=5.0, grid=mesh_dlr.dlr.τ)
G_t = MeshArray(mesh; dtype=ComplexF64)

G_t_correct = solve_syk_with_fixpoint_iter(mesh_dlr, G_t, 0.00, mix=0.1, verbose=false)
printG(mesh_dlr.dlr, G_t_correct)

printstyled("=====    Test Symmetrized and Unsymmetrized DLR solver for SYK model     =======\n", color=:yellow)

@printf("%30s%30s%30s%20s\n", "Euv", "sym_solver", "unsym_solver", "good or bad")
for Euv in LinRange(5.0, 10.0, 50)

    rtol = 1e-10
    # printstyled("=====     Symmetrized DLR solver for SYK model     =======\n", color = :yellow)
    mix = 0.01
    mesh_dlrph = MeshGrids.DLRFreq(β, isFermi; Euv=Euv, rtol=rtol, sym=:ph)   # Initialize DLR grid
    mesh1 = MeshGrids.ImTime(β, isFermi; Euv=Euv, grid=mesh_dlrph.dlr.τ)
    G_t1 = MeshArray(mesh1; dtype=ComplexF64)
    G_t_ph = solve_syk_with_fixpoint_iter(mesh_dlrph, G_t1, 0.00, mix=mix, verbose=verbose)
    printG(mesh_dlrph.dlr, G_t_ph)

    # printstyled("=====     Unsymmetrized DLR solver for SYK model     =======\n", color = :yellow)
    mix = 0.01
    mesh_dlrnone = MeshGrids.DLRFreq(β, isFermi; Euv=Euv, rtol=rtol, sym=:none)   # Initialize DLR grid
    mesh2 = MeshGrids.ImTime(β, isFermi; Euv=Euv, grid=mesh_dlrnone.dlr.τ)
    G_t2 = MeshArray(mesh2; dtype=ComplexF64)
    G_t_none = solve_syk_with_fixpoint_iter(mesh_dlrnone, G_t2, 0.00, mix=mix, verbose=verbose)
    printG(mesh_dlrnone.dlr, G_t_none)

    # printstyled("=====     Unsymmetrized versus Symmetrized DLR solver    =======\n", color = :yellow)
    # @printf("%15s%40s%40s%40s\n", "τ", "sym DLR (interpolated)", "unsym DLR", "difference")
    # G_x_interp = tau2tau(dsym_correct, G_x_correct, dnone.τ)
    # for i in 1:dnone.size
    #     if dnone.τ[i] <= dnone.β / 2
    #         @printf("%15.8f%40.15f%40.15f%40.15f\n", dnone.τ[i], real(G_x_interp[i]), real(G_x_none[i]), abs(real(G_x_interp[i] - G_x_none[i])))
    #     end
    # end

    G_t_interp_ph = dlr_to_imtime(to_dlr(G_t_correct, mesh_dlrph))
    G_t_interp_none = dlr_to_imtime(to_dlr(G_t_correct, mesh_dlrnone))
    d_ph = diff(G_t_interp_ph, G_t_ph)
    d_none = diff(G_t_interp_none, G_t_none)
    flag = (d_ph < 100rtol) && (d_none < 100rtol) ? "good" : "bad"

    @printf("%30.15f%30.15e%30.15e%20s\n", Euv, d_ph, d_none, flag)
    # println("symmetric Euv = $Euv maximumal difference: ", diff(G_x_interp, G_x_ph))
    # println("non symmetric Euv = $Euv maximumal difference: ", diff(G_x_interp, G_x_none))

end
