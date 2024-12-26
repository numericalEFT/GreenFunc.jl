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
using LinearAlgebra

const β = 1e4
const J = 1.0
const rtol = 1e-10

diff(a, b) = maximum(abs.(a - b)) # return the maximum deviation between a and b
distance(a, b) = norm(a - b, 2) # return the 1-norm distance between a and b

conformal_tau(τ, β) = π^(1 / 4) / sqrt(2β) * 1 / sqrt(sin(π * τ / β)) #analytic solution with the conformal invariance

const dlrmesh = DLRFreq(β, FERMION; Euv=5.0, rtol=rtol, symmetry=:ph)   # Initialize DLR grid

function selfenergy(Gt)
    ######### calculate sigma ###############
    minus_tau = β .- Gt.mesh[1] # Reversed imaginary time mesh point
    Gt_inv = dlr_to_imtime(to_dlr(Gt), minus_tau) # interpolate into minus_tau grid
    Σt = J .^ 2 .* Gt .^ 2 .* Gt_inv  # SYK self-energy in imaginary time
    return Σt |> to_dlr |> to_imfreq
end

function dyson(Gt)
    ########## sigma --> G  ################
    Σω = selfenergy(Gt)
    freq = matfreq(Σω.mesh[1]) * im
    Gω = 1im * imag.(-1 ./ (freq .+ Σω))

    return Gω |> to_dlr |> to_imtime # Gω --> Gτ

end

function nlsolve(G_t; tol=rtol, maxiter=1000, verbose=false, mix=0.1)
    for iter in 1:maxiter
        G_t_new = dyson(G_t)
        if verbose && (iter % (maxiter / 10) == 0)
            println("round $iter: change $(diff(G_t_new, G_t)), distance $(distance(G_t_new, G_t))")
        end
        if maximum(abs.(G_t_new - G_t)) < tol && iter > 10
            return G_t_new
        end
        G_t = mix .* G_t_new + (1 - mix) .* G_t # Linear mixing
    end
    return G_t
end

const G_t = MeshArray(ImTime(dlrmesh); dtype=ComplexF64)
fill!(G_t, 0.0)
G = nlsolve(G_t, verbose=true)

@printf("%15s%40s%40s%40s\n", "τ", "DLR imag", "DLR real", "asymtotically exact")
for (i, t) in enumerate(G.mesh[1])
    @printf("%15.8f%40.15f%40.15f%40.15f\n", t, imag(G[i]), real(G[i]), conformal_tau(t, β))
end
println()
