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

verbose = false
β = 1e4
isFermi = true
J = 1.0
rtol = 1e-10

diff(a, b) = maximum(abs.(a - b)) # return the maximum deviation between a and b
distance(a, b) = norm(a - b, 2) # return the 1-norm distance between a and b

conformal_tau(τ, β) = π^(1 / 4) / sqrt(2β) * 1 / sqrt(sin(π * τ / β))
reverseview(x) = view(x, reverse(axes(x, 1))) # reversed view of x

const dlrmesh = DLRFreq(β, isFermi; Euv=5.0, rtol=rtol, sym=:ph, rebuild=false)   # Initialize DLR grid

function selfenergy(Gt)
    ######### calculate sigma ###############
    minus_tau = reverse(β .- Gt.mesh[1]) # Reversed imaginary time mesh point
    Gt_inv = dlr_to_imtime(to_dlr(Gt, dlrmesh), minus_tau) # G at beta - tau
    Gmt = reverseview(Gt_inv)
    Σt = J .^ 2 .* Gt .^ 2 .* Gmt  # SYK self-energy in imaginary time
    return dlr_to_imfreq(to_dlr(Σt, dlrmesh))
end

function dyson(Gt)
    ########## sigma --> G  ################
    Σω = selfenergy(Gt)
    freq = matfreq(Σω.mesh[1]) * im
    Gω = 1im * imag.(-1 ./ (freq .+ Σω))

    #############  Gω --> Gτ  ###############
    # return dlr_to_imtime(to_dlr(Gω, dlrmesh))
    println(size(Gω))
    Gdlr = Gω |> to_dlr
    println(Gdlr.mesh[1].dlr)
    println(size(Gdlr))
    Gnew = Gdlr |> dlr_to_imtime
    println(size(Gnew))
    return Gnew
    # return Gω |> to_dlr |> dlr_to_imtime
    # return dlr_to_imtime(to_dlr(Gω, dlrmesh))

end

function nlsolve(G_t, tol=rtol, maxiter=1000, verbose=false, mix=0.1)
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
end

const mesh = ImTime(β, isFermi; Euv=5.0, grid=dlrmesh.dlr.τ)
const G_t = MeshArray(mesh; dtype=ComplexF64)
fill!(G_t, 0.0)
G = nlsolve(G_t)

@printf("%15s%40s%40s%40s\n", "τ", "DLR imag", "DLR real", "asymtotically exact")
for i in eachindex(G)
    # if dlr.τ[i] <= dlr.β / 2
    @printf("%15.8f%40.15f%40.15f%40.15f\n", mesh[i], imag(G[i]), real(G[i]), conformal_tau(mesh[i], mesh.β))
    # end
end
println()

