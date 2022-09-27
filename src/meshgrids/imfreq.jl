"""
    struct ImFreq{T<:Real, Grid} <: TemporalGrid{Int}

Imaginary-frequency grid for Green's functions. 

# Parameters
- `T`: type of `β` and `Euv`.
- `Grid`: type of 1D grid with Int as the grid point type.

# Members
- `grid`: 1D grid of time axis, with locate, volume, and AbstractArray interface implemented.
  It should be grid of Int for ImFreq, and DLRGrid for DLRFreq.
- `β`: inverse temperature.
- `Euv`:  the UV energy scale of the spectral density.
- `isFermi`: the statistics for particles is fermionic or not.
"""
struct ImFreq{T<:Real,Grid} <: TemporalGrid{Int}
    grid::Grid
    β::T
    Euv::T
    isFermi::Bool
end

"""
    function ImFreq(β, isFermi::Bool=false;
        dtype=Float64,
        Euv=1000 / β,
        grid::Union{AbstractGrid,AbstractVector,Nothing}=nothing
    )

Create a `ImFreq` struct.

# Arguments
- `β`: inverse temperature.
- `isFermi`: the statistics for particles is fermionic or not. False by default.
- `dtype`: type of `β` and `Euv`. By default, `dtype = Float64`.
- `Euv`: the UV energy scale of the spectral density. By default, `Euv = 1000 / β`.
- `grid`: 1D time grid as a AbstractVector or CompositeGrids.AbstractGrid. By default, a optimized grid built in DLR is used.
"""
function ImFreq(β, isFermi::Bool=false;
    dtype=Float64,
    Euv=1000 / β,
    rtol=1e-12,
    grid::Union{AbstractGrid,AbstractVector,Nothing}=nothing
)
    if isnothing(grid)
        dlr = DLRGrid(Euv, β, rtol, isFermi, :none)
        grid = SimpleG.Arbitrary{Int}(dlr.n)
    elseif (grid isa AbstractVector)
        grid = SimpleG.Arbitrary{Int}(Int.(grid))
    end
    return ImFreq{dtype,typeof(grid)}(grid, β, Euv, isFermi)
end

matfreq_to_int(tg::ImFreq, ωn) = tg.isFermi ? Int(round((ωn * tg.β / π - 1) / 2)) : Int(round((ωn * tg.β / π) / 2))
int_to_matfreq(tg::ImFreq, n::Int) = tg.isFermi ? (2n + 1) * π / tg.β : 2n * π / tg.β

"""
    getindex(g::ImFreq, I::Int)

Equivalent to `g[I]`, get the __real-valued__ Matsubara frequency of the Ith point in the grid. 
For fermion, return (2g[I]+1)π/β, for boson, return 2g[I]*π/β.

If you need the __integer-valued__ frequency, use `g.grid[I]` instead.
"""
Base.getindex(tg::ImFreq, I::Int) = int_to_matfreq(tg, tg.grid[I])

"""
    show(io::IO, tg::ImFreq)

Write a text representation of the Imaginary-frequency grid `tg` to the output stream `io`.
"""
Base.show(io::IO, tg::ImFreq) = print(io, "Matsubara frequency grid with $(length(tg)) points, inverse temperature = $(tg.β), UV Energy scale = $(tg.Euv), fermionic = $(tg.isFermi): $(_grid(tg.grid))")


"""
    locate(tg::ImFreq, n::Int)
    locate(tg::ImFreq, ωn)

Find the location in `tg.grid` for the Matsubara frequency `ωn` or the integer `n`.
"""
locate(tg::ImFreq, n::Int) = locate(tg.grid, n)
function locate(tg::ImFreq, ωn)
    n = matfreq_to_int(tg, ωn)
    return locate(tg.grid, n)
end
