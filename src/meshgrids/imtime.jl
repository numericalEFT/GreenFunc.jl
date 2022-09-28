"""
    struct ImTime{T, Grid} <: TemporalGrid{T}

Time grid for Green's functions.

# Parameters
- `T`: type of the `grid` point, `β` and `Euv`.
- `Grid`: type of 1D grid with `T` as the grid point type.

# Members
- `grid`: 1D grid of time axis, with locate, volume, and AbstractArray interface implemented.
  It should be grid of Int for ImFreq, and DLRGrid for DLRFreq.
- `β`: inverse temperature.
- `Euv`:  the UV energy scale of the spectral density.
- `isFermi`: the statistics for particles is fermionic or not.
"""
struct ImTime{T<:Real,Grid} <: TemporalGrid{T}
    grid::Grid
    β::T
    Euv::T
    isFermi::Bool
end

"""
    function ImTime(β, isFermi::Bool=false;
        dtype=Float64,
        rtol=1e-12,
        Euv=1000 / β,
        grid::Union{AbstractGrid,AbstractVector,Nothing}=nothing
    )

Create a `ImTime` struct.

# Arguments
- `β`: inverse temperature.
- `isFermi`: the statistics for particles is fermionic or not. False by default.
- `dtype`: type of the `grid` point. By default, `dtype = Float64`.
- `Euv`: the UV energy scale of the spectral density. By default, `Euv = 1000 / β`.
- `grid`: 1D time grid as a AbstractVector or CompositeGrids.AbstractGrid. By default, a optimized grid built in DLR is used.
"""
function ImTime(β, isFermi::Bool=false;
    dtype=Float64,
    rtol=1e-12,
    Euv=1000 / β,
    grid::Union{AbstractGrid,AbstractVector,Nothing}=nothing
)
    if isnothing(grid)
        dlr = DLRGrid(Euv, β, rtol, isFermi, :none)
        grid = SimpleG.Arbitrary{dtype}(dlr.τ)
        # grid = SimpleG.Uniform{dtype}([0, β], Int(round(β / resolution)))
        # grid = CompositeGrid.LogDensedGrid(:uniform, [0.0, β], [0.0, β], 8, 1 / Euv, 8) #roughly ~100 points if resolution = β/128
    elseif (grid isa AbstractVector)
        grid = SimpleG.Arbitrary{dtype}(grid)
    end
    @assert grid[1] >= 0 && grid[end] <= β "The grid should be in the range [0, β]."
    @assert issorted(grid) "The grid should be sorted."
    @assert eltype(grid) == dtype "The type of grid should be the same as dtype = $dtype"
    return ImTime{dtype,typeof(grid)}(grid, β, Euv, isFermi)
end

"""
    show(io::IO, tg::ImTime)

Write a text representation of the Imaginary-time grid `tg` to the output stream `io`.
"""
Base.show(io::IO, tg::ImTime) = print(io, "Imaginary Time grid with $(length(tg)) points, inverse temperature = $(tg.β), UV Energy scale = $(tg.Euv), fermionic = $(tg.isFermi): $(_grid(tg.grid))")
