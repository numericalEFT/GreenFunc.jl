"""
    struct ImTime{T, G, R} <: TemporalGrid{T}

Time grid for Green's functions.

# Parameters
- `T<:Real`: type of the `grid` point, `β` and `Euv`.
- `G<:AbstractGrid{T}`: type of 1D grid with `T` as the grid point type.
- `R`: type of the representation.
- REV: access the grid in reverse order or not.

# Members
- `grid`: 1D grid of time axis, with locate, volume, and AbstractArray interface implemented. Always in ascend order.
  It should be grid of Int for ImFreq, and DLRGrid for DLRFreq.
- `β`: inverse temperature.
- `Euv`:  the UV energy scale of the spectral density.
- `isFermi`: the statistics for particles is fermionic or not.
- `symmetry`: `:ph` for particle-hole symmetric, `:pha` for particle-hole symmetry, and `:none` for no symmetry. By default, `sym = :none`.
- `rtol`: relative tolerance
- `representation`: the representation of the Green's function.
"""
struct ImTime{T<:Real,G<:AbstractGrid{T},R,REV} <: TemporalGrid{T,REV}
    grid::G
    β::T
    Euv::T
    isFermi::Bool
    symmetry::Symbol
    rtol::T
    representation::R
end

"""
    function ImTime(β, isFermi::Bool=false;
        dtype=Float64,
        rtol=1e-12,
        Euv=1000 / β,
        symmetry=:none,
        grid::Union{AbstractGrid,AbstractVector,Nothing}=nothing
    )

Create a `ImTime` struct.

# Arguments
- `β`: inverse temperature.
- `isFermi`: the statistics for particles is fermionic or not. False by default.
- `dtype`: type of the `grid` point. By default, `dtype = Float64`.
- `Euv`: the UV energy scale of the spectral density. By default, `Euv = 1000 / β`.
- `symmetry`: `:ph` for particle-hole symmetric, `:pha` for particle-hole symmetry, and `:none` for no symmetry. By default, `sym = :none`.
- `grid`: 1D time grid as a AbstractVector or CompositeGrids.AbstractGrid. By default, a optimized grid built in DLR is used.
"""
function ImTime(β, isFermi::Bool=false;
    dtype=Float64,
    rtol=1e-12,
    Euv=1000 / β,
    symmetry=:none,
    grid::Union{AbstractGrid,AbstractVector,Nothing}=nothing
)
    dlr = DLRGrid(Euv, β, rtol, isFermi, symmetry)

    if isnothing(grid)
        rev = false
        grid = SimpleG.Arbitrary{dtype}(dlr.τ)
        # grid = SimpleG.Uniform{dtype}([0, β], Int(round(β / resolution)))
        # grid = CompositeGrid.LogDensedGrid(:uniform, [0.0, β], [0.0, β], 8, 1 / Euv, 8) #roughly ~100 points if resolution = β/128
    elseif (grid isa AbstractVector)
        # if !(grid isa AbstractGrid)
        #     rev = issorted(grid, rev=true)
        #     if rev
        #         grid = reverse(grid)
        #     end
        #     grid = SimpleG.Arbitrary{dtype}(grid)
        # end
        grid, rev = _to_AbstractGrid(grid, dtype)
    else
        error("Proper grid and basis are required!")
    end
    @assert eltype(grid) == dtype "The type of grid should be the same as dtype = $dtype"
    @assert issorted(grid) "The grid should be sorted."
    @assert grid[1] >= 0 && grid[end] <= β "The grid should be in the range [0, β]."
    return ImTime{dtype,typeof(grid),typeof(dlr),rev}(grid, β, Euv, isFermi, symmetry, rtol, dlr)
end

"""
    function ImTime(dlr::DLRGrid;
                    dtype=Float64,
                    grid::Union{AbstractGrid,AbstractVector}=SimpleG.Arbitrary{dtype}(dlr.τ))

Construct `ImTime` from a `DLRGrid`, with a given `grid`. By default, `grid` is the imaginary-time grid points from `DLRGrid`.
"""
function ImTime(dlr::DLRGrid;
    dtype=Float64,
    grid::Union{AbstractGrid,AbstractVector}=SimpleG.Arbitrary{dtype}(dlr.τ)
)
    rev = issorted(grid, rev=true)
    if rev
        grid = reverse(grid)
    end
    # if isnothing(grid)
    #     grid = SimpleG.Arbitrary{dtype}(dlr.τ)
    # end
    if (grid isa AbstractGrid) == false
        grid = SimpleG.Arbitrary{dtype}(grid)
    end
    @assert eltype(grid) == dtype "The type of grid should be the same as dtype = $dtype"
    @assert issorted(grid) "The grid should be sorted."
    return ImTime{dtype,typeof(grid),typeof(dlr),rev}(grid, dlr.β, dlr.Euv, dlr.isFermi, dlr.symmetry, dlr.rtol, dlr)
end
ImTime(dlrfreq::DLRFreq; kwargs...) = ImTime(dlrfreq.dlr; kwargs...)

"""
    show(io::IO, tg::ImTime)

Write a text representation of the Imaginary-time grid `tg` to the output stream `io`.
"""
Base.show(io::IO, tg::ImTime) = print(io, "Imaginary Time grid with $(length(tg)) points, inverse temperature = $(tg.β), UV Energy scale = $(tg.Euv), fermionic = $(tg.isFermi): $(_grid(tg.grid))")
