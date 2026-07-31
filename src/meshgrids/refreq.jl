"""
    ReFreq(grid; dtype=Float64)
    ReFreq(ωmin, ωmax, n_points; dtype=Float64)
    ReFreq(; window, n_points, dtype=Float64)

Real-frequency mesh. Explicit vectors may be ascending or descending; the
iteration order follows the input while the stored grid is always ascending.
The `n_w` keyword is accepted as an alias for `n_points`.
"""
struct ReFreq{T<:Real,G<:AbstractGrid{T},REV} <: TemporalGrid{T,REV}
    grid::G
end

function ReFreq(grid::Union{AbstractGrid,AbstractVector}; dtype=Float64)
    converted, rev = _to_realdomain_grid(grid, dtype, "real-frequency")
    return ReFreq{dtype,typeof(converted),rev}(converted)
end

function ReFreq(ωmin::Real, ωmax::Real, n_points::Integer; dtype=Float64)
    n_points >= 2 || throw(ArgumentError("n_points must be at least 2"))
    ωmin < ωmax || throw(ArgumentError("ωmin must be smaller than ωmax"))
    all(isfinite, (ωmin, ωmax)) || throw(ArgumentError("frequency window must be finite"))
    grid = SimpleG.Uniform{dtype}(dtype[ωmin, ωmax], Int(n_points))
    return ReFreq{dtype,typeof(grid),false}(grid)
end

ReFreq(window::Tuple{<:Real,<:Real}, n_points::Integer; kwargs...) =
    ReFreq(window[1], window[2], n_points; kwargs...)

function ReFreq(; window::Tuple{<:Real,<:Real},
    n_points::Union{Nothing,Integer}=nothing,
    n_w::Union{Nothing,Integer}=nothing,
    dtype=Float64)
    isnothing(n_points) == isnothing(n_w) ||
        return ReFreq(window, something(n_points, n_w); dtype=dtype)
    throw(ArgumentError("provide exactly one of n_points or n_w"))
end

Base.show(io::IO, mesh::ReFreq) = print(
    io,
    "Real frequency grid with $(length(mesh)) points: $(_grid(mesh))",
)
