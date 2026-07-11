"""
    ReTime(grid; dtype=Float64)
    ReTime(tmin, tmax, n_points; dtype=Float64)
    ReTime(; window, n_points, dtype=Float64)

Real-time mesh. Explicit vectors may be ascending or descending; the iteration
order follows the input while the stored grid is always ascending. The `n_t`
keyword is accepted as an alias for `n_points`.
"""
struct ReTime{T<:Real,G<:AbstractGrid{T},REV} <: TemporalGrid{T,REV}
    grid::G
end

function ReTime(grid::Union{AbstractGrid,AbstractVector}; dtype=Float64)
    converted, rev = _to_realdomain_grid(grid, dtype, "real-time")
    return ReTime{dtype,typeof(converted),rev}(converted)
end

function ReTime(tmin::Real, tmax::Real, n_points::Integer; dtype=Float64)
    n_points >= 2 || throw(ArgumentError("n_points must be at least 2"))
    tmin < tmax || throw(ArgumentError("tmin must be smaller than tmax"))
    all(isfinite, (tmin, tmax)) || throw(ArgumentError("time window must be finite"))
    grid = SimpleG.Uniform{dtype}(dtype[tmin, tmax], Int(n_points))
    return ReTime{dtype,typeof(grid),false}(grid)
end

ReTime(window::Tuple{<:Real,<:Real}, n_points::Integer; kwargs...) =
    ReTime(window[1], window[2], n_points; kwargs...)

function ReTime(; window::Tuple{<:Real,<:Real},
    n_points::Union{Nothing,Integer}=nothing,
    n_t::Union{Nothing,Integer}=nothing,
    dtype=Float64)
    isnothing(n_points) == isnothing(n_t) ||
        return ReTime(window, something(n_points, n_t); dtype=dtype)
    throw(ArgumentError("provide exactly one of n_points or n_t"))
end

Base.show(io::IO, mesh::ReTime) = print(
    io,
    "Real time grid with $(length(mesh)) points: $(_grid(mesh))",
)
