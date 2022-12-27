# basic AbstractArray implement
Base.length(tg::TemporalGrid) = length(tg.grid)
Base.size(tg::TemporalGrid) = size(tg.grid)
Base.size(tg::TemporalGrid, I::Int) = size(tg.grid, I)

Base.firstindex(tg::TemporalGrid) = 1
Base.lastindex(tg::TemporalGrid) = length(tg)

Base.getindex(tg::TemporalGrid{T,false}, I::Int) where {T} = tg.grid[I] # ascend order
Base.getindex(tg::TemporalGrid{T,true}, I::Int) where {T} = tg.grid[end-I+1] # descend order

# iterator
Base.iterate(tg::TemporalGrid) = (tg[1], 1)
Base.iterate(tg::TemporalGrid, state) = (state >= length(tg)) ? nothing : (tg[state+1], state + 1)
# Base.IteratorSize(tg)
Base.IteratorSize(::Type{TemporalGrid{GT,REV}}) where {GT,REV} = Base.HasLength()
Base.IteratorEltype(::Type{TemporalGrid{GT,REV}}) where {GT,REV} = Base.HasEltype()
Base.eltype(::Type{TemporalGrid{GT,REV}}) where {GT,REV} = eltype(GT)



# locate and volume could fail if tg.grid has no implementation
volume(tg::TemporalGrid) = volume(tg.grid)

#ascend order
volume(tg::TemporalGrid{T,false}, I::Int) where {T} = volume(tg.grid, I)
locate(tg::TemporalGrid{T,false}, pos) where {T} = locate(tg.grid, pos)

#descend order
volume(tg::TemporalGrid{T,true}, I::Int) where {T} = volume(tg.grid, length(tg) - I + 1)
locate(tg::TemporalGrid{T,true}, pos) where {T} = length(tg) - locate(tg.grid, pos) + 1 #TODO: how to implement?

"""
    Base.floor(tg::TemporalGrid{T,false}, pos) where {T} = floor(tg.grid, pos) #TODO: how to implement?
    Base.floor(tg::TemporalGrid{T,true}, pos) where {T} = length(tg) - floor(tg.grid, pos) #TODO: how to implement?

If the grid is in ascend order, then floor returns the largest index that the grid point is smaller than pos.
If the grid is in descend order, then floor returns the largest index that the grid point is larger than pos.

In both cases, the returned index is in the range [1, length(tg)-1]
"""
Base.floor(tg::TemporalGrid{T,false}, pos) where {T} = floor(tg.grid, pos) #TODO: how to implement?
Base.floor(tg::TemporalGrid{T,true}, pos) where {T} = length(tg) - floor(tg.grid, pos) #TODO: how to implement?

is_reverse(tg::TemporalGrid{T,REV}) where {T,REV} = REV

function _round(grid, sigdigits)
    if sigdigits <= 0
        return grid
    else
        return [round(x, sigdigits=sigdigits) for x in grid]
    end
end

# return a pretty print for grid. 
# If grid is very long, return [grid[1], grid[2], grid[3], ..., grid[end-2], grid[end-1], grid[end]]
# If grid is short, return the entire grid
function _grid(grid, n=3)
    if eltype(grid) <: Int
        # return join(io, ["[", join([grid[1:n]..., "...", grid[end-n:end]...], ", "), "]"])
        digits = 0
    else
        resolution = grid[2] - grid[1]
        digits = Int(round(log(grid[end] / resolution) / log(10))) + 3
        digits = digits < 5 ? 5 : digits
    end
    if length(grid) <= 2n + 3
        return join(["[", join(_round(grid, digits), ", "), "]"])
    else
        return join(["[", join([_round(grid[1:n], digits)..., "...", _round(grid[end-(n-1):end], digits)...], ", "), "]"])
    end
end

Base.show(io::IO, ::MIME"text/plain", tg::TemporalGrid) = Base.show(io, tg)
Base.show(io::IO, ::MIME"text/html", tg::TemporalGrid) = Base.show(io, tg)