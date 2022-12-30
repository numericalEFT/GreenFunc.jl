module MeshArrays
abstract type AbstractMeshArray{T,N} <: AbstractArray{T,N} end

import ..MeshGrids

isiterable(::Type{T}) where {T} = hasmethod(iterate, (T,))

include("dense.jl")

export MeshArray, MeshMatrix, MeshVector

########## Array Interface: https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-array #############

"""
    size(obj::AbstractMeshArray)

Return a tuple containing the dimensions of `obj.data` (`obj.dims`).
"""
Base.size(obj::AbstractMeshArray) = obj.dims

"""
    eltype(obj::AbstractMeshArray)

Return the type of the elements contained in `obj.data`.
"""
Base.eltype(::Type{AbstractMeshArray{T,N}}) where {T,N} = T

function Base.zero(obj::AbstractMeshArray{T,N}) where {T,N}
    _obj = similar(obj)
    _obj.data .= zero(T)
    return _obj
end

end