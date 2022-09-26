module MeshArrays
abstract type AbstractMeshArray{T,N} <: AbstractArray{T,N} end

import ..MeshGrids

isiterable(::Type{T}) where {T} = hasmethod(iterate, (T,))

include("dense.jl")

export MeshArray

end