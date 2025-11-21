"""
    struct MeshArray{T,N,MT} <: AbstractMeshArray{T,N}

Multi-dimensional array that is defined on a mesh.
The mesh is a tuple of meshgrid objects.
The mesh is stored in the field `mesh` and the data is stored in the field `data`.

# Parameters:
- `T`: type of data
- `MT`: type of mesh, e.g., `Tuple{MeshType1, MeshType2, ...}`
- `N`: number of dimensions

# Members:
- `mesh` (`MT`): the mesh is a tuple of meshes.
   The mesh should be an iterable object that contains an ordered list of grid points.
   Examples are the
   1. Meshes defined in the `MeshGrids` module.
   2. UnitRange such as `1:10`, etc.
   3. Product of meshes `MeshProduct` defined in the `MeshGrids` module.

   If a mesh is defined on a continuous manifold and supports the following methods, then one can perform interpolation, derivatives, etc. on the mesh:
    - `locate(mesh, value)`: find the index of the closest grid point for given value;
    - `volume(mesh, index)`: find the volume of grid space near the point at griven index.
    - `volume(mesh, gridpoint)`: locate the corresponding index of a given grid point and than find the volume spanned by the grid point.

- `data` (`Array{T,N}`): the data.
- `dims`: dimension of the data
"""
struct MeshArray{T,N,MT} <: AbstractMeshArray{T,N}
    mesh::MT
    data::Array{T,N}
    dims::NTuple{N,Int}
    function MeshArray{T,N,MT}(data::AbstractArray{T,N}, mesh::MT) where {T,N,MT}
        # do nothing constructor, so that it is fast with no additional allocation
        # but you need to make sure that the mesh and data make sense

        # @assert mesh isa Tuple "The mesh should be wrappered with a tuple."
        # for (i, s) in enumerate(size(data))
        #     @assert length(mesh[i]) == s "The size of data and the $(i)th mesh do not match."
        # end
        return new{T,N,MT}(mesh, data, size(data))
    end
end


"""
    MeshMatrix{T}
Alias for [`MeshArray{T,2,MT}`](@ref MeshArray).
"""
const MeshMatrix{T,MT} = MeshArray{T,2,MT}

"""
    MeshVector{T}
Alias for [`MeshArray{T,1,MT}`](@ref MeshArray).
"""
const MeshVector{T,MT} = MeshArray{T,1,MT}

# =====================================================
# Type stabilization helper functions (Issue #78 fix)
# =====================================================

# Type stabilization helper function (function barrier)
function _stabilize_mesh_type(mesh::Tuple)
    # Concretize the mesh tuple type
    if isconcretetype(typeof(mesh))
        return mesh
    # else
        # NOTE: This branch is unreachable in practice because typeof() always returns
        # the actual runtime type, which is concrete. Kept for defensive programming.
        # Create a new tuple while preserving the type of each element
        # return map(identity, mesh)
    end
end

function _stabilize_mesh_type(mesh)
    # Convert non-tuple to tuple
    return tuple(mesh...)
end

# Type-stable internal constructor (function barrier)
@inline function _create_mesharray_typed(data::AbstractArray{T,N}, mesh::MT, ::Type{T}, ::Val{N}) where {T,N,MT}
    return MeshArray{T,N,MT}(data, mesh)
end

function _create_mesharray_typed(data::AbstractArray{T,N}, mesh, dtype::Type, n::Int) where {T,N}
    # Handle type conversion when necessary
    if T != dtype
        data_converted = convert(Array{dtype,N}, data)
        return _create_mesharray_typed(data_converted, mesh, dtype, Val(n))
    else
        return _create_mesharray_typed(data, mesh, T, Val(N))
    end
end

# Generated function for improved type inference
@generated function _infer_mesh_type(mesh::MT) where {MT}
    return :(MT)
end

# =====================================================
# Main constructor (type-stabilized)
# =====================================================

"""
    function MeshArray(;
        mesh...;
        dtype = Float64,
        data::Union{Nothing,AbstractArray}=nothing) where {T}

Create a Green struct. Its memeber `dims` is setted as the tuple consisting of the length of all meshes.

# Arguments
- `mesh`: meshes of Green's function. See the docs of [`MeshArray`](@ref) for more details.
   Mesh could be any iterable object, examples are vector, tuple, array, number, UnitRange (say, 1:5).
- `dtype`: data type of Green's function's value.
- `data`: the data of the Green's function. By default, `data` is constructed to an unintialized Array with the `dims` size containing elements of `dtype`.
"""
function MeshArray(mesh...;
    dtype=Float64,
    data::Union{Nothing,AbstractArray}=nothing)

    @assert all(x -> isiterable(typeof(x)), mesh) "all meshes should be iterable"

    mesh = _stabilize_mesh_type(mesh)

    # NOTE: This check is unreachable in practice because typeof() always returns
    # the actual runtime type, which is concrete. Kept for defensive programming.
    # if isconcretetype(typeof(mesh)) == false
    #     @warn "Mesh type $(typeof(mesh)) is not concrete, it may cause performance issue."
    # end

    N = length(mesh)
    dims = tuple([length(v) for v in mesh]...)
    if isnothing(data) == false
        # @assert length(size(data)) == N
        @assert size(data) == dims "data size $(size(data)) should be the same as the mesh size $dims"
    else
        data = Array{dtype,N}(undef, dims...)
    end

    if dtype != eltype(data)
        data = convert(Array{dtype,N}, data)
    end
    return _create_mesharray_typed(data, mesh, dtype, N)
end
function MeshArray(; mesh::Union{Tuple,AbstractVector},
    dtype=Float64,
    data::Union{Nothing,AbstractArray}=nothing)

    @assert all(x -> isiterable(typeof(x)), mesh) "all meshes should be iterable"

    if mesh isa AbstractVector
        mesh = tuple(mesh...)
    end

    mesh = _stabilize_mesh_type(mesh)

    # NOTE: This check is unreachable in practice because typeof() always returns
    # the actual runtime type, which is concrete. Kept for defensive programming.
    # if isconcretetype(typeof(mesh)) == false
    #     @warn "Mesh type $(typeof(mesh)) is not concrete, it may cause performance issue."
    # end

    @assert mesh isa Tuple "mesh should be a tuple, now get $(typeof(mesh))"
    N = length(mesh)
    dims = tuple([length(v) for v in mesh]...)
    if isnothing(data) == false
        # @assert length(size(data)) == N
        @assert size(data) == dims "data size $(size(data)) should be the same as the mesh size $dims"
    else
        data = Array{dtype,N}(undef, dims...)
    end

    if dtype != eltype(data)
        data = convert(Array{dtype,N}, data)
    end
    return _create_mesharray_typed(data, mesh, dtype, N)
end

"""
    getindex(obj::MeshArray, inds...)

Return a subset of `obj`'s data as specified by `inds`, where each `inds` may be, for example, an Int, an AbstractRange, or a Vector.
"""
@inline Base.getindex(obj::MeshArray{T,N,MT}, inds::Vararg{Int,N}) where {T,MT,N} = @inbounds Base.getindex(obj.data, inds...)
# Base.getindex(obj::MeshArray, I::Int) = Base.getindex(obj.data, I)

"""
    setindex!(obj::MeshArray, v, inds...)
    obj[inds...] = v

Store values from array `v` within some subset of `obj.data` as specified by `inds`.
"""
@inline Base.setindex!(obj::MeshArray{T,N,MT}, v, inds::Vararg{Int,N}) where {T,MT,N} = @inbounds Base.setindex!(obj.data, v, inds...)
# Base.setindex!(obj::MeshArray, v, I::Int) = Base.setindex!(obj.data, v, I)

# IndexStyle(::Type{<:MeshArray}) = IndexCartesian() # by default, it is IndexCartesian

"""
    Base.similar(obj::MeshArray{T,N,MT}, ::Type{S}) where {T,MT,N,S}
    Base.similar(obj::MeshArray{T,N,MT}) where {T,MT,N} = Base.similar(obj, T)

# Return type:
- `Base.similar(obj::MeshArray)`: Return a new MeshArray with the same meshes, and the uninitialized data of the same type as `obj.data`.
- `Base.similar(obj::MeshArray, ::Type{S})`: Return a new MeshArray with the same meshes, but with the uninitialized data of type `S`.
"""
function Base.similar(obj::MeshArray{T,N,MT}, ::Type{S}) where {T,MT,N,S}
    return MeshArray(mesh=obj.mesh, dtype=S, data=similar(obj.data, S))
end
Base.similar(obj::MeshArray{T,N,MT}) where {T,MT,N} = Base.similar(obj, T)
#By default, the following functions will all call Base.similar(obj::MeshArray, ::Type{S}, inds)
#as explained in https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-array
#`Base.similar(obj::MeshArray, ::Type{S}, inds)`: Return a slice of obj.data.
#However, we don't want that since slice of GreeNew itself is not well defined with meshes.

################################ broadcast interface ###############################################
Base.BroadcastStyle(::Type{<:MeshArray}) = Broadcast.ArrayStyle{MeshArray}()

function Base.similar(bc::Base.Broadcast.Broadcasted{Broadcast.ArrayStyle{MeshArray}}, ::Type{ElType}) where {ElType}
    # Scan the inputs for the MeshArray:
    A = find_gf(bc)
    # Use other fields of A to create the output
    data = similar(Array{ElType}, axes(bc))
    m = MeshArray{ElType,length(axes(bc)),typeof(A.mesh)}(data, A.mesh)
    return m
end

find_gf(bc::Broadcast.Broadcasted) = find_gf(bc.args)
find_gf(args::Tuple) = find_gf(find_gf(args[1]), Base.tail(args))
find_gf(x) = x
find_gf(::Tuple{}) = nothing
find_gf(a::MeshArray, rest) = a
find_gf(::Any, rest) = find_gf(rest)

# Type-stable broadcast implementation
function Base.copyto!(dest::MeshArray{T,N,MT}, bc::Base.Broadcast.Broadcasted) where {T,N,MT}
    # without this function, inplace operation like g1 .+= g2 will make a lot of allocations
    # Please refer to the following posts for more details:
    # 1. manual on the interface: https://docs.julialang.org/en/v1/manual/interfaces/#extending-in-place-broadcast-2
    # 2. see the post: https://discourse.julialang.org/t/help-implementing-copyto-for-broadcasting/51204/3
    # 3. example from DataFrames.jl: https://github.com/JuliaData/DataFrames.jl/blob/main/src/other/broadcasting.jl#L193

    # Type stabilization: @simd and inlining
    indices = CartesianIndices(dest.data)
    bcf = Base.Broadcast.flatten(bc)

    # Call type-stable internal function (function barrier)
    _copyto_typed!(dest, bcf, indices)

    return dest
end

@inline function _copyto_typed!(dest::MeshArray{T,N,MT}, bcf, indices::CartesianIndices{N}) where {T,N,MT}
    # Fast copy with determined types
    @inbounds @simd for I in indices
        dest.data[I] = bcf[I]
    end
    return nothing
end

########### alternative approach ######################
# function Base.copyto!(dest::MeshArray{T, N, MT}, bc::Base.Broadcast.Broadcasted{Nothing}) where {T,MT,N}
#     _bcf = Base.Broadcast.flatten(bc)
#     bcf = Base.Broadcast.preprocess(dest, _bcf)
#     for I in CartesianIndices(dest)
#         dest[I] = bcf[I]
#     end
#     return dest
# end

function Base.show(io::IO, obj::MeshArray)
    print(io, "Meshed array with dims = $(obj.dims) and total length = $(length(obj.data))\n",
              "- Mesh: $(typeof(obj.mesh)) \n")
end

Base.show(io::IO, ::MIME"text/plain", obj::MeshArray) = Base.show(io, obj)
Base.show(io::IO, ::MIME"text/html", obj::MeshArray) = Base.show(io, obj)

"""
    function rank(obj::MeshArray)

Return the dimension of `obj.data` (`N`).
"""
rank(::Type{MeshArray{T,N,MT}}) where {T,MT,N} = N


#TODO:Following triqs design, we want the following two things:
# G<<init_function will initiate G with the given function.
# G<<g2 will copy g2 into G.


# """
#    def __lshift__(self, A):
#         A can be two things:
#           * G << any_init will init the GFBloc with the initializer
#           * G << g2 where g2 is a GFBloc will copy g2 into self
#         if isinstance(A, Gf):
#             if self is not A: # otherwise it is useless AND does not work !!
#                 assert self.mesh == A.mesh, "Green function meshes are not compatible:\n  %s\nand\n  %s" % (self.mesh, A.mesh)
#                 self.copy_from(A)
#         elif isinstance(A, lazy_expressions.LazyExpr): # A is a lazy_expression made of GF, scalars, descriptors
#             A2 = descriptors.convert_scalar_to_const(A)
#             def e_t (x):
#                 if not isinstance(x, descriptors.Base): return x
#                 tmp = self.copy()
#                 x(tmp)
#                 return tmp
#             self.copy_from (lazy_expressions.eval_expr_with_context(e_t, A2) )
#         elif isinstance(A, lazy_expressions.LazyExprTerminal): #e.g. g<< SemiCircular (...)
#             self << lazy_expressions.LazyExpr(A)
#         elif descriptors.is_scalar(A): #in the case it is a scalar ....
#             self << lazy_expressions.LazyExpr(A)
#         else:
#             raise NotImplemented
#                     return self
# """

"""
    function _check(objL::MeshArray, objR::MeshArray)

Check if the Green's functions `objL` and `objR` are on the same meshes. Throw an AssertionError if any check is false.
"""
function _check(objL::MeshArray, objR::MeshArray)
    # KUN: check --> __check
    # first:  check typeof(objL.tgrid)==typeof(objR.tgrid)
    # second: check length(objL.tgrid)
    # third:  hasmethod(objL.tgrid, isequal) --> assert
    # @assert objL.innerstate == objR.innerstate "Green's function innerstates are not inconsistent: $(objL.innerstate) and $(objR.innerstate)"
    @assert typeof(objL.mesh) == typeof(objR.mesh) "Green's function meshes' types are inconsistent: $(typeof(objL.mesh)) and $(typeof(objR.mesh))"
    @assert objL.dims == objR.dims "Green's function dims are inconsistent: $(objL.dims) and $(objR.dims)"

    return true
    # @assert objL.tgrid == objR.tgrid "Green's function time grids are not compatible:\n $(objL.tgrid)\nand\n $(objR.tgrid)"
    # @assert objL.mesh == objR.mesh "Green's function meshes are not compatible:\n $(objL.mesh)\nand\n $(objR.mesh)"
end