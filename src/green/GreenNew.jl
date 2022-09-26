# abstract type AbstractGreen{T,N,NINNER} <: AbstractArray{T,N} end

isiterable(::Type{T}) where {T} = hasmethod(iterate, (T,))
# isiterable(T) = hasmethod(iterate, (typeof(T),))

"""
    mutable struct GreenNew{T,N,MT} <: AbstractArray{T,N}

General Green's function on a multi-dimensional mesh.

# Parameters:
- `T`: data type of the Green's function's value.
- `N`: dimension of `data`.
- `MT`: type of mesh.

# Members:
- `mesh` (MT): a direct product of grids of all degrees of freedom of Green's function. The mesh has to support all standard Base functions of AbstractArray, plus the following two:
    - locate(`mesh`, value): find the index of the closest grid point for given value;
    - volume(`mesh`, index): find the volume of grid space near the point at griven index.
    - volume(`mesh`, gridvalue): locate the corresponding index of a given value and then find the volume of grid space. 
- `data` (Array{T,N}): the data of the Green's function.
- `dims` (NTuple{N,Int}): the tuple consisting of the lengths of all meshes. Its integer arguments must be equal to the lengths in each dimension of `data`.
"""
mutable struct GreenNew{T,N,MT} <: AbstractArray{T,N}
    #########   Mesh   ##############
    mesh::MT
    data::Array{T,N}
    dims::NTuple{N,Int}
end

"""
    function GreenNew(mesh...;
        dtype=Float64,
        data::Union{Nothing,AbstractArray}=nothing)
    
Create a Green struct. Its memeber `dims` is setted as the tuple consisting of the length of all meshes.

# Arguments
- `mesh`: meshes of Green's function. Mesh could be any iterable object, examples are vector, tuple, array, number, UnitRange (say, 1:5).
- `dtype`: data type of Green's function's value.
- `data`: the data of the Green's function. By default, `data` is constructed to an unintialized Array with the `dims` size containing elements of `dtype`.
"""
function GreenNew(mesh...;
    dtype=Float64,
    data::Union{Nothing,AbstractArray}=nothing)

    @assert all(x -> isiterable(typeof(x)), mesh) "all meshes should be iterable"

    N = length(mesh)
    dims = tuple([length(v) for v in mesh]...)

    if isnothing(data) == false
        # @assert length(size(data)) == N
        @assert size(data) == dims
    else
        data = Array{dtype,N}(undef, dims...)
    end

    return GreenNew{dtype,N,typeof(mesh)}(mesh, data, dims)
end

########## Array Interface: https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-array #############
"""
    size(obj::GreenNew)

Return a tuple containing the dimensions of `obj.data` (`obj.dims`).
"""
Base.size(obj::GreenNew) = obj.dims

"""
    eltype(obj::GreenNew)

Return the type of the elements contained in `obj.data`.
"""
Base.eltype(::Type{GreenNew{T,N,MT}}) where {T,N,MT} = T

"""
    getindex(obj::GreenNew, inds...)

Return a subset of `obj`'s data as specified by `inds`, where each `inds` may be, for example, an Int, an AbstractRange, or a Vector. 
"""
Base.getindex(obj::GreenNew{T,N,MT}, inds::Vararg{Int,N}) where {T,MT,N} = Base.getindex(obj.data, inds...)
# Base.getindex(obj::GreenNew, I::Int) = Base.getindex(obj.data, I)

"""
    setindex!(obj::GreenNew, v, inds...)
    obj[inds...] = v

Store values from array `v` within some subset of `obj.data` as specified by `inds`.
"""
Base.setindex!(obj::GreenNew{T,N,MT}, v, inds::Vararg{Int,N}) where {T,MT,N} = Base.setindex!(obj.data, v, inds...)
# Base.setindex!(obj::GreenNew, v, I::Int) = Base.setindex!(obj.data, v, I)

# IndexStyle(::Type{<:GreenNew}) = IndexCartesian() # by default, it is IndexCartesian

"""
    Base.similar(obj::GreenNew{T,N,MT}, ::Type{S}) where {T,MT,N,S}
    Base.similar(obj::GreenNew{T,N,MT}) where {T,MT,N} = Base.similar(obj, T)

# Return type:
- `Base.similar(obj::GreenNew)`: Return a new GreenNew with the same meshes, and the uninitialized data of the same type as `obj.data`.
- `Base.similar(obj::GreenNew, ::Type{S})`: Return a new GreenNew with the same meshes, but with the uninitialized data of type `S`.
"""
function Base.similar(obj::GreenNew{T,N,MT}, ::Type{S}) where {T,MT,N,S}
    return GreenNew(obj.mesh...; dtype=S, data=similar(obj.data, S))
end
Base.similar(obj::GreenNew{T,N,MT}) where {T,MT,N} = Base.similar(obj, T)
#By default, the following functions will all call Base.similar(obj::GreenNew, ::Type{S}, inds) as explained in https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-array
#`Base.similar(obj::GreenNew, ::Type{S}, inds)`: Return a slice of obj.data.
#However, we don't want that since slice of GreeNew itself is not well defined with meshes.

################################ broadcast interface ###############################################
Base.BroadcastStyle(::Type{<:GreenNew}) = Broadcast.ArrayStyle{GreenNew}()

function Base.similar(bc::Base.Broadcast.Broadcasted{Broadcast.ArrayStyle{GreenNew}}, ::Type{ElType}) where {ElType}
    # println("get called")
    # Scan the inputs for the GreenNew:
    A = find_gf(bc)
    # Use other fields of A to create the output
    GreenNew(A.mesh, similar(Array{ElType}, axes(bc)), A.dims)
end

find_gf(bc::Broadcast.Broadcasted) = find_gf(bc.args)
find_gf(args::Tuple) = find_gf(find_gf(args[1]), Base.tail(args))
find_gf(x) = x
find_gf(::Tuple{}) = nothing
find_gf(a::GreenNew, rest) = a
find_gf(::Any, rest) = find_gf(rest)

function Base.copyto!(dest, bc::Base.Broadcast.Broadcasted{GreenNew{T,N,MT}}) where {T,MT,N}
    # without this function, inplace operation like g1 .+= g2 will make a lot of allocations
    # Please refer to the following posts for more details:
    # 1. manual on the interface: https://docs.julialang.org/en/v1/manual/interfaces/#extending-in-place-broadcast-2
    # 2. see the post: https://discourse.julialang.org/t/help-implementing-copyto-for-broadcasting/51204/3
    # 3. example from DataFrames.jl: https://github.com/JuliaData/DataFrames.jl/blob/main/src/other/broadcasting.jl#L193

    ######## approach 2: use materialize ########
    bcf = Base.Broadcast.materialize(bc)
    for I in CartesianIndices(dest)
        dest[I] = bcf[I]
    end
    return dest
end

########### alternative approach ######################
# function Base.copyto!(dest::GreenNew{T, N, MT}, bc::Base.Broadcast.Broadcasted{Nothing}) where {T,MT,N}
#     _bcf = Base.Broadcast.flatten(bc)
#     bcf = Base.Broadcast.preprocess(dest, _bcf)
#     for I in CartesianIndices(dest)
#         dest[I] = bcf[I]
#     end
#     return dest
# end


# somehow, the following leads to stackoverflow due to some kind of infinite loop
# function Base.getproperty(obj::GreenNew{T,MT,N,Ninner}, sym::Symbol) where {T,MT,N,Ninner}
#     if sym === :N
#         return N
#     elseif sym === :Ninner
#         return Ninner
#     elseif sym === :dims
#         return obj.dims
#     else # fallback to getfield
#         return getfield(obj, sym)
#     end
# end

"""
    show(io::IO, obj::GreenNew)

Write a text representation of the Green's function `obj` to the output stream `io`.
"""
function Base.show(io::IO, obj::GreenNew)
    print(io, "Green's function with dims = $(obj.dims) and total length = $(length(obj.data))\n"
              *
              "- Mesh: $(typeof(obj.mesh)) \n"
    )
end

"""
    function rank(obj::GreenNew{T,N,MT})

Return the dimension of `obj.data` (`N`).
"""
rank(::Type{GreenNew{T,N,MT}}) where {T,MT,N} = N


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
    function _check(objL::GreenNew, objR::GreenNew)

Check if the Green's functions `objL` and `objR` are on the same meshes. Throw an AssertionError if any check is false.
"""
function _check(objL::GreenNew, objR::GreenNew)
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

 """
    <<(Obj::GreenDLR, objSrc::Py)
    Converts the green function from triqs to GreenNew.
"""
    
# """
#     <<(Obj::GreenDLR, objSrc::Expr)
#     Obj << objSrc

# Initiate the Green's function `Obj` with the given function expression `objSrc`.
# TODO: Add other behaviors including:
# 1. Convert other type of GreenFunc to GreenDLR
# 2. Assign a GreenDLR to another GreenDLR
# 3. Convert triqs Green's function to GreenDLR
# 4. First convert triqs object to triqs green's function, than to GreenDLR
# """
# function Base.:<<(Obj::GreenDLR, objSrc::Expr)
#     # init version of <<
#     # more general version needed
#     for (id, d) in enumerate(Obj)
#         inds = ind2sub_gen(size(Obj), id)
#         p, ωn, n, τ = NaN, NaN, NaN, NaN
#         G = d
#         β = Obj.β
#         if Obj.domain == ImFreq
#             n = Obj.tgrid[inds[3]]
#             if Obj.isFermi
#                 ωn = π * (2 * n + 1) / β
#             else
#                 ωn = π * 2 * n * β
#             end
#         elseif Obj.domain == ImTime
#             τ = tgrid[inds[3]]
#         end
#         p = Obj.mesh[inds[2]]

#         m = Dict(
#             :G => G,
#             :p => p,
#             :ωn => ωn,
#             :n => n,
#             :τ => τ,
#             :β => β
#         )
#         Obj[id] = DictParser.evalwithdict(objSrc, m)
#     end

#     return nothing
# end

# Return the single-particle density matrix of the Green's function `obj`.
# """
# function density(obj::GreenDLR; kwargs...)
#     G_ins = toTau(obj, [obj.β,]).data .* (-1)
#     return selectdim(G_ins, ndims(G_ins), 1)
# end

#rMatrix transform of the target space of a matrix valued Greens function.
#Sets the current Greens function :math:`g_{ab}` to the matrix transform of :math:`G_{cd}`
#using the left and right transform matrices :math:`L_{ac}` and :math:`R_{db}`.
#.. math::
#g_{ab} = \sum_{cd} L_{ac} G_{cd} R_{db}


# """
#     def from_L_G_R(self, L, G, R):
#         Parameters
#         ----------
#         L : (a, c) ndarray
#             Left side transform matrix.
#         G : Gf matrix valued target_shape == (c, d)
#             Greens function to transform.
#         R : (d, b) ndarray
#             Right side transform matrix.
#         Notes
#         -----
#         Only implemented for Greens functions with a single mesh.


#         assert self.rank == 1, "Only implemented for Greens functions with one mesh"
#         assert self.target_rank == 2, "Matrix transform only valid for matrix valued Greens functions"

#         assert len(L.shape) == 2, "L needs to be two dimensional"
#         assert len(R.shape) == 2, "R needs to be two dimensional"

#         assert L.shape[1] == G.target_shape[0], "Dimension mismatch between L and G"
#         assert R.shape[0] == G.target_shape[1], "Dimension mismatch between G and R"

#         assert L.shape[0] == self.target_shape[0], "Dimension mismatch between L and self"
#         assert R.shape[1] == self.target_shape[1], "Dimension mismatch between R and self"

#         if not L.strides == sorted(L.strides):
#             L = L.copy(order='C')

#         if not R.strides == sorted(R.strides):
#             R = R.copy(order='C')

#         wrapped_aux.set_from_gf_data_mul_LR(self.data, L, G.data, R)
# """
# function from_L_G_R(self,L,G::GreenDLR,R)
#     return 1
# end