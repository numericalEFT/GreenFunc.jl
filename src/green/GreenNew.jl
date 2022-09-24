# abstract type AbstractGreen{T,N,NINNER} <: AbstractArray{T,N} end

isiterable(::Type{T}) where {T} = hasmethod(iterate, (T,))
# isiterable(T) = hasmethod(iterate, (typeof(T),))

# """
#     mutable struct GreenDLR{T,Domain<:TimeDomain,TGT,MT,Ndata}

# General Green's function on a multi-dimensional mesh plus one in-built Discrete Lehmann Representation.

# # Parameters:
# - `T`: type of data
# - `Domain`: type of time domain, `Domain`<:`TimeDomain`.
# - `TGT`: type of time grid
# - `MT`: type of mesh
# - `N`: number of internal degrees of freedom
# - `Ndata`: rank of Green's function data, which always equals to N+2, 2 stands for the mesh and the extra dimension that has built-in DLR grid.

# # Members:
# - `DLR`: built-in DLR grid. Only one-dimensional DLR is available currently.
# - `tgrid` (TGT): the imaginary-time or Matsubara-frequency grid of dimension with built in DLR . If not provided by user, the optimized grid from DLR is used.
# - `mesh` (MT): the mesh is a direct product of grids of all other continuous degrees of freedom of Green's function, other than the one with DLR. The mesh has to support all standard Base functions of AbstractArray, plus the following two:
#     - locate(`mesh`, value): find the index of the closest grid point for given value;
#     - volume(`mesh`, index): find the volume of grid space near the point at griven index.
#     - volume(`mesh`, gridvalue): locate the corresponding index of a given value and than find the volume of grid space. 
# - `innerstate` (Tuple): innerstate saves the discrete inner dgrees of freedom of Green's function. 
# - `data` (Array{T,Ndata}): the data of the Green's function.
# """
# mutable struct GreenNew{T,MT,N} <: AbstractGreen{T,N}
mutable struct GreenNew{T,N,MT} <: AbstractArray{T,N}
    #########   Mesh   ##############
    mesh::MT
    data::Array{T,N}
    dims::NTuple{N,Int}
end

"""
    function GreenNew{T}(mesh...;
        innerstate::Union{AbstractVector{Int},Tuple{Vararg{Int}}}=(),
        data::Union{Nothing,AbstractArray}=nothing) where {T}
    
Create a Green struct. 

# Arguments
- `T`: data type of Green's function's value.
- `mesh`: meshes of Green's function. Mesh could be any iterable object, examples are vector, tuple, array, number, UnitRange (say, 1:5).
- `innerstate`: innerstate saves the discrete inner dgrees of freedom of Green's function. By default, `innerstate` = (1,).
- `data`: the data of the Green's function. By default, `data` = zeros(`datatype`, `Ndata`).
"""
function GreenNew(mesh...;
    dtype=Float64,
    data::Union{Nothing,AbstractArray}=nothing)

    @assert all(x -> isiterable(typeof(x)), mesh) "all meshes should be iterable"

    N = length(mesh)
    dims = tuple([length(v) for v in mesh]...)

    if isnothing(data) == false
        @assert length(size(data)) == N
    else
        data = Array{dtype,N}(undef, dims...)
    end

    return GreenNew{dtype,N,typeof(mesh)}(mesh, data, dims)
end

########## Array Interface: https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-array #############
Base.size(obj::GreenNew) = obj.dims

Base.eltype(::Type{GreenNew{T,N,MT}}) where {T,N,MT} = T

"""
    getindex(obj::GreenDLR, inds...)

Return a subset of `obj`'s data as specified by `inds`, where each `inds` may be, for example, an Int, an AbstractRange, or a Vector. 
"""
Base.getindex(obj::GreenNew{T,N,MT}, inds::Vararg{Int,N}) where {T,MT,N} = Base.getindex(obj.data, inds...)
# Base.getindex(obj::GreenNew, I::Int) = Base.getindex(obj.data, I)
Base.setindex!(obj::GreenNew{T,N,MT}, v, inds::Vararg{Int,N}) where {T,MT,N} = Base.setindex!(obj.data, v, inds...)
# Base.setindex!(obj::GreenNew, v, I::Int) = Base.setindex!(obj.data, v, I)

# IndexStyle(::Type{<:GreenNew}) = IndexCartesian() # by default, it is IndexCartesian

"""
    Base.similar(obj::GreenNew{T,MT,N,Ninner}, ::Type{S}) where {T,MT,N,Ninner,S}
    Base.similar(obj::GreenNew{T,MT,N,Ninner}) where {T,MT,N,Ninner} = Base.similar(obj, T)

# Return type:
- `Base.similar(obj::GreenNew)`: Return a new GreenNew with the same mesh and innerstate, and a new data of the same type as obj.data.
- `Base.similar(obj::GreenNew, ::Type{S})`: Return a new GreenNew with the same mesh and innerstate, but with a new data of type S.
- `Base.similar(obj::GreenNew, ::Type{S}, inds)`: Return a slice of obj.data. (slice of GreeNew itself is not well defined, because it has meshes) 
"""
function Base.similar(obj::GreenNew{T,N,MT}, ::Type{S}) where {T,MT,N,S}
    return GreenNew(obj.mesh...; dtype=S, data=similar(obj.data, S))
end
Base.similar(obj::GreenNew{T,N,MT}) where {T,MT,N} = Base.similar(obj, T)
#By default, the following functions will all call Base.similar(obj::GreenNew, ::Type{S}, inds) 
#as explained in https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-array
#However, we don't want that. What we want the following: 

################################ broadcast interface ###############################################
Base.BroadcastStyle(::Type{<:GreenNew}) = Broadcast.ArrayStyle{GreenNew}()

function Base.similar(bc::Base.Broadcast.Broadcasted{Broadcast.ArrayStyle{GreenNew}}, ::Type{ElType}) where {ElType}
    # println("get called")
    # Scan the inputs for the ArrayAndChar:
    A = find_gf(bc)
    # Use the char field of A to create the output
    GreenNew(A.mesh, similar(Array{ElType}, axes(bc)), A.dims)
end

find_gf(bc::Broadcast.Broadcasted) = find_gf(bc.args)
find_gf(args::Tuple) = find_gf(find_gf(args[1]), Base.tail(args))
find_gf(x) = x
find_gf(::Tuple{}) = nothing
find_gf(a::GreenNew, rest) = a
find_gf(::Any, rest) = find_gf(rest)

function Base.copyto!(dest::GreenNew{T,N,MT}, bc::Base.Broadcast.Broadcasted{Nothing}) where {T,MT,N}
    # don't why, but this function make no allocations anymore
    # see the post https://discourse.julialang.org/t/help-implementing-copyto-for-broadcasting/51204/3
    bcf = Base.Broadcast.flatten(bc)
    bcf2 = Base.Broadcast.preprocess(dest, bcf)
    for I in CartesianIndices(dest)
        dest[I] = bcf2[I]
    end
    return dest
end


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
    show(io::IO, obj::GreenDLR)

Write a text representation of the Green's function `obj` to the output stream `io`.
"""
function Base.show(io::IO, obj::GreenNew)
    print(io, "Green's function with dims = $(obj.dims) and innerstate = $(obj.innerstate), total length = $(length(obj.data))\n"
              *
              "- Mesh: $(typeof(obj.mesh)), \n"
    )
end

# """
#     similar(obj::GreenDLR)

# Create a data-uninitialized GreenDLR with the element type and size, based upon the given `obj`.
# Note that the elements `innerstate`, `tgrid`, `mesh`, and `DLR` are copied from `obj`.
# """
# function Base.similar(obj::GreenDLR)
#     new = GreenDLR(obj.β)
#     new.innerstate = obj.innerstate
#     new.tgrid = obj.tgrid
#     new.mesh = obj.mesh
#     new.DLR = obj.DLR
#     new.data = similar(obj.data)
#     return new
# end

"""
    function rank(obj::GreenDLR)

Return the rank of Green's function data, which always equals to N+2. 
N stands for the number of internal degrees of freedom; 2 stands for the mesh and the extra dimension that has built-in DLR grid.
"""
rank(obj::GreenNew{T,N,MT}) where {T,MT,N} = N


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
    function _check(objL::GreenDLR, objR::GreenDLR)

Check if the Green's functions `objL` and `objR` are on the same `innerstate`, `tgrid`, and `mesh`. Throw an AssertionError if any check is false.
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


# """
#     function toTau(obj::Green2DLR, targetGrid =  obj.DLR.τ)
# Convert Green's function to τ space by Fourier transform.
# If green is already in τ space then it will be interpolated to the new grid.

# #Arguements
# - 'green': Original Green's function
# - 'targetGrid': Grid of outcome Green's function. Default: DLR τ grid
# """
# function toimtime(obj::Green2DLR, targetGrid = obj.DLR.τ)

#     if targetGrid isa AbstractVector
#         targetGrid = CompositeGrids.SimpleG.Arbitrary{eltype(targetGrid)}(targetGrid)
#     end

#     # do nothing if the domain and the grid remain the same
#     if obj.domain == ImTime && length(obj.tgrid.grid) ≈ length(targetGrid.grid) && obj.tgrid.grid ≈ targetGrid.grid
#         return obj
#     end
#     if isempty(obj.dynamic) # if dynamic data has not yet been initialized, there is nothing to do
#         return obj
#     end


#     if (obj.domain == ImTime)
#         dynamic = tau2tau(obj.DLR, obj.dynamic, targetGrid.grid, obj.tgrid.grid; axis = 4)
#     elseif (obj.domain == ImFreq)
#         dynamic = matfreq2tau(obj.DLR, obj.dynamic, targetGrid.grid, obj.tgrid.grid; axis = 4)
#     elseif (obj.domain == DLRFreq)
#         dynamic = dlr2tau(obj.DLR, obj.dynamic, targetGrid.grid; axis = 4)
#     end

#     return Green2DLR{eltype(dynamic)}(
#         green.name, IMTIME, green.β, green.isFermi, green.DLR.Euv, green.spaceGrid, green.color;
#         tsym = green.tsym, tgrid = targetGrid, rtol = green.DLR.rtol,
#         dynamic = dynamic, instant = green.instant)
# end

# """
#     function toMatFreq(green::Green2DLR, targetGrid =  green.DLR.n)
# Convert Green's function to matfreq space by Fourier transform.
# If green is already in matfreq space then it will be interpolated to the new grid.

# #Arguements
# - 'green': Original Green's function
# - 'targetGrid': Grid of outcome Green's function. Default: DLR n grid
# """
# function toimfreq(green::Green2DLR, targetGrid = green.DLR.n)

#     if targetGrid isa AbstractVector
#         targetGrid = CompositeGrids.SimpleG.Arbitrary{eltype(targetGrid)}(targetGrid)
#     end

#     # do nothing if the domain and the grid remain the same
#     if green.domain == ImFreq && length(green.tgrid.grid) ≈ length(targetGrid.grid) && green.tgrid.grid ≈ targetGrid.grid
#         return green
#     end
#     if isempty(green.dynamic) # if dynamic data has not yet been initialized, there is nothing to do
#         return green
#     end


#     if (green.domain == ImFreq)
#         dynamic = matfreq2matfreq(green.DLR, green.dynamic, targetGrid.grid, green.tgrid.grid; axis = 4)
#     elseif (green.domain == ImTime)
#         dynamic = tau2matfreq(green.DLR, green.dynamic, targetGrid.grid, green.tgrid.grid; axis = 4)
#     elseif (green.domain == DLRFreq)
#         dynamic = dlr2matfreq(green.DLR, green.dynamic, targetGrid.grid; axis = 4)
#     end

#     return Green2DLR{eltype(dynamic)}(
#         green.name, IMFREQ, green.β, green.isFermi, green.DLR.Euv, green.spaceGrid, green.color;
#         tsym = green.tsym, tgrid = targetGrid, rtol = green.DLR.rtol,
#         dynamic = dynamic, instant = green.instant)

# end

# """
#     function toDLR(green::Green2DLR)
# Convert Green's function to dlr space.

# #Arguements
# - 'green': Original Green's function
# """
# function toDLR(green::Green2DLR)

#     # do nothing if the domain and the grid remain the same
#     if green.domain == DLRFreq
#         return green
#     end
#     if isempty(green.dynamic) # if dynamic data has not yet been initialized, there is nothing to do
#         return green
#     end


#     if (green.domain == ImTime)
#         dynamic = tau2dlr(green.DLR, green.dynamic, green.tgrid.grid; axis = 4)
#     elseif (green.domain == ImFreq)
#         dynamic = matfreq2dlr(green.DLR, green.dynamic, green.tgrid.grid; axis = 4)
#     end

#     return Green2DLR{eltype(dynamic)}(
#         green.name, DLRFREQ, green.β, green.isFermi, green.DLR.Euv, green.spaceGrid, green.color;
#         tsym = green.tsym, tgrid = green.DLR.ω, rtol = green.DLR.rtol,
#         dynamic = dynamic, instant = green.instant)

# end



# """
#     function dynamic(green::Union{Green2DLR{DT,Domain,TGT,MT},GreenSym2DLR{DT,Domain,TGT,MT}}, time, space, color1::Int, color2::Int, timeMethod::TM , spaceMethod::SM) where {DT,Domain,TGT<:CompositeGrids.AbstractGrid,MT<:CompositeGrids.AbstractGrid,TM,SM}
# Find value of Green's function's dynamic part at given color and k/x by interpolation.
# Interpolation method is by default depending on the grid, but could also be chosen to be linear.
# #Argument
# - 'green': Green's function
# - 'time': Target τ/ω_n point
# - 'space': Target k/x point
# - 'color1': Target color1
# - 'color2': Target color2
# - 'timeMethod': Method of interpolation for time
# - 'spaceMethod': Method of interpolation for space 
# """
# function _interpolation(TIM, SIM; green::Union{Green2DLR{DT,Domain,TGT,MT},GreenSym2DLR{DT,Domain,TGT,MT}}, time, space, color1::Int, color2::Int) where {DT,Domain,TGT<:CompositeGrids.AbstractGrid,MT<:CompositeGrids.AbstractGrid}
#     # for double composite
#     if isempty(green.data)
#         error("Dynamic Green's function can not be empty!")
#     else
#         spaceNeighbor = CompositeGrids.Interp.findneighbor(SIM, green.spaceGrid, space)
#         println(TIM)
#         if green.domain == ImFreq && TIM != DLRInterp
#             tgrid = (green.tgrid.grid * 2 .+ 1) * π / green.β
#             comTimeGrid = CompositeGrids.SimpleG.Arbitrary{eltype(tgrid)}(tgrid)
#             comTime = (2*time+1)*π/green.β
#         else
#             comTimeGrid = green.tgrid
#             comTime = time
#         end

#         timeNeighbor = CompositeGrids.Interp.findneighbor(TIM, comTimeGrid, comTime)
#         data_slice = view(green.data, color1, color2, spaceNeighbor.index, timeNeighbor.index)
#         data_slice_xint = CompositeGrids.Interp.interpsliced(spaceNeighbor,data_slice, axis=1)
#         result = CompositeGrids.Interp.interpsliced(timeNeighbor,data_slice_xint, axis=1)
#     end
#     return result
# end

# function _interpolation( ::LinearInterp , ::LinearInterp
#     ;green::Union{Green2DLR{DT,Domain,TGT,MT},GreenSym2DLR{DT,Domain,TGT,MT}}, time, space, color1::Int, color2::Int,
#     ) where {DT,Domain,TGT<:CompositeGrids.AbstractGrid,MT<:CompositeGrids.AbstractGrid}
#     # for double composite and double linear
#     if isempty(green.data)
#         error("Dynamic Green's function can not be empty!")
#     else
#         if green.domain == ImFreq
#             tgrid = (green.tgrid.grid * 2 .+ 1) * π / green.β
#             comTimeGrid = CompositeGrids.SimpleG.Arbitrary{eltype(tgrid)}(tgrid)            
#             comTime = (2*time+1)*π/green.β
#         else
#             comTimeGrid = green.tgrid
#             comTime = time
#         end
#         data_slice = view(green.data, color1, color2, :,:)
#         result = CompositeGrids.Interp.linear2D(data_slice, green.spaceGrid, comTimeGrid,space,comTime)
#     end
#     return result
# end


# function _interpolation(::DLRInterp, SIM; green::Union{Green2DLR{DT,Domain,TGT,MT},GreenSym2DLR{DT,Domain,TGT,MT}}, time, space, color1::Int, color2::Int) where {DT,Domain,TGT<:CompositeGrids.AbstractGrid,MT<:CompositeGrids.AbstractGrid}
#     # for composite space and dlr time
#     if isempty(green.data)
#         error("Dynamic Green's function can not be empty!")
#     else
#         spaceNeighbor = CompositeGrids.Interp.findneighbor(SIM, green.spaceGrid, space)
#         data_slice = view(green.data, color1, color2, spaceNeighbor.index,:)
#         data_slice_xint = CompositeGrids.Interp.interpsliced(spaceNeighbor,data_slice, axis=1)
#         if green.domain == ImFreq
#             result = (matfreq2matfreq(green.DLR, data_slice_xint, [time,], green.tgrid.grid))[1]
#         elseif green.domain == ImTime
#             result = (tau2tau(green.DLR, data_slice_xint, [time,], green.tgrid.grid))[1]
#         end
#     end
#     return result
# end

# interpolation(;timeMethod::TM, spaceMethod::SM ,green::Union{Green2DLR{DT,Domain,TGT,MT},GreenSym2DLR{DT,Domain,TGT,MT}}, time, space , color1::Int, color2::Int) where {TM,SM,TGT<:CompositeGrids.AbstractGrid,MT<:CompositeGrids.AbstractGrid,DT,Domain} = _interpolation(InterpMethod(TGT,TM), InterpMethod(MT, SM); green, time, space, color1, color2)

# function interpolation(green::Union{Green2DLR{DT,Domain,TGT,MT},GreenSym2DLR{DT,Domain,TGT,MT}}, time, space, color1::Int=1, color2::Int=color1, timeMethod::TM=DEFAULTINTERP , spaceMethod::SM=DEFAULTINTERP) where {DT,Domain,TGT<:CompositeGrids.AbstractGrid,MT<:CompositeGrids.AbstractGrid,TM,SM}
#     return  interpolation(; timeMethod = timeMethod, spaceMethod = spaceMethod, green = green, time=time, space=space, color1=color1, color2 =color2)
# end

# function interpolation(green::Union{Green2DLR{DT,Domain,TGT,MT},GreenSym2DLR{DT,Domain,TGT,MT}}, time, space, timeMethod::TM, spaceMethod::SM) where {DT,Domain,TGT<:CompositeGrids.AbstractGrid,MT<:CompositeGrids.AbstractGrid,TM,SM}
#     return  interpolation(; timeMethod = timeMethod, spaceMethod = spaceMethod, green = green, time=time, space=space, color1=1, color2 =1)
# end
