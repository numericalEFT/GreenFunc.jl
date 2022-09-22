include("./DictParser.jl")
# using .DictParser

"""
    mutable struct GreenDLR{T,Domain<:TimeDomain,TGT,MT,Ndata}

General Green's function on a multi-dimensional mesh plus one in-built Discrete Lehmann Representation.

# Parameters
- `T`: type of data
- `Domain`: type of time domain, `Domain`<:`TimeDomain`.
- `TGT`: type of time grid
- `MT`: type of mesh
- `N`: number of internal degrees of freedom
- `Ndata`: rank of Green's function data, which always equals to N+2, 2 stands for the mesh and the extra dimension that has built-in DLR grid.

# Members
- `DLR`: built-in DLR grid. Only one-dimensional DLR is available currently.
- `domain` (Domain): domain of time grid.
- `tgrid` (TGT): the imaginary-time or Matsubara-frequency grid of dimension with built in DLR . If not provided by user, the optimized grid from DLR is used.
- `mesh` (MT): the mesh is a direct product of grids of all other continuous degrees of freedom of Green's function, other than the one with DLR. The mesh has to support all standard Base functions of AbstractArray, plus the following two:
    - locate(`mesh`, value): find the index of the closest grid point for given value;
    - volume(`mesh`, index): find the volume of grid space near the point at griven index.
    - volume(`mesh`, gridvalue): locate the corresponding index of a given value and than find the volume of grid space. 
- `innerstate` (Tuple): innerstate saves the discrete inner dgrees of freedom of Green's function. 
- `data` (Array{T,Ndata}): the data of the Green's function.
"""
mutable struct GreenDLR{T,Domain<:TimeDomain,TGT,MT,Ninner,Ndata}
    DLR::DLRGrid
    domain::Domain

    #########   Mesh   ##############

    tgrid::TGT
    mesh::MT
    innerstate::NTuple{Ninner,Int}
    ###########     data   ###########
    data::Array{T,Ndata}

    function GreenDLR{T}(; domain::Domain, DLR::DLRGrid, tgrid::TGT, mesh::MT, innerstate, data::Array{T,Ndata}) where {T,Domain<:TimeDomain,TGT,MT,Ndata}
        @assert Ndata == length(innerstate) + 2 "ndims of data must be equal to length(innerstate) + 2"
        @assert Tuple(size(data)[1:length(innerstate)]) == innerstate
        @assert size(data)[end] == length(tgrid)
        @assert size(data)[end-1] == length(mesh)
        if tgrid isa AbstractVector
            tgrid = CompositeGrids.SimpleG.Arbitrary{eltype(tgrid)}(tgrid)
        elseif !(tgrid isa CompositeGrids.AbstractGrid)
            error("Input tgrid has to be Vector or Composite grid")
        end

        gnew = new{T,Domain,typeof(tgrid),typeof(mesh),Ndata - 2,Ndata}(
            DLR,
            domain,
            tgrid,
            mesh, Tuple(innerstate),
            data)
        return gnew
    end
end

"""
    function GreenDLR(β::Real;
        domain::Union{ImTime,ImFreq}=ImFreq(), 
        datatype=(domain == ImTime() ? Float64 : ComplexF64), 
        mesh=[0.0, 1.0], isFermi=true, Euv=100.0, rtol=1e-10, tsym=:none, 
        innerstate=(1,),
        kwargs...
    )
    
Create a GreenDLR struct. 

# Arguments
- `β`: inverse temperature. Must be input by the user.
- `domain`: domain of time grid, Domain<:TimeDomain. By default, `domain` = IMFREQ.
- `datatype`: data type of Green's function's value. By default, `datatype` = (domain == IMTIME ? Float64 : ComplexF64)`.
- `mesh`: direct product of grids of all other continuous degrees of freedom of Green's function. By default, `mesh` ...
- `isFermi` (Bool): is the particle fermion or boson? By default, `isFermi` = true.
- `Euv`: the UV energy scale of the spectral density. By default, `Euv` = 100.0.
- `rtol`: tolerance absolute error. By default, `rtol` = 1e-10.
- `tsym`: whether the Green's function has particle-hole symmetry (:ph), anti-particle-hole symmetry (:pha) or none of them (:none). By default, `tsym` = :none.
- `innerstate`: innerstate saves the discrete inner dgrees of freedom of Green's function. By default, `innerstate` = (1,).

# Optional Arguments
- `tgrid`: time grid as a AbstractVector or CompositeGrids.AbstractGrid. By default, a optimized grid built in DLR is used.
- `data`: the data of the Green's function. By default, `data` = zeros(`datatype`, `Ndata`).
"""
# function GreenDLR(; β::Real,
function GreenDLR(β::Real;
    domain::Union{ImTime,ImFreq}=ImFreq(), #....
    datatype=(domain == ImTime() ? Float64 : ComplexF64), #...
    mesh=[0.0, 1.0], isFermi=true, Euv=100.0, rtol=1e-10, tsym=:none, #...
    innerstate=(1,), #...
    kwargs...
)
    @assert tsym == :ph || tsym == :pha || tsym == :none "tsym must be :ph, :pha, or :none"
    # if :datatype in keys(kwargs)
    #     datatype = kwargs[:datatype]
    # else
    #     datatype = ComplexF64
    # end

    # if :domain in keys(kwargs)
    #     domain = kwargs[:domain]
    # else
    #     domain = IMFREQ
    # end
    #TODO: after support from mesh is done, the default value of mesh has to be a structure with locate and volume function instead of an array.
    # if :isFermi in keys(kwargs)
    #     isFermi = kwargs[:isFermi]
    # else
    #     isFermi = true
    # end
    # if :Euv in keys(kwargs)
    #     Euv = kwargs[:Euv]
    # else
    #     Euv = 100.0
    # end
    # if :rtol in keys(kwargs)
    #     rtol = kwargs[:rtol]
    # else
    #     rtol = 1e-10
    # end

    # if :tsym in keys(kwargs)
    #     tsym = kwargs[:tsym]
    # else
    #     tsym = :none
    # end

    DLR = DLRGrid(Euv, β, rtol, isFermi, tsym)
    Domain = typeof(domain)
    if :tgrid in keys(kwargs)
        tgrid = kwargs[:tgrid]
    else
        if Domain == ImFreq
            bareTimeGrid = DLR.n
        elseif Domain == ImTime
            bareTimeGrid = DLR.τ
        elseif Domain == DLRFreq
            bareTimeGrid = DLR.ω
        else
            error("$TimeType is not supported!")
        end
        tgrid = CompositeGrids.SimpleG.Arbitrary{eltype(bareTimeGrid)}(bareTimeGrid)
    end

    if :data in keys(kwargs)
        data = kwargs[:data]
    else
        # N = length(innerstate)+2
        # data = Array{datatype, N}[]
        datasize = collect(innerstate)
        datasize = push!(datasize, length(mesh))
        datasize = push!(datasize, length(tgrid))
        datasize = tuple(datasize...)

        data = zeros(datatype, datasize)    # data is initilized with zeros if not provided by user
    end


    return GreenDLR{datatype}(; domain=domain, DLR=DLR, tgrid=tgrid, mesh=mesh, innerstate=innerstate, data=data)

end

"""
    getindex(obj::GreenDLR, inds...)

Return a subset of `obj`'s data as specified by `inds`, where each `inds` may be, for example, an Int, an AbstractRange, or a Vector. 
"""
Base.getindex(obj::GreenDLR, inds...) = Base.getindex(obj.data, inds...)
Base.getindex(obj::GreenDLR, I::Int) = Base.getindex(obj.data, I)
Base.firstindex(obj::GreenDLR) = 1
Base.lastindex(obj::GreenDLR) = length(obj)

"""
    iterate(obj::GreenDLR, state)

Return a 2-tuple of the next element and the new iteration state. If no elements remain, `nothing` will be returned.
"""
Base.iterate(obj::GreenDLR) = (obj[1], 1)
Base.iterate(obj::GreenDLR, state) = (state >= length(obj)) ? nothing : (obj[state+1], state + 1)
# Base.IteratorSize(::Type{MeshProduct{MT,N}}) where {MT,N} = Base.HasLength()
# Base.IteratorEltype(::Type{MeshProduct{MT,N}}) where {MT,N} = Base.HasEltype()
# Base.eltype(::Type{MeshProduct{MT,N}}) where {MT,N} = tuple(eltype.(fieldtypes(MT))...) # 

@generated function sub2ind_gen(dims::NTuple{N}, I::Integer...) where {N}
    ex = :(I[$N] - 1)
    for i = (N-1):-1:1
        ex = :(I[$i] - 1 + dims[$i] * $ex)
    end
    return :($ex + 1)
end
@generated function ind2sub_gen(dims::NTuple{N}, I::Integer) where {N}
    inds, quotient = :((I - 1) % dims[1] + 1), :((I - 1) ÷ dims[1])
    for i = 2:N-1
        inds, quotient = :($inds..., $quotient % dims[$i] + 1), :($quotient ÷ dims[$i])
    end
    inds = :($inds..., $quotient + 1)
    return :($inds)
end

"""
    setindex!(obj::GreenDLR, X, inds...)
    obj[inds...] = X

Store values from array `X` within some subset of `obj.data` as specified by `inds`.
"""
function Base.setindex!(obj::GreenDLR, X, inds...)
    obj.data[inds...] = X
end
function Base.setindex!(obj::GreenDLR, X, I::Int)
    obj.data[I] = X
end

"""
    view(obj::GreenDLR, inds...)

Return a lightweight array that is effectively a _view_ into the parent array `obj.data` at the given index or indices `inds` instead of eagerly extracting elements or constructing a copied subset.
"""
Base.view(obj::GreenDLR, inds...) = Base.view(obj.data, inds...)

function Base.getproperty(obj::GreenDLR{T,Domain,TGT,MT}, sym::Symbol) where {T,Domain,TGT,MT}
    if sym === :isFermi
        return obj.DLR.isFermi
    elseif sym === :β
        return obj.DLR.β
        # elseif sym === :domain
        #     return Domain
    elseif sym === :tsym
        return obj.DLR.symmetry
    else # fallback to getfield
        return getfield(obj, sym)
    end
end

"""
    size(obj::GreenDLR)

Return a tuple containing the dimensions of `obj.data`.
"""
Base.size(obj::GreenDLR) = size(obj.data)

"""
    length(obj::GreenDLR)

Return the number of elements in `obj.data`.
"""
Base.length(obj::GreenDLR) = length(obj.data)

"""
    show(io::IO, obj::GreenDLR)

Write a text representation of the Green's function `obj` to the output stream `io`.
"""
function Base.show(io::IO, obj::GreenDLR)
    if obj.tsym == :ph
        sym = "particle-hole"
    elseif obj.tsym == :pha
        sym = "anti-particle-hole"
    else
        sym = "none"
    end
    print(io, (obj.isFermi ? "Fermionic " : "Bosonic ")
              * "Green's function with beta = $(obj.β) and innerstate = $(obj.innerstate) \n"
              * "- mesh: $(typeof(obj.mesh)), shape = $(size(obj.mesh)), length = $(length(obj.mesh))\n"
              * "- time grid: $(typeof(obj.tgrid)), domain = $(obj.domain), symmetry = " * sym * ", length = $(length(obj.tgrid))\n"
              * "- data: $(typeof(obj.data)), shape = $(size(obj))\n"
    )
end

"""
    similar(obj::GreenDLR)

Create a data-uninitialized GreenDLR with the element type and size, based upon the given `obj`.
Note that the elements `domain`, `innerstate`, `tgrid`, `mesh`, and `DLR` are copied from `obj`.
"""
function Base.similar(obj::GreenDLR)
    type = eltype(obj.data)
    return GreenDLR{type}(; domain=obj.domain, DLR=obj.DLR, tgrid=obj.tgrid, mesh=obj.mesh, innerstate=obj.innerstate, data=zeros(type, size(obj.data)))
end
# function Base.similar(obj::GreenDLR)
#     new = GreenDLR(β=obj.β)
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
rank(obj::GreenDLR) = length(obj.innerstate) + 2


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

Check if the Green's functions `objL` and `objR` are on the same `innerstate`, `domain`, `tgrid`, and `mesh`. Throw an AssertionError if any check is false.
"""
function _check(objL::GreenDLR, objR::GreenDLR)
    # KUN: check --> __check
    # first:  check typeof(objL.tgrid)==typeof(objR.tgrid) 
    # second: check length(objL.tgrid)
    # third:  hasmethod(objL.tgrid, isequal) --> assert
    @assert objL.innerstate == objR.innerstate "Green's function innerstates are not inconsistent: $(objL.innerstate) and $(objR.innerstate)"
    @assert objL.domain == objR.domain "Green's function domains are inconsistent: $(objL.domain) and $(objR.domain)"
    @assert typeof(objL.tgrid) == typeof(objR.tgrid) "Green's function time grids' types are inconsistent: $(typeof(objL.tgrid)) and $(typeof(objR.tgrid))"
    @assert typeof(objL.mesh) == typeof(objR.mesh) "Green's function meshes' types are inconsistent: $(typeof(objL.mesh)) and $(typeof(objR.mesh))"
    @assert objL.tgrid.size == objR.tgrid.size "Green's function time grids' length are inconsistent: $(objL.tgrid.size) and $(objR.tgrid.size)"
    @assert size(objL.mesh) == size(objR.mesh) "Green's function meshes' shapes are inconsistent: $(size(objL.mesh)) and $(size(objR.mesh))"

    # @assert objL.tgrid == objR.tgrid "Green's function time grids are not compatible:\n $(objL.tgrid)\nand\n $(objR.tgrid)"
    # @assert objL.mesh == objR.mesh "Green's function meshes are not compatible:\n $(objL.mesh)\nand\n $(objR.mesh)"
end

"""
TODO:Reload Base.Broadcast.broadcast(f,obj...)
Apply function f over the data of obj. There could be several situations:
obj...  includes
1. two green's function
2. green's function and a scalar
3. one Green's function
"""
Base.Broadcast.broadcast(f, obj::GreenDLR, I::Number) = Base.Broadcast.broadcast(f, obj.data, I::Number)
Base.Broadcast.broadcast(f, obj::GreenDLR) = Base.Broadcast.broadcast(f, obj.data)
Base.Broadcast.broadcast(f, objL::GreenDLR, objR::GreenDLR) = Base.Broadcast.broadcast(f, objL.data, objR.data)




"""
    <<(Obj::GreenDLR, objSrc::Expr)
    Obj << objSrc

Initiate the Green's function `Obj` with the given function expression `objSrc`.
TODO: Add other behaviors including:
1. Convert other type of GreenFunc to GreenDLR
2. Assign a GreenDLR to another GreenDLR
3. Convert triqs Green's function to GreenDLR
4. First convert triqs object to triqs green's function, than to GreenDLR
"""
function Base.:<<(Obj::GreenDLR, objSrc::Expr)
    # init version of <<
    # more general version needed
    for (id, d) in enumerate(Obj)
        inds = ind2sub_gen(size(Obj), id)
        p, ωn, n, τ = NaN, NaN, NaN, NaN
        G = d
        β = Obj.β
        if Obj.domain == IMFREQ
            n = Obj.tgrid[inds[3]]
            if Obj.isFermi
                ωn = π * (2 * n + 1) / β
            else
                ωn = π * 2 * n * β
            end
        elseif Obj.domain == IMTIME
            τ = tgrid[inds[3]]
        end
        p = Obj.mesh[inds[2]]

        m = Dict(
            :G => G,
            :p => p,
            :ωn => ωn,
            :n => n,
            :τ => τ,
            :β => β
        )
        Obj[id] = DictParser.evalwithdict(objSrc, m)
    end

    return nothing
end

"""
    -(obj::GreenDLR)

Map elements of `obj.data` to their additive inverses.
"""
function Base.:-(obj::GreenDLR)
    new = similar(obj)
    new.data = -obj.data
    return new
end

"""
    +(objL::GreenDLR, objR::GreenDLR)
    objL + objR

Perform addition between `objL.data` and `objR.data`.
"""
function Base.:+(objL::GreenDLR, objR::GreenDLR)
    _check(objL, objR)
    new = similar(objL)
    new.data = objL.data + objR.data
    return new
end

"""
    -(objL::GreenDLR, objR::GreenDLR)
    objL - objR

Perform subtraction between `objL.data` and `objR.data`.
"""
function Base.:-(objL::GreenDLR, objR::GreenDLR)
    _check(objL, objR)
    new = similar(objL)
    new.data = objL.data - objR.data
    return new

end

"""
    *(objL::GreenDLR, objR::GreenDLR)
    objL * objR

Perform multiplication between `objL.data` and `objR.data`.
"""
function Base.:*(objL::GreenDLR, objR::GreenDLR)
    _check(objL, objR)
    new = similar(objL)
    new.data = objL.data .* objR.data
    return new
end


#TODO:return density matrix of the Green's function

# """
# def density(self, *args, **kwargs):
#     rCompute the density matrix of the Greens function
#         Parameters
#         ----------
#         beta : float, optional
#             Used for finite temperature density calculation with ``MeshReFreq``.
#         Returns
#         -------
#         density_matrix : ndarray
#             Single particle density matrix with shape ``target_shape``.
#         Notes
#         -----
#         Only works for single mesh Greens functions with a, Matsubara,
#         real-frequency, or Legendre mesh.
# return gf_fnt.density(self, *args, **kwargs)
# """
"""
    function density(obj::GreenDLR; kwargs...)

Return the single-particle density matrix of the Green's function `obj`.
"""
function density(obj::GreenDLR; kwargs...)
    G_ins = toTau(obj, [obj.β,]).data .* (-1)
    return selectdim(G_ins, ndims(G_ins), 1)
end


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
