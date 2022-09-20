"""
General Green's function. 
"""

"""
Green's function on a multi-dimensional mesh plus one in-built Discrete Lehmann Representation.
#Parameters:
- 'T': type of data
- 'Domain': type of time domain, Domain<:TimeDomain
- 'TGT': type of time grid
- 'MT': type of mesh
-'N': number of internal degrees of freedom
-'Ndata': rank of Green's function data, which always equals to N+2, 2 stands for the mesh and the extra dimension that has built-in DLR grid.

#Members:
- 'DLR': built-in DLR grid. Only one-dimensional DLR is available currently.
- 'tgrid': the imaginary time/ Matsubara frequency grid of dimension with built in DLR . If not provided by user, the optimized grid from DLR is used.
- 'mesh': the mesh is a direct product of grids of all other continuous dgrees of freedom of Green's function, other than the one with DLR. The mesh has to support all standard Base functions of AbstractArray, plus the following two:
    - locate(mesh,value): find the index of the closest grid point for given value;
    - volume(mesh,index): find the volume of grid space near the point at griven index.
    - volume(mesh,gridvalue): locate the corresponding index of a given value and than find the volume of grid space. 
- 'innerstate': innerstate saves the discrete inner dgrees of freedom of Green's function. 
- 'data': the data of Green's function.
"""
mutable struct GreenDLR{T,Domain<:TimeDomain,TGT,MT,Ndata}
    DLR::DLRGrid

    #########   Mesh   ##############

    tgrid::TGT
    mesh::MT
    innerstate::Tuple
    ###########     data   ###########
    data::Array{T,Ndata}

    """

    """
    function GreenDLR{T}(; domain::Domain, DLR, tgrid::TGT, mesh::MT, β, isFermi, Euv, rtol, tsym, innerstate::Tuple, data::Array{T,Ndata}) where {T,Domain<:TimeDomain,TGT,MT,Ndata}
        @assert tsym == :ph || tsym == :pha || tsym == :none        # data is initilized with zeros if not provided by user
        @assert Ndata == length(innerstate) + 2

        gnew = new{T,Domain,typeof(tgrid),typeof(mesh),Ndata}(
            DLR,
            tgrid,
            mesh, innerstate,
            data)
        return gnew
    end
end

function GreenDLR(; kwargs...)
    if :datatype in keys(kwargs)
        datatype = kwargs[:datatype]
    else
        datatype = ComplexF64
    end

    if :domain in keys(kwargs)
        domain = kwargs[:domain]
    else
        domain = IMFREQ
    end
    #TODO: after support from mesh is done, the default value of mesh has to be a structure with locate and volume function instead of an array.
    if :mesh in keys(kwargs)
        mesh = kwargs[:mesh]
    else
        mesh = [0.0, 1.0]
    end
    #TODO: in principle \beta and EUV should be in proper energy unit of EF, which comes from parameter. Should we make it globel?
    if :β in keys(kwargs)
        β = kwargs[:β]
    else
        β = 100.0
    end
    if :isFermi in keys(kwargs)
        isFermi = kwargs[:isFermi]
    else
        isFermi = true
    end
    if :Euv in keys(kwargs)
        Euv = kwargs[:Euv]
    else
        Euv = 100.0
    end
    if :rtol in keys(kwargs)
        rtol = kwargs[:rtol]
    else
        rtol = 1e-10
    end

    if :tsym in keys(kwargs)
        tsym = kwargs[:tsym]
    else
        tsym = :none
    end

    DLR = DLRGrid(Euv, β, rtol, isFermi, tsym)
    Domain = typeof(domain)
    if :tgrid in keys(kwargs)
        givenTimeGrid = kwargs[:tgrid]
        if givenTimeGrid isa AbstractVector
            tgrid = CompositeGrids.SimpleG.Arbitrary{eltype(givenTimeGrid)}(givenTimeGrid)
        elseif givenTimeGrid isa CompositeGrids.AbstractGrid
            tgrid = givenTimeGrid
        else
            error("Input tgrid has to be Vector or Composite grid")
        end
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

    if :innerstate in keys(kwargs)
        innerstate = kwargs[:innerstate]
    else
        innerstate = (1,)
    end
    if :data in keys(kwargs)
        data = kwargs[:data]
    else
        # N = length(innerstate)+2
        # data = Array{datatype, N}[]
        datasize = collect(innerstate)
        datasize = push!(datasize, length(mesh))
        datasize = push!(datasize, DLR.size)
        datasize = tuple(datasize...)

        data = zeros(datatype, datasize)
    end


    return GreenDLR{datatype}(; domain=domain, DLR=DLR, tgrid=tgrid, mesh=mesh, β=β, isFermi=isFermi, Euv=Euv, rtol=rtol, tsym=tsym, innerstate=innerstate, data=data)

end
"""


"""
Base.getindex(obj::GreenDLR, inds...) = Base.getindex(obj.data, inds...)


"""

"""
function Base.setindex!(obj::GreenDLR, X, inds...)
    obj.data[inds...] = X
end

Base.view(obj::GreenDLR, inds...) = Base.view(obj.data, inds...)

function Base.getproperty(obj::GreenDLR{T,Domain,TGT,MT}, sym::Symbol) where {T,Domain,TGT,MT}
    if sym === :isFermi
        return obj.DLR.isFermi
    elseif sym === :β
        return obj.DLR.β
    elseif sym === :domain
        return Domain
    elseif sym === :tsym
        return obj.DLR.symmetry
    else # fallback to getfield
        return getfield(obj, sym)
    end
end

function Base.size(obj::GreenDLR)
    return size(obj.data)
end
#TODO:nice print
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
              * "- Mesh: shape = $(size(obj.mesh)), length = $(length(obj.mesh))\n"
              * "- timeGrid: domain = $(obj.domain), symmetry = " * sym * ", length = $(size(obj.tgrid.grid)[1])\n"
    )
end

function rank(obj::GreenDLR)
    return length(obj.innerstate) + 2
end


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

function check(objL::GreenDLR, objR::GreenDLR)
    @assert objL.tgrid == objR.tgrid "Green's function time grids are not compatible:\n $(objL.tgrid)\nand\n $(objR.tgrid)"
    @assert objL.mesh == objR.mesh "Green's function meshes are not compatible:\n $(objL.mesh)\nand\n $(objR.mesh)"
    @assert objL.innerstate == objR.innerstate "Green's function innerstates are not compatible:\n $(objL.innerstate) and $(objR.innerstate)"
end

function Base.:<<(obj::GreenDLR, objSrc::GreenDLR)
    check(obj, objSrc)
    obj = deepcopy(objSrc)
    return obj
end

function Base.:<<(Obj::GreenDLR, objSrc::Expr)
    return 1
end

"""
TODO:Need to check objL and objR are on the same grid and has same innerstate
"""
function Base.:-(obj::GreenDLR)
    new = obj
    new.data = -new.data
    return new
end


function Base.:+(objL::GreenDLR, objR::GreenDLR)
    # check(objL, objR)
    new = objL
    new.data = objL.data + objR.data
    return new
end

function Base.:-(objL::GreenDLR, objR::GreenDLR)
    # check(objL, objR)
    new = objL
    new.data = objL.data - objR.data
    return new

end

function Base.:*(objL::GreenDLR, objR::GreenDLR)
    # check(objL, objR)
    new = objL
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

#Density calculation only supports homogeneous system now.
function density(obj::GreenDLR; kwargs...)
    G_ins = toTau(obj, [obj.β,]).data .* (-1)

    grid = zeros(length(obj.mesh))
    for (ki, ind) in enumerate(linear_to_index(obj.mesh))
        grid[ki] = norm(getindex(obj.mesh))
    end
    kgrid = CompositeGrids.AbstractGrid(grid)

    d = rank(obj.mesh)
    Ωd = d == 2 ? 1 : √π / 2
    integrand = real(G_ins) .* kgrid.grid .^ d

    return CompositeGrids.Interp.integrate1D(integrand, kgrid) * Ωd / 2 / (2π)^d
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
