"""
    Base.:<<(objL::MeshArray, objR::MeshArray)

DLR Fourier transform of functions that has exactly one TemporalGrid(ImTime, ImFreq or DLRFreq) among the meshes. 
If objL and objR have identical TemporalGrid, objL<<objR assign objR to objL.
If objL and objR have different TemporalGrid, one of them has to be in DLR space.
If objL is in DLR space, objL<<objR calculates the DLR spectral density of data in objR
if objR is in DLR space, objL<<objR calculates the corresponding data from the DLR spectral density in objR.
"""

function Base.:<<(objL::MeshArray, objR::MeshArray)
    # init version of <<
    # more general version needed
    axes = []
    for (mi, mesh) in enumerate(objL.mesh)
        if (typeof(mesh) <: TemporalGrid)
            append!(axes, mi)
        else
            #TODO:add hashtable for mesh to support == operation 
            @assert typeof(mesh) == typeof(objR.mesh[mi]) "Meshes not involved in Fourier transform have to be identical" #should assert mesh == objR.mesh[mi] when == is defined
        end
    end
    @assert length(axes) == 1 "Only one temporal mesh with built in DLR grid is allowed"
    typeL = typeof(objL.mesh[axes[1]])
    typeR = typeof(objR.mesh[axes[1]])
    meshL = objL.mesh[axes[1]]
    meshR = objR.mesh[axes[1]]

    if typeL == typeR
        error("Assignment still work in progress") #need == operation for mesh
    # if meshL == meshR
    #     objL.data=deepcopy(objR.data)
    # else
    #     error("Green's function can only be assigned to another with the same mesh")
    # end
    elseif typeL <: MeshGrids.DLRFreq
        if typeR <: MeshGrids.ImFreq
            objL.data .= matfreq2dlr(meshL.dlr, objR.data, meshR.grid; axis=axes[1])
        elseif typeR <: MeshGrids.ImTime
            objL.data .= tau2dlr(meshL.dlr, objR.data, meshR.grid; axis=axes[1])
        end
    elseif typeR <: MeshGrids.DLRFreq
        if typeL <: MeshGrids.ImFreq
            objL.data .= dlr2matfreq(meshR.dlr, objR.data, meshL.grid; axis=axes[1])
        elseif typeL <: MeshGrids.ImTime
            objL.data .= dlr2tau(meshR.dlr, objR.data, meshL.grid; axis=axes[1])
        end
    else
        error("One of the Grren's function has to be in DLRfreq space to do Fourier transform")
    end
end


"""
    function dlr_to_imtime(obj::MeshArray, tgrid=nothing; dim=nothing)

Transform a Green's function in DLR space to the imaginary-time space. 
#Arguements
- 'obj': Function in DLR space
- `tgrid`: The imaginary-time grid which the function transforms into. Default value is the imaginary-time grid from the `DLRFreq` constructor.
- `dim`: The dimension of the temporal mesh. Default value is the first ImTime mesh.
"""
function dlr_to_imtime(obj::MeshArray{T,N,D,MT}, tgrid=nothing; dim::Union{Nothing,Int}=nothing) where {T,N,MT,D}
    ########################## generic interface #################################
    if isnothing(dim)
        dim = findfirst(x -> (x isa MeshGrids.DLRFreq), obj.mesh)
        @assert isnothing(dim) == false "No temporal can be transformed to imtime."
    end

    mesh = obj.mesh[dim]
    @assert mesh isa MeshGrids.DLRFreq "DLRFreq is expect for the dim = $dim."

    if tgrid isa MeshGrids.ImTime
        @assert tgrid.β ≈ mesh.β "Target grid has to have the same inverse temperature as the source grid."
        @assert tgrid.isFermi ≈ mesh.isFermi "Target grid has to have the same statistics as the source grid."
        # @assert tgrid.Euv ≈ mesh.Euv "Target grid has to have the same Euv as the source grid."
    elseif tgrid === nothing
        tgrid = MeshGrids.ImTime(mesh.β, mesh.isFermi; grid=mesh.dlr.τ, Euv=mesh.Euv)
    else
        tgrid = MeshGrids.ImTime(mesh.β, mesh.isFermi; grid=tgrid, Euv=mesh.Euv)
    end

    mesh_new = _replace_mesh(obj.mesh, mesh, tgrid)
    # mesh_new = (obj.mesh[1:dim-1]..., tgrid, obj.mesh[dim+1:end]...)
    data = dlr2tau(mesh.dlr, obj.data, tgrid.grid; axis=dim)
    return MeshArray(mesh=mesh_new, dtype=eltype(data), data=data)
end

@generated function replace_mesh_tuple(mesh::MT, N::Int, dim::Int, mesh_new::M) where {MT,M}
    if dim == 1
        m = :(mesh_new)
        for i in 2:N
            m = :($m, mesh[$i])
        end
    else
        m = :(mesh[1])
        for i in 2:N
            m = i == dim ? :($m, mesh_new) : :($m, mesh[$i])
        end
    end
    return :($m)
end

@generated function _replace_mesh(meshes::MT, mesh_old::OldM, mesh_new::NewM) where {MT,OldM,NewM}
    # return a new mesh tuple where mesh_old is replaced by mesh_new, if not found, the original meshes will be returned

    types = fieldtypes(MT)
    if types[1] == OldM
        m = :(mesh_new)
    else
        m = :(meshes[1])
    end

    for (i, t) in enumerate(types)
        # Core.println(t, ", ", M)
        if i <= 1
            continue
        else
            if t == OldM
                m = :($m, mesh_new)
            else
                m = :($m, meshes[$i])
            end
        end
    end
    # return :($m)

    # the following will always return a tuple
    if length(types) == 1
        return :($m,)
    else
        return :($m)
    end
end

#TODO: we need a version with the type of mesh as the second argument
@generated function _find_mesh(meshs::MT, mesh::M) where {MT,M}
    #find the first M type mesh in obj.mesh, if not found, return 0
    # Core.println(MT, ", ", M)
    for (i, t) in enumerate(fieldtypes(MT))

        # type equality is implemented as t<:M and M<:t, 
        # see https://discourse.julialang.org/t/how-to-test-type-equality/14144/6?u=mrbug
        if t == M
            return :($i)
        end
    end
    return :(0)
end


"""
    function dlr_to_imfreq(obj::MeshArray, ngrid=nothing; dim=nothing)

Transform a Green's function in DLR space to Matsubara frequency space. 
#Arguements
- 'obj': Function in DLR space
- `ngrid`: The Matsubara-frequency grid which the function transforms into. Default value is the Matsubara-frequency grid from the `DLRFreq` constructor.
- `dim`: The dimension of the temporal mesh. Default value is the first ImFreq mesh.
"""
# function dlr_to_imfreq(obj::MeshArray{T,N,MT}, ngrid=nothing; dim::Union{Nothing,Int}=nothing) where {T,N,MT}
function dlr_to_imfreq(obj::MeshArray{T,N,D,MT},
    ngrid=nothing;
    dim::Int=-1) where {T,N,MT,D}
    ########################## generic interface #################################
    if dim <= 0
        ind = findfirst(x -> (x isa MeshGrids.DLRFreq), obj.mesh)
        dim = isnothing(ind) ? -1 : ind
    end
    # dim = _find_mesh(obj.mesh, Type{MeshGrids.DLRFreq})
    # println(dim)
    @assert dim > 0 "No temporal can be transformed to imfreq."

    mesh = obj.mesh[dim]::MeshGrids.DLRFreq
    @assert mesh isa MeshGrids.DLRFreq "DLRFreq is expect for the dim = $dim."

    if ngrid isa MeshGrids.ImFreq
        @assert ngrid.β ≈ mesh.β "Target grid has to have the same inverse temperature as the source grid."
        @assert ngrid.isFermi ≈ mesh.isFermi "Target grid has to have the same statistics as the source grid."
        # @assert ngrid.Euv ≈ mesh.Euv "Target grid has to have the same Euv as the source grid."
    elseif isnothing(ngrid)
        ngrid = MeshGrids.ImFreq(mesh.β, mesh.isFermi; grid=mesh.dlr.n, Euv=mesh.Euv)
    else
        ngrid = MeshGrids.ImFreq(mesh.β, mesh.isFermi; grid=ngrid, Euv=mesh.Euv)
    end

    mesh_new = _replace_mesh(obj.mesh, mesh, ngrid)
    # mesh_new = (obj.mesh[1:dim-1]..., ngrid, obj.mesh[dim+1:end]...)
    data = dlr2matfreq(mesh.dlr, obj.data, ngrid.grid.grid; axis=dim)
    return MeshArray(mesh=mesh_new, dtype=Base.eltype(data), data=data)
end

"""
    function imfreq_to_dlr(obj::MeshArray; dim=nothing, rtol=1e-12, sym=:none)

Calculate the DLR sepctral density of a Matsubara-frequency Green's function.
#Arguements
- 'obj': Function in the Matsubara-frequency space.
- `dlrgrid`: The DLR grid which the function transforms into. Default value is the DLR grid constructed from the `ImFreq` members and the optional arguments.
- `dim`: The dimension of the mesh to be transformed. Default value is the first dimension with mesh type ImFreq.
- `rtol`: The relative tolerance of the DLR transform. Default value is 1e-12.
- `sym`: The symmetry of the Green's function, :none, :ph or :pha. Default value is :none.
"""
function imfreq_to_dlr(obj::MeshArray{T,N,D,MT}, dlrgrid::Union{Nothing,DLRFreq}=nothing; dim::Union{Nothing,Int}=nothing, rtol=1e-12, sym=:none) where {T,N,D,MT}
    if isnothing(dim)
        dim = findfirst(x -> (x isa MeshGrids.ImFreq), obj.mesh)
        @assert isnothing(dim) == false "No temporal can be transformed to dlr."
    end

    mesh = obj.mesh[dim]
    @assert mesh isa MeshGrids.ImFreq "ImFreq is expect for the dim = $dim."
    if dlrgrid === nothing
        dlrgrid = MeshGrids.DLRFreq(mesh.β, mesh.isFermi; Euv=mesh.Euv, rtol=rtol, sym=sym)
    else
        @assert dlrgrid.β ≈ mesh.β "Target grid has to have the same inverse temperature as the source grid."
        @assert dlrgrid.isFermi ≈ mesh.isFermi "Target grid has to have the same statistics as the source grid."
        # @assert dlrgrid.Euv ≈ mesh.Euv "Target grid has to have the same Euv as the source grid."
    end

    mesh_new = _replace_mesh(obj.mesh, mesh, dlrgrid)
    # mesh_new = (obj.mesh[1:dim-1]..., dlrgrid, obj.mesh[dim+1:end]...)
    data = matfreq2dlr(dlrgrid.dlr, obj.data, mesh.grid.grid; axis=dim) # should be mesh.grid.grid here
    return MeshArray(mesh=mesh_new, dtype=eltype(data), data=data)
end

"""
    function imtime_to_dlr(obj::MeshArray; dim=nothing, rtol=1e-12, sym=:none)

Calculate the DLR sepctral density of an imaginary-time Green's function.
#Arguements
- 'obj': Function in the imaginary-time space.
- `dlrgrid`: The DLR grid which the function transforms into. Default value is the DLR grid constructed from the `ImTime` members and the optional arguments.
- `dim`: The dimension of the mesh to be transformed. Default value is the first dimension with mesh type ImTime.
- `rtol`: The relative tolerance of the DLR transform. Default value is 1e-12.
- `sym`: The symmetry of the Green's function, :none, :ph or :pha. Default value is :none.
"""
function imtime_to_dlr(obj::MeshArray{T,N,D,MT}, dlrgrid::Union{Nothing,DLRFreq}=nothing; dim::Union{Nothing,Int}=nothing, rtol=1e-12, sym=:none) where {T,N,D,MT}
    if isnothing(dim)
        dim = findfirst(x -> (x isa MeshGrids.ImTime), obj.mesh)
        @assert isnothing(dim) == false "No temporal can be transformed to imtime."
    end

    mesh = obj.mesh[dim]
    @assert mesh isa MeshGrids.ImTime "ImTime is expect for the dim = $dim."
    if dlrgrid === nothing
        dlrgrid = MeshGrids.DLRFreq(mesh.β, mesh.isFermi; Euv=mesh.Euv, rtol=rtol, sym=sym)
    else
        @assert dlrgrid.β ≈ mesh.β "Target grid has to have the same inverse temperature as the source grid."
        @assert dlrgrid.isFermi ≈ mesh.isFermi "Target grid has to have the same statistics as the source grid."
        # @assert dlrgrid.Euv ≈ mesh.Euv "Target grid has to have the same Euv as the source grid."
    end

    mesh_new = _replace_mesh(obj.mesh, mesh, dlrgrid)
    # mesh_new = (obj.mesh[1:dim-1]..., dlrgrid, obj.mesh[dim+1:end]...)

    # println(mesh_new)
    data = tau2dlr(dlrgrid.dlr, obj.data, mesh.grid.grid; axis=dim)
    return MeshArray(mesh=mesh_new, dtype=eltype(data), data=data)
end

"""
    function to_dlr(obj::MeshArray; dim=nothing, rtol=1e-12, sym=:none)

Calculate the DLR sepctral density of an imaginary-time or Matsubara-frequency Green's function.
#Arguements
- 'obj': Function in the imaginary-time space or in the Matsubara-frequency
- `dim`: The dimension of the mesh to be transformed. Default value is the first dimension with mesh type ImTime or ImFreq.
- `rtol`: The relative tolerance of the DLR transform. Default value is 1e-12.
- `sym`: The symmetry of the Green's function, :none, :ph or :pha. Default value is :none.
"""
function to_dlr(obj::MeshArray{T,N,D,MT}, dlrgrid::Union{Nothing,DLRFreq}=nothing; dim::Union{Nothing,Int}=nothing, rtol=1e-12, sym=:none) where {T,N,D,MT}
    if isnothing(dim)
        dim = findfirst(x -> (x isa MeshGrids.ImTime || x isa MeshGrids.ImFreq), obj.mesh)
        @assert isnothing(dim) == false "No temporal can be transformed to imtime."
    end

    if obj.mesh[dim] isa MeshGrids.ImTime
        return imtime_to_dlr(obj, dlrgrid; dim=dim, rtol=rtol, sym=sym)
    elseif obj.mesh[dim] isa MeshGrids.ImFreq
        return imfreq_to_dlr(obj, dlrgrid; dim=dim, rtol=rtol, sym=sym)
    else
        error("ImTime or ImFreq is expect for the dim = $dim.")
    end
end

# function Base.:<<(Obj::MeshArray, objSrc::Expr)
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
# function density(obj::MeshArray; kwargs...)
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
# function from_L_G_R(self,L,G::MeshArray,R)
#     return 1
# end