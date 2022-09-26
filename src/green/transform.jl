"""
    Base.:<<(objL::ManifoldArray, objR::ManifoldArray)

DLR Fourier transform of functions that has exactly one TemporalGrid(ImTime, ImFreq or DLRFreq) among the meshes. 
If objL and objR have identical TemporalGrid, objL<<objR assign objR to objL.
If objL and objR have different TemporalGrid, one of them has to be in DLR space.
If objL is in DLR space, objL<<objR calculates the DLR spectral density of data in objR
if objR is in DLR space, objL<<objR calculates the corresponding data from the DLR spectral density in objR.
"""

function Base.:<<(objL::ManifoldArray, objR::ManifoldArray)
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
            objL.data = matfreq2dlr(meshL.dlr, objR.data, meshR.grid; axis=axes[1])
        elseif typeR <: MeshGrids.ImTime
            objL.data = tau2dlr(meshL.dlr, objR.data, meshR.grid; axis=axes[1])
        end
    elseif typeR <: MeshGrids.DLRFreq
        if typeL <: MeshGrids.ImFreq
            objL.data = dlr2matfreq(meshR.dlr, objR.data, meshL.grid; axis=axes[1])
        elseif typeL <: MeshGrids.ImTime
            objL.data = dlr2tau(meshR.dlr, objR.data, meshL.grid; axis=axes[1])
        end
    else
        error("One of the Grren's function has to be in DLRfreq space to do Fourier transform")
    end
end


"""
    function dlr_to_imtime(obj::ManifoldArray; kwargs...)

Convert sepctral density in DLR space to imaginary time space for function with exactly one DLRFreq grid among meshes.
#Arguements
#- 'obj': Function in DLR space
#- 'targetGrid': The imaginary time grid which the function transforms into. Default value is the imaginary time space from the DLR grid in obj.
"""
function dlr_to_imtime(obj::ManifoldArray; kwargs...)
    # init version of <<
    # more general version needed
    axes = []
    for (mi, mesh) in enumerate(obj.mesh)
        if (typeof(mesh) <: TemporalGrid)
            @assert typeof(mesh) <: MeshGrids.DLRFreq "Green's function has to be in DLRFreq space"
            append!(axes, mi)
        end
    end
    @assert length(axes) == 1 "Only one temporal mesh with built in DLR grid is allowed"
    mesh = obj.mesh[axes[1]]
    if :targetGrid in keys(kwargs)
        tgrid = kwargs[:targetGrid]
        @assert typeof(tgrid) == MeshGrids.ImTime "Target grid has to be ImTime type"
    else
        tgrid = MeshGrids.ImTime(mesh.β, mesh.statistics)
    end

    mesh_copy = Tuple(i == axes[1] ? tgrid : obj.mesh[i] for i in 1:length(obj.dims))
    data = dlr2tau(mesh.dlr, obj.data, tgrid.grid; axis=axes[1])
    datatype = typeof(data[1])
    return ManifoldArray(mesh_copy...; dtype=datatype, data=data)
end

"""
    function dlr_to_imfreq(obj::ManifoldArray; kwargs...)

Convert sepctral density in DLR space to matsubara frequency space for function with exactly one DLRFreq grid among meshes.
#Arguements
#- 'obj': Function in DLR space
#- 'targetGrid': The matsubara frequency grid which the function transforms into. Default value is the imaginary time space from the DLR grid in obj.
"""

function dlr_to_imfreq(obj::ManifoldArray; kwargs...)
    # init version of <<
    # more general version needed
    axes = []
    for (mi, mesh) in enumerate(obj.mesh)
        if (typeof(mesh) <: TemporalGrid)
            @assert typeof(mesh) <: MeshGrids.DLRFreq "Green's function has to be in DLRFreq space"
            append!(axes, mi)
        end
    end
    @assert length(axes) == 1 "Only one temporal mesh with built in DLR grid is allowed"
    mesh = obj.mesh[axes[1]]
    if :targetGrid in keys(kwargs)
        tgrid = kwargs[:targetGrid]
        @assert typeof(tgrid) == MeshGrids.ImFreq "Target grid has to be MeshGrids.ImFreq type"
    else
        tgrid = MeshGrids.ImFreq(mesh.β, mesh.statistics)
    end
    mesh_copy = Tuple(i == axes[1] ? tgrid : obj.mesh[i] for i in 1:length(obj.dims))
    data = dlr2matfreq(mesh.dlr, obj.data, tgrid.grid; axis=axes[1])
    datatype = typeof(data[1])
    return ManifoldArray(mesh_copy...; dtype=datatype, data=data)
end

"""
    function to_dlr(obj::ManifoldArray; kwargs...)

Calculate the sepctral density of a function in imaginary time/matsubara frequency space, with exactly one DLRFreq grid among meshes.
#Arguements
#- 'obj': Function in imaginary time/matsubara frequency space
#- 'targetGrid': The DLRFreq grid which the function transforms into. Default value is a DLRFreq gridwith the same beta and statistics of the temporal grid in obj, with default rtol=1e-12 and Euv = 1000/beta. 
"""

function to_dlr(obj::ManifoldArray; kwargs...)
    # init version of <<
    # more general version needed
    axes = []
    for (mi, mesh) in enumerate(obj.mesh)
        if (typeof(mesh) <: TemporalGrid)
            @assert typeof(mesh) <: MeshGrids.ImTime || typeof(mesh) <: MeshGrids.ImFreq "Green's function has to be in ImFreq or ImTime space"
            append!(axes, mi)
        end
    end
    @assert length(axes) == 1 "Only one temporal mesh with built in DLR grid is allowed"
    mesh = obj.mesh[axes[1]]
    if :targetGrid in keys(kwargs)
        tgrid = kwargs[:targetGrid]
        @assert typeof(tgrid) == MeshGrids.DLRFreq "Target grid has to be DLRFreq type"
    else
        tgrid = MeshGrids.DLRFreq(mesh.β, mesh.statistics)
    end
    mesh_copy = Tuple(i == axes[1] ? tgrid : obj.mesh[i] for i in 1:length(obj.dims))
    if typeof(mesh) <: MeshGrids.ImFreq
        data = matfreq2dlr(tgrid.dlr, obj.data, mesh.grid; axis=axes[1])
    else
        typeof(mesh) <: MeshGrids.ImTime
        data = tau2dlr(tgrid.dlr, obj.data, mesh.grid; axis=axes[1])
    end
    datatype = typeof(data[1])
    return ManifoldArray(mesh_copy...; dtype=datatype, data=data)
end

# """
#     <<(Obj::GreenDLR, objSrc::Py)
#     Converts the green function from triqs to ManifoldArray.
# """

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
# function density(obj::ManifoldArray; kwargs...)
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
# function from_L_G_R(self,L,G::ManifoldArray,R)
#     return 1
# end