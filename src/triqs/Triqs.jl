module Triqs

import ..MeshArray
import ..MeshGrids
import ..BrillouinZoneMeshes

using PythonCall
# gf = pyimport("triqs.gf")
# sys = pyimport("sys")

function _get_statistics(mesh)
    statis = pyconvert(String, mesh.statistic)
    if (statis == "Fermion")
        stat = true
    else # Boson or other
        stat = false
    end
    return stat
end

function _get_mesh_from_triqs(triqs_mesh)
    gf = pyimport("triqs.gf")
    ############ load time/freq mesh #############
    if pyisinstance(triqs_mesh, gf.meshes.MeshImTime)
        β = pyconvert(Float64, triqs_mesh.beta)
        stat = _get_statistics(triqs_mesh)
        grid_t = pyconvert(Vector, triqs_mesh.values())
        tgrid = MeshGrids.ImTime(β, stat, grid=grid_t)
        return tgrid
    elseif pyisinstance(triqs_mesh, gf.meshes.MeshImFreq)
        β = pyconvert(Float64, triqs_mesh.beta)
        stat = _get_statistics(triqs_mesh)
        grid_first_index = pyconvert(Int32, triqs_mesh.first_index())
        grid_last_index = pyconvert(Int32, triqs_mesh.last_index())
        grid_t = [grid_first_index:grid_last_index;]
        tgrid = MeshGrids.ImFreq(β, stat, grid=grid_t)
        return tgrid
    elseif pyisinstance(triqs_mesh, gf.meshes.MeshBrZone)
        # when import from python, linear_index remain consistent
        # while cartesian index and order of lattice vector reversed
        mkdims = pyconvert(Array, triqs_mesh.dims)
        mkunits = pyconvert(Array, triqs_mesh.units)
        DIM = count(i -> (i != 1), mkdims) # actual dimension of the grid
        nk = mkdims[1]
        latvec = pyconvert(Array, mkunits)[1:DIM, 1:DIM]' .* nk
        latvec = reverse(latvec, dims=2)
        umesh = BrillouinZoneMeshes.BaseMesh.UniformMesh{DIM,nk,BrillouinZoneMeshes.BaseMesh.EdgedMesh}([0.0, 0.0], latvec)
        return umesh
    elseif pyisinstance(triqs_mesh, gf.mesh_product.MeshProduct)
        # expand MeshProduct to a tuple
        rank = pyconvert(Int, triqs_mesh.rank)
        mlist = Tuple(_get_mesh_from_triqs(triqs_mesh[i-1]) for i in rank:-1:1)
        return mlist
    else
        error("Unknown mesh type:", pytype(triqs_mesh))
    end
end

function MeshArray(objSrc::Py)

    innerstate = pyconvert(Tuple, objSrc.target_shape)
    # @assert innerstate == dims[1:length(innerstate)] "Inner state dimensions do not match!"
    mesh = (1:state for state in innerstate[end:-1:1])
    triqmesh = _get_mesh_from_triqs(objSrc.mesh)
    if isa(triqmesh, Tuple)
        mesh = (mesh..., triqmesh...)
    else
        mesh = (mesh..., triqmesh)
    end
    _data = PyArray(objSrc.data, copy=false) #no copy is made, but PyArray will be in column-major 
    g = MeshArray(mesh...; dtype=Float64)
    for i in 1:length(g)
        g.data[i] = unsafe_load(_data.ptr, i) #read data from pointer
    end
    return g
end

function Base.:<<(obj::MeshArray{T,N,MT}, objSrc::Py) where {T,MT,N}
    @assert obj.dims[end:-1:1] == pyconvert(Tuple, objSrc.data.shape) "Dimensions do not match!"
    _data = PyArray(objSrc.data, copy=false) #no copy is made, but PyArray will be in column-major 
    for i in 1:length(obj)
        obj.data[i] = unsafe_load(_data.ptr, i) #read data from pointer
    end
    return obj
end

# function Base.:<<(obj::MeshArray{T,N,MT}, objSrc::Py) where {T,MT,N}

#     innerstate = pyconvert(Tuple, objSrc.target_shape)
#     @assert innerstate == dims[1:length(innerstate)] "Inner state dimensions do not match!"

#     β = pyconvert(Float64, objSrc.mesh.beta)
#     statis = pyconvert(String, objSrc.mesh.statistic)
#     if (statis == "Fermion")
#         stat = MeshGrids.FERMI
#     elseif (statis == "Boson")
#         stat = MeshGrids.BOSE
#     else
#         stat = MeshGrids.UNKNOWN
#     end
#     if pyisinstance(objSrc.mesh, gf.meshes.MeshImTime)
#         grid_t = pyconvert(Vector, objSrc.mesh.values())
#         tgrid = GreenFunc.MeshGrids.ImTime(β, stat, grid=grid_t)
#     elseif pyisinstance(objSrc.mesh, gf.meshes.MeshImFreq)
#         grid_first_index = pyconvert(Int32, objSrc.mesh.first_index())
#         grid_last_index = pyconvert(Int32, objSrc.mesh.last_index())
#         grid_t = [grid_first_index:grid_last_index;]
#         tgrid = GreenFunc.MeshGrids.ImFreq(β, stat, grid=grid_t)
#     else
#         error("Unknown mesh type")
#     end
#     data_t = pyconvert(Array, objSrc.data)
#     obj.mesh = ()
#     return Obj = GreenFunc.MeshArray(tgrid, 1:tar_sh[1], 1:tar_sh[2], dtype=eltype(data_t), data=data_t)
# end

end
