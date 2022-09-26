module Triqs

import ..ManifoldArray
import ..MeshGrids

using PythonCall
# gf = pyimport("triqs.gf")
# sys = pyimport("sys")

function _get_statistics(mesh)
    statis = pyconvert(String, mesh.statistic)
    if (statis == "Fermion")
        stat = MeshGrids.FERMI
    elseif (statis == "Boson")
        stat = MeshGrids.BOSE
    else
        stat = MeshGrids.UNKNOWN
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
    else
        error("Unknown mesh type")
    end
end

function ManifoldArray(objSrc::Py)

    innerstate = pyconvert(Tuple, objSrc.target_shape)
    # @assert innerstate == dims[1:length(innerstate)] "Inner state dimensions do not match!"
    mesh = (1:state for state in innerstate[end:-1:1])
    mesh = (mesh..., _get_mesh_from_triqs(objSrc.mesh))
    _data = PyArray(objSrc.data, copy=false) #no copy is made, but PyArray will be in column-major 
    g = ManifoldArray(mesh...; dtype=Float64)
    for i in 1:length(g)
        g.data[i] = unsafe_load(_data.ptr, i) #read data from pointer
    end
    return g
end

function Base.:<<(obj::ManifoldArray{T,N,MT}, objSrc::Py) where {T,MT,N}
    @assert obj.dims[end:-1:1] == pyconvert(Tuple, objSrc.data.shape) "Dimensions do not match!"
    _data = PyArray(objSrc.data, copy=false) #no copy is made, but PyArray will be in column-major 
    for i in 1:length(obj)
        obj.data[i] = unsafe_load(_data.ptr, i) #read data from pointer
    end
    return obj
end

# function Base.:<<(obj::ManifoldArray{T,N,MT}, objSrc::Py) where {T,MT,N}

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
#     return Obj = GreenFunc.ManifoldArray(tgrid, 1:tar_sh[1], 1:tar_sh[2], dtype=eltype(data_t), data=data_t)
# end

end
