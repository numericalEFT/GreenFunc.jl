using PythonCall
using GreenFunc

gf = pyimport("triqs.gf")
lat = pyimport("triqs.lattice")
tb = pyimport("triqs.lattice.tight_binding")
np = pyimport("numpy")

BL = lat.BravaisLattice(units=((1, 0, 0), (0, 1, 0))) #square lattice
# the following also works
# BL = lat.BravaisLattice(units=pylist([(1, 0, 0), (0, 1, 0)])) #square lattice
# but this will not work (you can not directly pass a julia vector as a python-list argument
# BL = lat.BravaisLattice(units=[(1, 0, 0), (0, 1, 0)]) #square lattice

nk = 20
mk = gf.MeshBrillouinZone(lat.BrillouinZone(BL), nk)
miw = gf.MeshImFreq(beta=1.0, S="Fermion", n_max=100) #grid number : 200
mprod = gf.MeshProduct(mk, miw)

G_w = gf.GfImFreq(mesh=miw, target_shape=[1, 1]) #G_w.data.shape will be [200, 1, 1]
G_k_w = gf.GfImFreq(mesh=mprod, target_shape=[1, 1]) #G_k_w.data.shape will be [400, 200, 1, 1]

t = 1.0

####### fill the Green's function with data ################
## only make small number of allocations: 0.534508 seconds (28.01 k allocations: 522.039 KiB)
for (ik, k) in enumerate(G_k_w.mesh[0])
    G_w << gf.inverse(gf.iOmega_n - 2 * t * (np.cos(k[0]) + np.cos(k[1])))
    G_k_w.data[ik-1, pyslice(0, nk^2), 0, 0] = G_w.data[pyslice(0, nk^2), 0, 0]
end

Gkw = from_triqs(G_k_w)
println(size(Gkw)) #should be (1, 1, 200, 400)

