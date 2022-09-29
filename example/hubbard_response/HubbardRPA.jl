using GreenFunc
using GreenFunc.Triqs.PythonCall

#TODO: load rpa kernel of hubbard model from triqs and triqs_tprf

# load python script(import from file not allowed in PythonCall)
f = open("tprf_rpa.py", "r")
# exec python script
pyexec(read(f, String), Main)
gamma_rpa = pyeval(Py, "gamma_rpa", Main)

function Gamma(; norb=1, t=1.0, nk=32, dim=2, beta=10, n_max=100, mu=0, U=1.0)
    gamma = gamma_rpa(norb, t, nk, dim, beta, n_max, mu, U)
    return GreenFunc.MeshArray(gamma)
end

