from triqs_tprf.tight_binding import create_square_lattice
from triqs.gf import MeshImFreq
from triqs_tprf.lattice import lattice_dyson_g0_wk
from triqs_tprf.lattice_utils import imtime_bubble_chi0_wk
import numpy as np
from triqs_tprf.lattice import solve_rpa_PH
from triqs_tprf.lattice import construct_phi_wk
from triqs_tprf.eliashberg import construct_gamma_singlet_rpa

def gamma_rpa(norb, t, nk, dim, beta, n_max, mu, U):
    square_lattice = create_square_lattice(norb=norb, t=t)
    kmesh = square_lattice.get_kmesh(n_k=[nk] * dim + [1] * (3 - dim))
    e_k = square_lattice.fourier(kmesh)
    wmesh = MeshImFreq(beta=beta, S='Fermion', n_max=n_max)
    g0_wk = lattice_dyson_g0_wk(mu=mu, e_k=e_k, mesh=wmesh)
    chi0_wk = imtime_bubble_chi0_wk(g0_wk, nw=n_max)
    UU = U * np.ones(shape=(1, 1, 1, 1), dtype=complex)
    chi_d_wk = solve_rpa_PH(chi0_wk, -UU) # Minus here for correct density RPA equation
    chi_m_wk = solve_rpa_PH(chi0_wk, UU)
    phi_d_wk = construct_phi_wk(chi_d_wk, UU)
    phi_m_wk = construct_phi_wk(chi_m_wk, UU)
    gamma_singlet = construct_gamma_singlet_rpa(UU, UU, phi_d_wk, phi_m_wk)
    return gamma_singlet, g0_wk

# gamma_rpa(norb=1,t=1.0,nk=32,dim=2,beta=100,n_max=1000,mu=0,U=1.0)