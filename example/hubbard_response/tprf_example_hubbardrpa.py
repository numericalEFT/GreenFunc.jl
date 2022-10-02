from triqs_tprf.tight_binding import create_square_lattice

square_lattice = create_square_lattice(norb=1, t=1.0)

nk = 32
dim = 2

kmesh = square_lattice.get_kmesh(n_k=[nk] * dim + [1] * (3 - dim))
e_k = square_lattice.fourier(kmesh)

from triqs.gf import MeshImFreq
from triqs_tprf.lattice import lattice_dyson_g0_wk

wmesh = MeshImFreq(beta=10, S='Fermion', n_max=100)
g0_wk = lattice_dyson_g0_wk(mu=0, e_k=e_k, mesh=wmesh)

from triqs_tprf.lattice_utils import imtime_bubble_chi0_wk

chi0_wk = imtime_bubble_chi0_wk(g0_wk, nw=100)

import numpy as np
from triqs_tprf.lattice import solve_rpa_PH

U = 1.0 * np.ones(shape=(1, 1, 1, 1), dtype=complex)

chi_d_wk = solve_rpa_PH(chi0_wk, -U) # Minus here for correct density RPA equation
chi_m_wk = solve_rpa_PH(chi0_wk, U)

from triqs_tprf.lattice import construct_phi_wk

phi_d_wk = construct_phi_wk(chi_d_wk, U)
phi_m_wk = construct_phi_wk(chi_m_wk, U)

from triqs_tprf.eliashberg import construct_gamma_singlet_rpa

gamma_singlet = construct_gamma_singlet_rpa(U, U, phi_d_wk, phi_m_wk)

import functools
from triqs_tprf.symmetries import enforce_symmetry

variables = ["frequency", "momentum"]
symmetries = ["even", "even"]

symmetrize_freq_even_mom_even = functools.partial(enforce_symmetry, variables=variables, symmetries=symmetries)

symmetries = ["odd", "odd"]

symmetrize_freq_odd_mom_odd = functools.partial(enforce_symmetry, variables=variables, symmetries=symmetries)

from triqs_tprf.eliashberg import solve_eliashberg

lambdas_freq_even_mom_even, deltas_freq_even_mom_even = solve_eliashberg(gamma_singlet, g0_wk, symmetrize_fct=symmetrize_freq_even_mom_even, k=1)
lambdas_freq_odd_mom_odd, deltas_freq_odd_mom_odd = solve_eliashberg(gamma_singlet, g0_wk, symmetrize_fct=symmetrize_freq_odd_mom_odd, k=1)
