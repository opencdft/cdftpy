#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Main RISM module
"""
import os
import sys
import time
from types import SimpleNamespace

import numpy as np

from cdftpy import __version__
from cdftpy.cdft1d.coulomb import compute_long_range_coul_pot_kspace, compute_coulomb_potential
from cdftpy.cdft1d.coulomb import compute_long_range_coul_pot_rspace
from cdftpy.cdft1d.coulomb import compute_short_range_coul_pot_rspace
from cdftpy.cdft1d.diis import diis_session
from cdftpy.cdft1d.exceptions import ConvergenceError
from cdftpy.cdft1d.io_utils import print_banner, print_solute, print_parameters
from cdftpy.cdft1d.io_utils import print_simulation
from cdftpy.cdft1d.potential import compute_lj_potential
from cdftpy.cdft1d.solvent import solvent_model_locate, Solvent
from cdftpy.utils.units import R

HEADER = F"""
==================================
1D RISM PROGRAM

version {__version__}
Marat Valiev and Gennady Chuev
==================================
"""

kb = R
PI = np.pi

DEFAULT_PARAMS = dict(diis_iterations=2, tol=1.0e-9, output_rate=10, max_iter=200, rcoul=1.25)

# DG_STRING = "|\u0394\u03B3|"
DG_STRING = "d_g"

DEFAULT_PRINT_LEVEL=frozenset(['header','parameters','solvent','solute'])


def rism_1d(solute, solvent, params=None, quiet=False, gr_guess=None,
            print_level = None):

    if print_level is None:
        print_level = DEFAULT_PRINT_LEVEL

    if quiet:
        sys.stdout = open(os.devnull, 'w')

    qs = solute["charge"]
    sig_s = solute["sigma"]
    eps_s = solute["eps"]

    rho_0 = solvent.density
    sig_v = solvent.sigma
    eps_v = solvent.eps
    qv = solvent.charge

    if params is None:
        params = DEFAULT_PARAMS
    params = {**DEFAULT_PARAMS, **params}

    ndiis = int(params["diis_iterations"])
    rcoul = float(params["rcoul"])
    tol = float(params["tol"])
    max_iter = int(params["max_iter"])
    output_rate = int(params["output_rate"])

    if "temp" not in params:
        params["temp"] = solvent.temp

    temp = float(params["temp"])
    beta = 1.0 / (kb * temp)

    if "rmax" in params:
        rmax = params["rmax"]
        solvent.extend(rmax)
    else:
        rmax = solvent.ifft.rgrid[-1]

    if 'header' in print_level:
        print(HEADER)
    if 'solute' in print_level:
        print_solute(solute)
    if 'parameters' in print_level:
        print_parameters(params)
    if 'solvent' in print_level:
        solvent.report()

    # initialize fft
    ifft = solvent.ifft
    rgrid = ifft.rgrid
    kgrid = ifft.kgrid

    # calculate lj potential
    vlj_r = beta * compute_lj_potential(sig_s, eps_s, sig_v, eps_v, rgrid)

    # coulomb long and short range
    vl_k = beta * compute_long_range_coul_pot_kspace(qs, qv, kgrid, r_s=rcoul)
    vl_r = beta * compute_long_range_coul_pot_rspace(qs, qv, rgrid, r_s=rcoul)
    vcs_r = beta * compute_short_range_coul_pot_rspace(qs, qv, rgrid, r_s=rcoul)

    # total short range potential
    vs_r = vcs_r + vlj_r

    # structure factor
    s_k = solvent.s_k
    delta = np.diagflat(np.ones(s_k.shape[0]))
    hbar = s_k - np.expand_dims(delta, axis=2)

    # compute long range part of h as -beta S*v_l
    gl_k = -np.einsum("abn,bn->an", s_k, vl_k)
    gl_r = np.apply_along_axis(ifft.to_rspace, 1, gl_k)

    if gr_guess is None:
        g_r = np.zeros(hbar[0].shape)
    else:
        g_r = gr_guess

    print("")
    print_banner("   Self-consistent cycle     ")

    converged = False

    diis_update = diis_session()

    for it in range(max_iter):

        # calculate h
        h_r = np.exp(-vs_r + g_r) - 1.0

        # calculate short range c
        cs_r = h_r - g_r
        cs_k = np.apply_along_axis(ifft.to_kspace, 1, cs_r)

        # calculate gamma
        hbar_cs_k = np.einsum("abn,bn->an", hbar, cs_k)
        gn_r = np.apply_along_axis(ifft.to_rspace, 1, hbar_cs_k) + gl_r

        # compute error
        dg_r = gn_r - g_r
        err = np.sum(dg_r ** 2)
        err = np.sqrt(err / gn_r.size)

        converged = err < tol

        if it % output_rate == 0 or converged:
            fe, _ = compute_free_energy(beta, rho_0, ifft, vl_r, g_r, h_r, cs_r)
            if it == 0:
                print(f"{'iter':<5} {DG_STRING:<12}{'Free Energy':<10} ")
            print(f"{it:<5} {err:>.2e}   {fe:<.7f}")

        if converged:
            print(f"\nReached convergence, {DG_STRING} < {tol}")
            break

        g_r = diis_update(ndiis, gn_r, dg_r)

    if not converged:
        print(
            f"Could not reach specified convergence criteria after {max_iter} iterations"
        )
        raise ConvergenceError

    print("\n")
    print(f"{'Total Free Energy ':<30} {fe:>12.6f}")

    fe, fe_extra = compute_free_energy(beta, rho_0, ifft, vl_r, g_r, h_r, cs_r)

    fe_gas, fe_inter = fe_extra
    h_k = np.apply_along_axis(ifft.to_kspace, 1, h_r)
    epot_r, epot_k = compute_coulomb_potential(rho_0, qv, ifft, h_k)

    sys.stdout = sys.__stdout__

    return SimpleNamespace(
        method=F"rism v{__version__}",
        solute=SimpleNamespace(**solute),
        solvent=solvent,
        ifft=ifft,
        gl_r=gl_r,
        gl_k=gl_k,
        vl_r=vl_r,
        vl_k=vl_k,
        vcs_r=vcs_r,
        vs_r=vs_r,
        h_r=h_r,
        h_k=h_k,
        epot_r=epot_r,
        epot_k=epot_k,
        g_r=g_r,
        fe_tot=fe,
        fe_gas=fe_gas,
        fe_inter=fe_inter,
        phi_r=-(g_r + vl_r) / beta,
        beta=beta,
        rmax=rmax
    )


def compute_free_energy(beta, rho_0, ifft, vl_r, g_r, h_r, c_r):
    phi_r = -(g_r + vl_r) / beta

    fe0 = -rho_0 * ifft.integrate_rspace(c_r) / beta
    fe1 = 0.5 * ifft.integrate_rspace(rho_0 * phi_r * h_r)
    fe_tot = fe0 - fe1
    return fe_tot, (fe_tot - fe1, fe1)


if __name__ == "__main__":
    # load solvent model
    solvent_name = "s2"
    filename = solvent_model_locate(solvent_name)
    solvent = Solvent.from_file(filename, rism_patch=True)
    params = dict(diis_iterations=2, tol=1.0e-7, max_iter=600, rmax=500)

    solute = dict(name="Cl", charge=-1.0, sigma=4.83, eps=0.05349244)

    rmax = solvent.ifft.rgrid[-1]

    start_time = time.process_time()
    sim = rism_1d(solute, solvent, params=params)
    t = time.process_time() - start_time
