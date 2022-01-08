#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Legacy RISM module
"""
import os
import sys
from types import SimpleNamespace

import numpy as np

from cdftpy.cdft1d.coulomb import compute_coulomb_energy
from cdftpy.cdft1d.coulomb import compute_coulomb_potential
from cdftpy.cdft1d.coulomb import compute_long_range_coul_pot_kspace
from cdftpy.cdft1d.coulomb import compute_long_range_coul_pot_rspace
from cdftpy.cdft1d.coulomb import compute_short_range_coul_pot_rspace
from cdftpy.cdft1d.diis import diis_session
from cdftpy.cdft1d.io_utils import print_banner
from cdftpy.cdft1d.io_utils import print_simulation
from cdftpy.cdft1d.potential import compute_lj_potential
from cdftpy.cdft1d.rdf import analyze_rdf_peaks_sim
from cdftpy.cdft1d.rdf import write_rdf_sim
from cdftpy.cdft1d.solvent import solvent_model_locate, Solvent
from cdftpy.utils.units import R

from cdftpy.cdft1d._version import __version__

HEADER = F"""
==================================
1D RISM PROGRAM

version {__version__}
Marat Valiev and Gennady Chuev
==================================
"""

kb = R
PI = np.pi

DEFAULT_PARAMS = dict(diis_iterations=2, tol=1.0e-9, output_rate=10, max_iter=400, rcoul=1.25)


# DG_STRING = "|\u0394\u03B3|"
DG_STRING = "d_g"


def rism_1d(solute, solvent_model, params=None, verbose=True):

    if not verbose:
        sys.stdout = open(os.devnull, 'w')

    qs = solute["charge"]
    sig_s = solute["sigma"]
    eps_s = solute["eps"]

    rho_0 = solvent_model.density
    sig_v = solvent_model.sigma
    eps_v = solvent_model.eps
    qv = solvent_model.charge
    s_k = solvent_model.s_k

    if s_k is None:
        print("missing structure factor in solvent model")
        print(
            f"check if solvent file {solvent_model.filename} is appropriate for RISM calculations"
        )
        sys.exit(1)
    if params is None:
        params = DEFAULT_PARAMS
    params = {**DEFAULT_PARAMS, **params}

    ndiis = int(params["diis_iterations"])
    rcoul= float(params["rcoul"])
    tol = float(params["tol"])
    max_iter = int(params["max_iter"])
    output_rate = int(params["output_rate"])

    if "temp" in params:
        temp = float(params["temp"])
    else:
        temp = solvent_model.temp
    beta = 1.0 / (kb * temp)

    print(HEADER)
    print_simulation(solute, solvent_model, params)

    # initialize fft
    ifft = solvent_model.ifft
    rgrid = ifft.rgrid
    kgrid = ifft.kgrid

    # calculate lj potential
    vlj_r = beta * compute_lj_potential(sig_s, eps_s, sig_v, eps_v, rgrid)

    # coulomb long and short range
    vl_k  = beta * compute_long_range_coul_pot_kspace(qs, qv, kgrid, r_s=rcoul)
    vl_r  = beta * compute_long_range_coul_pot_rspace(qs, qv, rgrid, r_s=rcoul)
    vcs_r = beta * compute_short_range_coul_pot_rspace(qs, qv, rgrid, r_s=rcoul)

    # total short range potential
    vs_r = vcs_r + vlj_r

    # compute long range part of h as -beta S*v_l
    hl_k = -np.einsum("abn,bn->an", s_k, vl_k)
    hl_r = np.apply_along_axis(ifft.to_rspace, 1, hl_k)

    # initial guess for gamma
    g_r = np.zeros(s_k[0].shape)

    print("")
    print_banner("   Self-consistent cycle     ")

    converged = False
    diis_update = diis_session()

    for it in range(max_iter):

        c_r = np.exp(-vs_r + g_r) - g_r - 1.0
        c_k = np.apply_along_axis(ifft.to_kspace, 1, c_r)

        hs_k = np.einsum("abn,bn->an", s_k, c_k)
        hs_r = np.apply_along_axis(ifft.to_rspace, 1, hs_k)

        # new guess for gamma
        gn_r = hl_r + hs_r - c_r

        h_r = hl_r + hs_r
        h_k = hl_k + hs_k
        # compute error
        dg_r = gn_r - g_r
        err = np.sum(dg_r ** 2)
        err = np.sqrt(err / gn_r.size)

        if it % output_rate == 0 or converged:
            fe = compute_free_energy(
                beta, rho_0, ifft, qs, qv, c_k, h_r, h_k, vl_r, g_r
            )
            if it == 0:
                print(f"{'iter':<5} {DG_STRING:<11}{'Free Energy':<10} ")
            print(f"{it:<5} {err:>.2e}   {fe.total:<.7f}")

        converged = err < tol
        if converged:
            print(f"Reached specified convergence criteria of {DG_STRING} < {tol}")
            break

        g_r = diis_update(ndiis, gn_r, dg_r)

    if not converged:
        print(
            f"Could not reach specified convergence criteria after {max_iter} iterations"
        )
        return

    fe = compute_free_energy(beta, rho_0, ifft, qs, qv, c_k, h_r, h_k, vl_r, g_r)

    print("\n")
    print(f"{'Free Energy Total':<30} {fe.total:>12.6f}")
    print(f"{'':-<30}-{'':-<12}")
    print(f"{'Volume contribution':<30} {fe.volume:>12.6f}")
    print(f"{'Surface contribution':<30} {fe.surface:>12.6f}")
    print(f"{'Coulomb contribution':<30} {fe.coulomb:>12.6f}")

    epot_r, epot_k = compute_coulomb_potential(rho_0, qv, ifft, h_k)

    sys.stdout = sys.__stdout__
    return SimpleNamespace(
        solute=SimpleNamespace(**solute),
        solvent=solvent_model,
        ifft=ifft,
        vl_r=vl_r,
        vl_k=vl_k,
        vcs_r=vcs_r,
        h_r=h_r,
        h_k=h_k,
        hl_r=hl_r,
        hl_k=hl_k,
        g_r=g_r,
        fe_tot=fe.total,
        epot_r=epot_r,
        epot_k=epot_k,
    )


def compute_free_energy(beta, rho_0, ifft, qs, qv, c_k, h_r, h_k, vl_r, g_r):

    # volume contribution
    f_vol = -0.5 * np.sum(c_k[:, 0] + h_k[:, 0]) * rho_0 / beta

    # surface + coulomb
    gtot_r = g_r + vl_r
    rho_r = h_r + 1.0
    f_surf_coul = 0.5 * rho_0 * ifft.integrate_rspace(gtot_r * rho_r) / beta

    epot_r, epot_k = compute_coulomb_potential(rho_0, qv, ifft, h_k)

    f_coul = compute_coulomb_energy(qs, epot_r)
    f_surf = f_surf_coul - f_coul

    f_total = f_vol + f_surf_coul

    return SimpleNamespace(total=f_total, surface=f_surf, volume=f_vol, coulomb=f_coul)

