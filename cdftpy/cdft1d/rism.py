#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Main RISM module
"""
import sys


import numpy as np

from cdftpy.cdft1d._version import __version__

from cdftpy.cdft1d.diis import diis_session
from cdftpy.cdft1d.exceptions import ConvergenceError
from cdftpy.cdft1d.io_utils import get_banner

from cdftpy.cdft1d.loggers import get_stream_logger

HEADER = F"""
==================================
1D RISM PROGRAM

version {__version__}
Marat Valiev and Gennady Chuev
==================================
"""

PI = np.pi

DEFAULT_PARAMS = dict(diis_iterations=2, tol=1.0e-9, output_rate=10, max_iter=200, rcoul=1.25)

# DG_STRING = "|\u0394\u03B3|"
DG_STRING = "d_g"

def rism_1d(sim, quiet=True, capture=True):

    logger, streamer = get_stream_logger(__name__,
                                         capture=capture,
                                         debug=True,
                                         quiet=quiet)

    rho_0 = sim.rho_0

    ndiis = sim.ndiis
    tol = sim.tol
    max_iter = sim.max_iter
    output_rate = sim.output_rate

    beta = sim.beta

    # initialize fft
    ifft = sim.ifft

    vl_k = sim.vl_k
    vl_r = sim.vl_r
    vs_r = sim.vs_r


    # structure factor
    s_k = sim.s_k
    hbar = sim.hbar

    # compute long range part of h as -beta S*v_l
    gl_k = -np.einsum("abn,bn->an", s_k, vl_k)
    gl_r = np.apply_along_axis(ifft.to_rspace, 1, gl_k)

    # g_r = np.zeros(hbar[0].shape)
    g_r = np.array(gl_r)

    logger.info(get_banner("   Self-consistent cycle     "))

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
                logger.info(f"{'iter':<5} {DG_STRING:<12}{'Free Energy':<10} ")
            logger.info(f"{it:<5} {err:>.2e}   {fe:<.7f}")

        if converged:
            logger.info(f"\nReached convergence, {DG_STRING} < {tol}")
            break

        g_r = diis_update(ndiis, gn_r, dg_r)

    if not converged:
        logger.info(
            f"Could not reach specified convergence criteria after {max_iter} iterations"
        )
        if capture:
            print(streamer.getvalue())
        raise ConvergenceError(F"RISM calculation.\n Try increasing ndiis")

    logger.info("\n")
    logger.info(f"{'Total Free Energy ':<30} {fe:>12.6f} kj/mol")

    fe, fe_extra = compute_free_energy(beta, rho_0, ifft, vl_r, g_r, h_r, cs_r)

    sim.fe_tot = fe
    sim.h_r = h_r
    sim.g_r = g_r

    if streamer is not None:
        sim.log = streamer.getvalue()
    else:
        sim.log = None

    return fe

def compute_free_energy(beta, rho_0, ifft, vl_r, g_r, h_r, c_r):
    phi_r = -(g_r + vl_r) / beta

    fe0 = -rho_0 * ifft.integrate_rspace(c_r) / beta
    fe1 = 0.5 * ifft.integrate_rspace(rho_0 * phi_r * h_r)
    fe_tot = fe0 - fe1
    return fe_tot, (fe_tot - fe1, fe1)


