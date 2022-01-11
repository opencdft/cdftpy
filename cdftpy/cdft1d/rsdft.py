#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Main RSDFT module
"""
import logging
import time

from cdftpy.cdft1d.exceptions import ConvergenceError
from cdftpy.cdft1d.loggers import get_stream_logger

logging.getLogger('matplotlib').setLevel(logging.WARNING)

import numpy as np

from cdftpy.cdft1d.diis import diis_session
from cdftpy.cdft1d.io_utils import get_banner

from cdftpy.cdft1d._version import __version__


PI = np.pi

# DG_STRING = "|\u0394\u03B3|"
DG_STRING = "d_g"

HEADER = F"""
==================================
1D RSDFT PROGRAM

version {__version__}
Marat Valiev and Gennady Chuev
==================================
"""


def rsdft_1d(sim, quiet=True, capture=True, bridge=False):

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

    ifft = sim.ifft

    hbar = sim.hbar
    sm = sim.sm

    vl_k = sim.vl_k
    vl_r = sim.vl_r
    vs_r = sim.vs_r
    zeta = sim.zeta

    gl_r, gl_k = compute_gamma_long_range(ifft, hbar, sm, vl_k, vl_r, zeta)

    # initial guess for gamma
    # g_r = np.zeros(hbar[0].shape)
    g_r = np.array(gl_r)


    logger.info(get_banner("   Self-consistent cycle     "))

    converged = False

    diis_update = diis_session()

    for it in range(max_iter):

        f_r, f_k = compute_mayer_function(ifft, vs_r, gl_r, gl_k, g_r)

        delta_h_r = compute_delta_h(ifft, sm, f_r, f_k)

        hm_r = compute_h_mol(f_r)

        h_r = hm_r + delta_h_r

        cs_k = compute_c_k(ifft, sm, h_r, g_r)

        gn_r = compute_g_r(ifft, hbar, cs_k, gl_r)

        # compute error
        dg_r = gn_r - g_r
        err = np.sum(dg_r ** 2)
        err = np.sqrt(err / gn_r.size)

        converged = err < tol

        if it % output_rate == 0 or converged:
            fe_tot, _ = compute_free_energy(beta, rho_0, ifft, vl_r, g_r, h_r, cs_k)
            if it == 0:
                logger.info(f"{'iter':<5} {DG_STRING:<12}{'Free Energy':<10} ")
            logger.info(f"{it:<5} {err:>.2e}   {fe_tot:<.7f}")

        if converged:
            logger.info(f"\nReached convergence, {DG_STRING} < {tol}")
            g_r = gn_r
            break

        g_r = diis_update(ndiis, gn_r, dg_r)

    if not converged:
        logger.info(
            f"Could not reach specified convergence criteria after {max_iter} iterations"
        )
        raise ConvergenceError

    logger.info("\n")
    logger.info(f"{'Total Free Energy ':<30} {fe_tot:>12.6f} kj/mol")

    fe_tot, fe_extra = compute_free_energy(beta, rho_0, ifft, vl_r, g_r, h_r, cs_k)

    sim.h_r = h_r
    sim.fe_tot = fe_tot

    if streamer is not None:
        sim.log = streamer.getvalue()
    else:
        sim.log = None
    return fe_tot

def compute_mayer_function(ifft, vs_r, gl_r, gl_k, g_r):
    f_r = np.exp(-vs_r + g_r) - 1.0

    delta_f_r = f_r - gl_r

    delta_f_k = np.apply_along_axis(ifft.to_kspace, 1, delta_f_r)

    f_k = delta_f_k + gl_k

    return f_r, f_k


def compute_h_mol(f_r):
    f1 = f_r + 1.0
    hm_r = np.prod(f1, axis=0) - 1

    return hm_r


def compute_delta_h(ifft, sm, f_r, f_k):
    dd = sm - 1.0

    dd_fk = np.einsum("abn,bn->an", dd, f_k)

    dd_fr = np.apply_along_axis(ifft.to_rspace, 1, dd_fk)

    dh_r = dd_fr * (1.0 + f_r)

    return dh_r


def compute_correlation_hole(ifft, sm, f_r, f_k):
    dd = sm - 1.0

    dd_fk = np.einsum("abn,bn->an", dd, f_k)

    dd_fr = np.apply_along_axis(ifft.to_rspace, 1, dd_fk)

    xi_r = np.sum(f_r, axis=0) + dd_fr - f_r

    return xi_r


def compute_c_k(ifft, sm, h_r, g_r):
    h_k = np.apply_along_axis(ifft.to_kspace, 1, h_r)
    g_k = np.apply_along_axis(ifft.to_kspace, 1, g_r)

    c_k = h_k - np.einsum("abn,bn->an", sm, g_k)

    return c_k


def compute_g_r(ifft, hb, c_k, gl_r):
    hb_c_k = np.einsum("abn,bn->an", hb, c_k)

    hb_c_r = np.apply_along_axis(ifft.to_rspace, 1, hb_c_k)

    g_r = gl_r + hb_c_r

    return g_r


def compute_corr_pot(beta, vl_r, g_r):
    return -(g_r + vl_r) / beta


def compute_gamma_long_range(ifft, hb, sm, vl_k, vl_r, zeta):
    """

    Args:
        ifft: FFT instance
        hb: renormalized correlation function
        sm: molecular structure factor
        vl_k: long range coulomb potential in k-space
        vl_r: long range coulomb potential in k-space
        zeta: long range coulomb potential in k-space

    Returns:

    """

    hb_sm = np.einsum("abn,bcn->acn", hb, sm)

    gl_k = -np.einsum("abn,bn->an", hb_sm, vl_k) - vl_k

    delta_gl_k = gl_k + zeta * vl_k

    gl_r = np.apply_along_axis(ifft.to_rspace, 1, delta_gl_k)

    gl_r = gl_r - zeta * vl_r

    return gl_r, gl_k


def compute_free_energy(beta, rho_0, ifft, vl_r, g_r, h_r, c_k):
    nsites = g_r.shape[0]
    if nsites != 2:
        raise ValueError

    c_r = np.apply_along_axis(ifft.to_rspace, 1, c_k)
    fe0 = -rho_0 * ifft.integrate_rspace(c_r) / nsites / beta

    # alternative way of calculating fe0
    # fe0 = -rho_0 * np.sum(c_k[:, 0]) / nsites / beta

    phi_r = -(g_r + vl_r) / beta

    fe1 = 0.5 * rho_0 * ifft.integrate_rspace(phi_r * h_r)
    fe_tot = fe0 - fe1
    return fe_tot, (fe_tot - fe1, fe1)

