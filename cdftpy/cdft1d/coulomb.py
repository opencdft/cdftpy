#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
This module provides collection
of routines for the calculations of Coulomb interactions
"""
from math import pi

import numpy as np
from scipy.special import erf

from cdftpy.utils.units import bohr_2_ang
from cdftpy.utils.units import epot_2_kjmol
from cdftpy.utils.units import hartree_2_kjmol

_unit_factor = {"hartree": 1.0 / bohr_2_ang, "kj/mol": hartree_2_kjmol / bohr_2_ang}


def compute_coulomb_potential(rho_0, qv, ifft, htot_k):
    kgrid = ifft.kgrid
    epot_k = rho_0 * np.einsum("an,n->an", htot_k, 4.0 * np.pi / kgrid ** 2)
    epot_r = np.apply_along_axis(ifft.to_rspace, 1, epot_k)
    epot_r = np.einsum("a,an->an", qv, epot_r)
    return epot_r, epot_k


def compute_coulomb_energy(qs, epot_r):
    return 0.5 * qs * epot_2_kjmol * np.sum(epot_r[:, 0])


def compute_gauss_pot_kspace(kgrid, r_s):
    r"""
    Generates Coulomb potential of gaussian smeared density

    .. math::
        u(k) = \frac{4\pi}{k^2} e^{-r_s^2k^2/4}

    Args:
        kgrid: kspace grid

        r_s: smear_factor

    Returns:
        :math:`u(k)` - 1D numpy array of potential

    """
    l2 = r_s ** 2
    k2 = kgrid ** 2
    v = 4.0 * pi * np.exp(-l2 * k2 / 4) / k2

    return v


def compute_gauss_pot_rspace(rgrid, r_s):
    r"""
    Generates Coulomb potential of gaussian smeared density

    .. math::
        u(r) = \frac{erf(r/r_s)}{r}

    Args:
        rgrid: rspace grid

        r_s: smear_factor

    Returns:
        :math:`u(r)` - 1D numpy array of potential

    """
    if r_s != 0:
        v = erf(rgrid / r_s) / rgrid
    else:
        v = 1.0 / rgrid

    return v


def compute_compl_gauss_pot_rspace(rgrid, r_s):
    r"""
    Generates Coulomb potential difference generated
    by point charge and gaussian smeared density

    .. math::
        u(r) = \left(\frac{1}{r} - \frac{erf(r/r_s)}{r}\right)

    Args:
        rgrid: rspace grid

        r_s: smear_factor

    Returns:
        :math:`u(r)` - 1D numpy array of potential

    """
    v = 1 / rgrid - erf(rgrid / r_s) / rgrid

    return v


def compute_long_range_coul_pot_kspace(qs, qv, kgrid, r_s=1.25, units="kj/mol"):
    r"""
    Generates long range Coulomb potential in reciprocal space created
    by **single** solute charge

    .. math::
        u_\alpha(k) =  q_{v\alpha}\frac{4\pi q_s}{k^2} e^{-r_s^2k^2/4}


    Args:
        qs: solute charge

        qv: array of solvent charges

        kgrid: kspace grid

        r_s: smear_factor

        units: **kj/mol**  or hartree
    Returns:
        :math:`u_\alpha(k)` - 2D numpy array of potential

    """

    # matrix of solute-solvent charges
    zab = float(qs) * qv

    v = compute_gauss_pot_kspace(kgrid, r_s)
    ul = np.multiply.outer(zab, v) * _unit_factor[units]

    return ul


def compute_long_range_coul_pot_rspace(qs, qv, rgrid, r_s=1.25, units="kj/mol"):
    r"""
    Generates long range Coulomb potential in real space created
    by **single** solute charge

    .. math::
        u_\alpha(r) =  q_{v\alpha} \frac{q_s erf(r/r_s)}{r}


    Args:
        qs: solute charge

        qv: array of solvent charges

        rgrid: rspace grid

        r_s: smear_factor

        units: **kj/mol**  or hartree
    Returns:
        :math:`u_\alpha(r)` - 2D numpy array of potential

    """

    # matrix of solute-solvent charges
    zab = float(qs) * qv

    v = compute_gauss_pot_rspace(rgrid, r_s)
    ul = np.multiply.outer(zab, v) * _unit_factor[units]

    return ul


def compute_short_range_coul_pot_rspace(qs, qv, rgrid, r_s=1.25, units="kj/mol"):
    r"""
    Generates short range Coulomb potential in real space created
    by **single** solute charge

    .. math::
        u_\alpha(r) =  q_{v\alpha} \left(\frac{q_s}{r} - \frac{q_s erf(r/r_s)}{r}\right)


    Args:
        qs: solute charge

        qv: array of solvent charges

        rgrid: rspace grid

        r_s: smear_factor

        units: **kj/mol**  or hartree
    Returns:
        :math:`u_\alpha(r)` - 2D numpy array of potential

    """

    # matrix of solute-solvent charges
    zab = float(qs) * qv

    v = compute_compl_gauss_pot_rspace(rgrid, r_s)
    ul = np.multiply.outer(zab, v) * _unit_factor[units]

    return ul
