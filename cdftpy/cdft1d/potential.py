#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
This module provides collection
of routines for Lennard-Jones potential calculations
"""
import numpy as np


def compute_lj_potential(sig_solute, eps_solute, sig_solvent, eps_solvent, rgrid):

    eps_solute = float(eps_solute)
    sig_solute = float(sig_solute)

    eps = eps_solute * eps_solvent
    eps = np.sqrt(eps)

    sig = -0.5 * (sig_solute + sig_solvent)

    sig6 = sig ** 6

    r_m6 = 1.0 / rgrid ** 6
    r_m12 = 1.0 / rgrid ** 12
    sig12 = sig ** 12

    lj_repuls = 4.0 * np.multiply.outer(eps * sig12, r_m12)
    lj_attrac = -4.0 * np.multiply.outer(eps * sig6, r_m6)

    lj = lj_repuls + lj_attrac
    return lj


def compute_lj_potential_mod(sig_solute, eps_solute, sig_solvent, eps_solvent, rgrid):

    eps_solute = float(eps_solute)
    sig_solute = float(sig_solute)

    eps = eps_solute * eps_solvent
    eps = np.sqrt(eps)

    sig = -0.5 * (sig_solute + sig_solvent)

    sig6 = sig ** 6
    sig6[1] = 0

    r_m6 = 1.0 / rgrid ** 6
    r_m12 = 1.0 / rgrid ** 12
    sig12 = sig ** 12

    lj_repuls = 4.0 * np.multiply.outer(eps * sig12, r_m12)
    lj_attrac = -4.0 * np.multiply.outer(eps * sig6, r_m6)

    lj = lj_repuls + lj_attrac
    return lj
