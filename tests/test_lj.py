#!/usr/bin/python
# -*- coding: utf-8 -*-
""" Module for testing structure factor"""

import numpy as np
import pytest

from cdftpy.cdft1d.potential import compute_lj_potential, compute_lj_potential_mod
from cdftpy.utils.rad_fft import fft_rgrid_iv


def test_compute_lj_matrix():

    sig_s = np.array([2.16])
    sig_v = np.array([3.16, 1 ])
    eps_s = np.array([1.4755])
    eps_v = np.array([0.648954, 0.19260])

    dr = 0.01
    ngrid = 8192
    rgrid = fft_rgrid_iv(dr, ngrid)

    lj = compute_lj_potential(sig_s, eps_s, sig_v, eps_v, rgrid)

    lj_0_100 = 461273.34383968828

    assert lj[0,100] == pytest.approx(lj_0_100, abs=1e-8)

    lj = compute_lj_potential_mod(sig_s, eps_s, sig_v, eps_v, rgrid)

    lj_1_100 = 486.1273348004878

    assert lj[1,100] == pytest.approx(lj_1_100, abs=1e-8)