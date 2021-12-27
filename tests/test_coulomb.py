#!/usr/bin/python
# -*- coding: utf-8 -*-
""" Module for testing structure factor"""

import numpy as np
import pytest

from cdftpy.cdft1d.coulomb import compute_long_range_coul_pot_kspace, compute_long_range_coul_pot_rspace, \
    compute_short_range_coul_pot_rspace
from cdftpy.utils.rad_fft import fft_rgrid_iv, fft_kgrid_iv


def test_long_range_coul_pot_kspace():
    qs = np.array([1.0])
    qv = np.array([-0.38, 0.38])
    dr = 0.01
    ngrid = 8192
    rgrid = fft_rgrid_iv(dr, ngrid)
    kgrid = fft_kgrid_iv(dr, ngrid)
    ul = compute_long_range_coul_pot_kspace(qs, qv, kgrid)
    ul_0_0 = -18041957.83284725
    ul_1_1 = 2002360.0013005158
    assert ul[0, 0] == pytest.approx(ul_0_0, abs=1e-12)
    assert ul[1, 1] == pytest.approx(ul_1_1, abs=1e-9)


def test_long_range_coul_pot_rspace():
    qs = np.array([1.0])
    qv = np.array([-0.38, 0.38])
    dr = 0.01
    ngrid = 8192
    rgrid = fft_rgrid_iv(dr, ngrid)
    kgrid = fft_kgrid_iv(dr, ngrid)
    ul = compute_long_range_coul_pot_rspace(qs, qv, rgrid)
    ul_0_0 = -476.58324543962715
    ul_1_1 = 476.5629120887197
    assert ul[0, 0] == pytest.approx(ul_0_0, abs=1e-12)
    assert ul[1, 1] == pytest.approx(ul_1_1, abs=1e-9)


def test_short_range_coul_pot_rspace():
    qs = np.array([1.0])
    qv = np.array([-0.38, 0.38])
    dr = 0.01
    ngrid = 8192
    rgrid = fft_rgrid_iv(dr, ngrid)
    kgrid = fft_kgrid_iv(dr, ngrid)
    ul = compute_short_range_coul_pot_rspace(qs, qv, rgrid)
    ul_0_100 = -134.23515614184834
    assert ul[0, 100] == pytest.approx(ul_0_100, abs=1e-12)
