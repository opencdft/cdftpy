#!/usr/bin/python
# -*- coding: utf-8 -*-
""" Module for testing structure factor"""

import numpy as np
import pytest

from cdftpy.cdft1d.coulomb import compute_long_range_coul_pot_kspace, compute_long_range_coul_pot_rspace, \
    compute_short_range_coul_pot_rspace
from cdftpy.cdft1d.rad_fft import fft_rgrid_iv, fft_kgrid_iv


def test_long_range_coul_pot_kspace():
    qs = np.array([1.0])
    qv = np.array([-0.38, 0.38])
    dr = 0.01
    ngrid = 8192
    rgrid = fft_rgrid_iv(dr, ngrid)
    kgrid = fft_kgrid_iv(dr, ngrid)
    ul = compute_long_range_coul_pot_kspace(qs, qv, kgrid)
    ul_0_0 = -18041984.928852376
    ul_1_1 = 2002363.0085105628
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
    ul_0_0 = -476.58396118800556
    ul_1_1 = 476.5636278065608
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
    ul_0_100 = -134.23535774061676
    assert ul[0, 100] == pytest.approx(ul_0_100, abs=1e-12)
