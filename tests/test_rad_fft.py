# -*- coding: utf-8 -*-
from cdftpy.utils.rad_fft import RadFFT, fft_kgrid_spacing

""" Module for testing discrete sine transform"""

import math
import numpy as np
import pytest


def fft_pair(rg, kg, a=1.0):
    pi = math.pi

    fr = np.exp(-a * rg ** 2)

    fk = math.sqrt(pi / a) ** 3 * np.exp(-kg ** 2 / (4.0 * a))

    return fr, fk


def test_fft_grid():
    dr = 0.02
    ntot = 1000

    dk = fft_kgrid_spacing(dr, ntot)

    ifft = RadFFT(dr, ntot)

    # note that this is a mid-point based grid r(i) = dr(i-0.5), i=1,N
    rg_max_ref = (ntot - 0.5) * dr
    kg_max_ref = (ntot - 0.5) * dk

    assert ifft.rgrid[-1] == pytest.approx(rg_max_ref, rel=1e-10)
    assert ifft.kgrid[-1] == pytest.approx(kg_max_ref, rel=1e-10)


def test_fft_to_kspace():
    dr = 0.01
    rmax = 40

    ntot = rmax / dr

    ifft = RadFFT(dr, ntot)

    fr, fk0 = fft_pair(ifft.rgrid, ifft.kgrid, a=2.0)

    fk = ifft.to_kspace(fr)

    assert fk == pytest.approx(fk0, rel=1e-12)


def test_fft_to_rspace():
    dr = 0.01
    rmax = 40

    ntot = rmax / dr

    ifft = RadFFT(dr, ntot)

    fr0, fk = fft_pair(ifft.rgrid, ifft.kgrid, a=2.0)

    fr = ifft.to_rspace(fk)

    assert fr == pytest.approx(fr0, rel=1e-12)
