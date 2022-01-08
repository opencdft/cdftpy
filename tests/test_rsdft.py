#!/usr/bin/python
# -*- coding: utf-8 -*-
""" Module for testing structure factor"""

import pytest

from cdftpy import IonSolvation
from cdftpy.cdft1d.rdf import analyze_rdf_peaks_sim

def test_rsdft_cation():

    solute = dict(name="Na", charge=1.0, sigma=2.16, eps=1.4755)

    params = dict(ndiis=2, tol=1.0E-7, max_iter=200)

    sim = IonSolvation(**solute, **params, solvent="s2", method="rsdft")
    sim.cdft()


    fe = sim.fe_tot
    # fe_ref = -318.3690861027192 old value when using kspace evaluation
    fe_ref = -318.810422

    assert fe == pytest.approx(fe_ref, rel=1e-4)

    pk = analyze_rdf_peaks_sim(sim)

    o_peak = pk["O"]['first_peak']
    o_peak_value = round(o_peak[0],2)
    o_peak_pos = round(o_peak[1],2)

    assert o_peak_value == 6.27
    assert o_peak_pos == 2.46
    h_peak = pk["H"]['first_peak']
    h_peak_value = round(h_peak[0],2)
    h_peak_pos = round(h_peak[1],2)

    assert h_peak_value == 2.99
    assert h_peak_pos == 3.34

def test_rsdft_cation_rmax():

    solute = dict(name="Na", charge=1.0, sigma=2.16, eps=1.4755)

    params = dict(ndiis=2, tol=1.0E-7, max_iter=200, rmax=500)

    sim = IonSolvation(**solute, **params, solvent="s2", method="rsdft")
    sim.cdft()

    fe = sim.fe_tot
    fe_ref = -319.111030

    assert fe == pytest.approx(fe_ref, rel=1e-4)

    pk = analyze_rdf_peaks_sim(sim)

    o_peak = pk["O"]['first_peak']
    o_peak_value = round(o_peak[0],2)
    o_peak_pos = round(o_peak[1],2)

    assert o_peak_value == 6.27
    assert o_peak_pos == 2.46
    h_peak = pk["H"]['first_peak']
    h_peak_value = round(h_peak[0],2)
    h_peak_pos = round(h_peak[1],2)

    assert h_peak_value == 2.99
    assert h_peak_pos == 3.34

def test_rsdft_anion():

    solute = dict(name="Cl", charge=-1.0, sigma=4.83, eps=0.05349244)

    params = dict(ndiis=2, tol=1.0E-7, max_iter=200)

    sim = IonSolvation(**solute, **params, solvent="s2", method="rsdft")
    sim.cdft()

    fe = sim.fe_tot

    # fe_ref = -296.719999 old value with bridge
    # fe_ref = -296.3113532110815 old value when using k-space evaluation
    fe_ref = -301.50062568437266

    assert fe == pytest.approx(fe_ref, abs=1e-3)

    pk = analyze_rdf_peaks_sim(sim)

    o_peak = pk["O"]['first_peak']
    o_peak_value = round(o_peak[0],2)
    o_peak_pos = round(o_peak[1],2)

    # assert o_peak_value == 4.68 old data value using bridge
    assert o_peak_value == 4.78
    # assert o_peak_pos == 3.13 old data value using bridge
    assert o_peak_pos == 3.12
    h_peak = pk["H"]['first_peak']
    h_peak_value = round(h_peak[0],2)
    h_peak_pos = round(h_peak[1],2)

    # assert h_peak_value == 8.0 old data value using bridge
    assert h_peak_value == 8.29
    # assert h_peak_pos == 2.16 old data value using bridge
    assert h_peak_pos == 2.15

    pass

def test_rsdft_anion_rmax():

    solute = dict(name="Cl", charge=-1.0, sigma=4.83, eps=0.05349244)

    params = dict(ndiis=2, tol=1.0E-7, max_iter=200, rmax=500)

    sim = IonSolvation(**solute, **params, solvent="s2", method="rsdft")
    sim.cdft()

    fe = sim.fe_tot

    # fe_ref = -297.0210376503721
    fe_ref = -301.8017110482896

    assert fe == pytest.approx(fe_ref, abs=1e-3)

    pk = analyze_rdf_peaks_sim(sim)

    o_peak = pk["O"]['first_peak']
    o_peak_value = round(o_peak[0],2)
    o_peak_pos = round(o_peak[1],2)

    assert o_peak_value == 4.78
    assert o_peak_pos == 3.12
    h_peak = pk["H"]['first_peak']
    h_peak_value = round(h_peak[0],2)
    h_peak_pos = round(h_peak[1],2)

    assert h_peak_value == 8.29
    assert h_peak_pos == 2.15

    pass