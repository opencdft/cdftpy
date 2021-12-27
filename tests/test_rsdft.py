#!/usr/bin/python
# -*- coding: utf-8 -*-
""" Module for testing structure factor"""

import pytest

from cdftpy.cdft1d.rdf import analyze_rdf_peaks_sim
from cdftpy.cdft1d.rsdft import rsdft_1d

from cdftpy.cdft1d.solvent import Solvent, solvent_model_locate


def test_rsdft_cation():

    # load solvent model
    solvent_name = "s2"

    filename = solvent_model_locate(solvent_name)

    solvent = Solvent.from_file(filename)

    solute = dict(name="Na", charge=1.0, sigma=2.16, eps=1.4755)

    params = dict(diis_iterations=2, tol=1.0E-7, max_iter=200)

    sim = rsdft_1d(solute, solvent, params=params)

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

    # load solvent model
    solvent_name = "s2"

    filename = solvent_model_locate(solvent_name)

    solvent = Solvent.from_file(filename)

    solute = dict(name="Na", charge=1.0, sigma=2.16, eps=1.4755)

    params = dict(diis_iterations=2, tol=1.0E-7, max_iter=200, rmax=500)

    sim = rsdft_1d(solute, solvent, params=params)

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

    # load solvent model
    solvent_name = "s2"

    filename = solvent_model_locate(solvent_name)

    solv = Solvent.from_file(filename)

    solute = dict(name="Cl", charge=-1.0, sigma=4.83, eps=0.05349244)

    params = dict(diis_iterations=2, tol=1.0E-7, max_iter=200)

    sim = rsdft_1d(solute, solv, params=params)

    fe = sim.fe_tot

    fe_ref = -296.719999
    # fe_ref = -296.3113532110815 old value when using k-space evaluation

    assert fe == pytest.approx(fe_ref, abs=1e-3)

    pk = analyze_rdf_peaks_sim(sim)

    o_peak = pk["O"]['first_peak']
    o_peak_value = round(o_peak[0],2)
    o_peak_pos = round(o_peak[1],2)

    assert o_peak_value == 4.68
    assert o_peak_pos == 3.13
    h_peak = pk["H"]['first_peak']
    h_peak_value = round(h_peak[0],2)
    h_peak_pos = round(h_peak[1],2)

    assert h_peak_value == 8.0
    assert h_peak_pos == 2.16

    pass

def test_rsdft_anion_rmax():

    # load solvent model
    solvent_name = "s2"

    filename = solvent_model_locate(solvent_name)

    solv = Solvent.from_file(filename)

    solute = dict(name="Cl", charge=-1.0, sigma=4.83, eps=0.05349244)

    params = dict(diis_iterations=2, tol=1.0E-7, max_iter=200, rmax=500)

    sim = rsdft_1d(solute, solv, params=params)

    fe = sim.fe_tot

    fe_ref = -297.0210376503721

    assert fe == pytest.approx(fe_ref, abs=1e-3)

    pk = analyze_rdf_peaks_sim(sim)

    o_peak = pk["O"]['first_peak']
    o_peak_value = round(o_peak[0],2)
    o_peak_pos = round(o_peak[1],2)

    assert o_peak_value == 4.68
    assert o_peak_pos == 3.13
    h_peak = pk["H"]['first_peak']
    h_peak_value = round(h_peak[0],2)
    h_peak_pos = round(h_peak[1],2)

    assert h_peak_value == 8.0
    assert h_peak_pos == 2.16

    pass