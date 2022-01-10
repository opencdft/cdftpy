#!/usr/bin/python
# -*- coding: utf-8 -*-
""" Module for testing structure factor"""

import pytest

from cdftpy import IonSolvation
from cdftpy.cdft1d.rdf import analyze_rdf_peaks_sim

def test_rism_cation():

    solute = dict(name="Na", charge=1.0, sigma=2.16, eps=1.4755)
    params = dict(ndiis=2, tol=1.0E-7, max_iter=200)

    sim = IonSolvation(**solute, **params, solvent="s2", method="rism")
    sim.cdft()

    fe_ref = -315.6560905213441

    assert sim.fe_tot == pytest.approx(fe_ref, abs=1e-3)

    pk = analyze_rdf_peaks_sim(sim)

    o_peak = pk["O"]['first_peak']
    o_peak_value = round(o_peak[0],2)
    o_peak_pos = round(o_peak[1],2)

    assert o_peak_value == 6.04
    assert o_peak_pos == 2.44
    h_peak = pk["H"]['first_peak']
    h_peak_value = round(h_peak[0],2)
    h_peak_pos = round(h_peak[1],2)

    assert h_peak_value == 1.66
    assert h_peak_pos == 3.34

def test_rism_cation_rmax():

    solute = dict(name="Na", charge=1.0, sigma=2.16, eps=1.4755)
    params = dict(ndiis=2, tol=1.0E-7, max_iter=200, rmax=500)

    sim = IonSolvation(**solute, **params, solvent="s2", method="rism")
    sim.cdft()

    fe_ref = -315.65718693067294

    assert sim.fe_tot == pytest.approx(fe_ref, abs=1e-3)

    pk = analyze_rdf_peaks_sim(sim)

    o_peak = pk["O"]['first_peak']
    o_peak_value = round(o_peak[0],2)
    o_peak_pos = round(o_peak[1],2)

    assert o_peak_value == 6.04
    assert o_peak_pos == 2.44
    h_peak = pk["H"]['first_peak']
    h_peak_value = round(h_peak[0],2)
    h_peak_pos = round(h_peak[1],2)

    assert h_peak_value == 1.66
    assert h_peak_pos == 3.34

def test_rism_anion():

    solute = dict(name="Cl", charge=-1.0, sigma=4.83, eps=0.05349244)
    params = dict(ndiis=2, tol=1.0E-7, max_iter=200)

    sim = IonSolvation(**solute, **params, solvent="s2", method="rism")
    sim.cdft()
    fe_ref = -265.72289940523075

    assert sim.fe_tot == pytest.approx(fe_ref, abs=1e-3)

    pk = analyze_rdf_peaks_sim(sim)

    o_peak = pk["O"]['first_peak']
    o_peak_value = round(o_peak[0],2)
    o_peak_pos = round(o_peak[1],2)

    assert o_peak_value == 3.11
    assert o_peak_pos == 3.4
    h_peak = pk["H"]['first_peak']
    h_peak_value = round(h_peak[0],2)
    h_peak_pos = round(h_peak[1],2)

    assert h_peak_value == 6.16
    assert h_peak_pos == 2.25

def test_rism_anion_rmax():

    solute = dict(name="Cl", charge=-1.0, sigma=4.83, eps=0.05349244)
    params = dict(ndiis=2, tol=1.0E-7, max_iter=200, rmax=500)

    sim = IonSolvation(**solute, **params, solvent="s2", method="rism")
    sim.cdft()
    fe_ref = -265.72446784728345

    assert sim.fe_tot == pytest.approx(fe_ref, abs=1e-3)

    pk = analyze_rdf_peaks_sim(sim)

    o_peak = pk["O"]['first_peak']
    o_peak_value = round(o_peak[0], 2)
    o_peak_pos = round(o_peak[1], 2)

    assert o_peak_value == 3.11
    assert o_peak_pos == 3.4
    h_peak = pk["H"]['first_peak']
    h_peak_value = round(h_peak[0], 2)
    h_peak_pos = round(h_peak[1], 2)

    assert h_peak_value == 6.16
    assert h_peak_pos == 2.25