#!/usr/bin/python
# -*- coding: utf-8 -*-
""" Module for testing structure factor"""

import pytest

from cdftpy.cdft1d.rdf import analyze_rdf_peaks_sim
from cdftpy.cdft1d.rism_legacy import rism_1d as rism_1d_legacy
from cdftpy.cdft1d.rism import rism_1d as rism_1d
from cdftpy.cdft1d.solvent import solvent_model_locate, Solvent


def test_rism_legacy():

    # load solvent model
    solvent_name = "s2"
    filename = solvent_model_locate(solvent_name)

    solvent = Solvent.from_file(filename, rism_patch=True)

    solute = dict(name="Na", charge=1.0, sigma=2.16, eps=1.4755)
    params = dict(diis_iterations=2, tol=1.0E-7, max_iter=200)

    sim = rism_1d_legacy(solute, solvent, params=params)
    fe_ref = -315.64395375311017

    assert sim.fe_tot == pytest.approx(fe_ref, rel=1e-9)

def test_rism_cation():

    # load solvent model
    solvent_name = "s2"
    filename = solvent_model_locate(solvent_name)

    solv = Solvent.from_file(filename, rism_patch=True)

    solute = dict(name="Na", charge=1.0, sigma=2.16, eps=1.4755)
    params = dict(diis_iterations=2, tol=1.0E-7, max_iter=200)

    sim = rism_1d(solute, solv, params=params)
    fe_ref = -315.65631175156324

    assert sim.fe_tot == pytest.approx(fe_ref, rel=1e-12)

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

    # load solvent model
    solvent_name = "s2"
    filename = solvent_model_locate(solvent_name)

    solv = Solvent.from_file(filename, rism_patch=True)

    solute = dict(name="Na", charge=1.0, sigma=2.16, eps=1.4755)
    params = dict(diis_iterations=2, tol=1.0E-7, max_iter=200, rmax=500)

    sim = rism_1d(solute, solv, params=params)
    fe_ref = -315.657032

    assert sim.fe_tot == pytest.approx(fe_ref, rel=1e-8)

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

    # load solvent model
    solvent_name = "s2"
    filename = solvent_model_locate(solvent_name)

    solv = Solvent.from_file(filename, rism_patch=True)

    solute = dict(name="Cl", charge=-1.0, sigma=4.83, eps=0.05349244)
    params = dict(diis_iterations=2, tol=1.0E-7, max_iter=200)

    sim = rism_1d(solute, solv, params=params)
    fe_ref = -265.72210193862884

    assert sim.fe_tot == pytest.approx(fe_ref, rel=1e-12)

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
    # load solvent model
    solvent_name = "s2"
    filename = solvent_model_locate(solvent_name)

    solv = Solvent.from_file(filename, rism_patch=True)

    solute = dict(name="Cl", charge=-1.0, sigma=4.83, eps=0.05349244)
    params = dict(diis_iterations=2, tol=1.0E-7, max_iter=200, rmax=500)

    sim = rism_1d(solute, solv, params=params)
    fe_ref = -265.724375

    assert sim.fe_tot == pytest.approx(fe_ref, rel=1e-8)

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