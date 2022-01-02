#!/usr/bin/python
# -*- coding: utf-8 -*-
""" Module for testing structure factor"""

import pytest

from cdftpy.cdft1d.workflow import cdft1d_single_point, cdft1d_multi_solute


def test_single_point():
    input_file = "data/cl.dat"
    method = "rsdft"
    solvent_model = "s2"
    fe = cdft1d_single_point(input_file, method, solvent_model)
    fe_ref = -296.7836464403522
    assert fe == pytest.approx(fe_ref, abs=1e-6)


def test_single_point_adjust():
    input_file = "data/cl.dat"
    method = "rsdft"
    solvent_model = "s2"
    adjust = [("charge", "0")]
    fe = cdft1d_single_point(input_file, method, solvent_model, adjust=adjust)
    fe_ref = 71.16765356431705
    assert fe == pytest.approx(fe_ref, abs=1e-6)


def test_multi_point():
    input_file = "data/cl.dat"
    method = "rsdft"
    solvent_model = "s2"
    var = "charge"

    fe = cdft1d_multi_solute(input_file, method, solvent_model,
                        var,
                        values=[-0.5,0],
                        )

    fe_ref = [-14.061937854150358,71.16765356431705]
    assert fe == pytest.approx(fe_ref, abs=1e-6)

    fe = cdft1d_multi_solute(input_file, method, solvent_model,
                        var,
                        values=[-0.5],
                        )

    fe_ref = [-14.061937854150358]
    assert fe == pytest.approx(fe_ref, abs=1e-6)

def test_multi_point_triplet():

    input_file = "data/cl.dat"
    method = "rsdft"
    solvent_model = "s2"
    var = "charge"

    fe = cdft1d_multi_solute(input_file, method, solvent_model,
                             var,
                             start=-1,
                             stop=0,
                             nsteps=2
                             )

    fe_ref = [-296.7836464403522, 71.16765356431705]
    assert fe == pytest.approx(fe_ref, abs=1e-6)