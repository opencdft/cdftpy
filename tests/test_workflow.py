#!/usr/bin/python
# -*- coding: utf-8 -*-
""" Module for testing structure factor"""
import json

import pytest

from cdftpy.cdft1d.workflow import cdft1d_single_point, cdft1d_multi_solute

input_file_data = """
# Cl- solvation in s2 solvent
<solute>
# site   sigma        eps(kj/mol)    charge(e)   x   y   z
Cl       4.83         0.05349244     -1.0         0.0 0.0 0.0
<simulation>
tol 1.0E-7
max_iter 500
rmax 100
#optional analysis
<analysis>
rdf_peaks
<output>
rdf
"""

@pytest.fixture
def input_file(tmpdir):

    filename = tmpdir.join('cl.dat')

    with open(filename, 'w') as fh:
        fh.write(input_file_data)
    return str(filename)


def test_single_point(input_file):
    method = "rsdft"
    solvent_model = "s2"
    fe = cdft1d_single_point(input_file, method, solvent_model)
    fe_ref = -296.7836464403522
    assert fe == pytest.approx(fe_ref, abs=1e-6)


def test_single_point_adjust(input_file):
    method = "rsdft"
    solvent_model = "s2"
    adjust = [("charge", "0")]
    fe = cdft1d_single_point(input_file, method, solvent_model, adjust=adjust)
    fe_ref = 71.16765356431705
    assert fe == pytest.approx(fe_ref, abs=1e-6)


def test_multi_point(input_file):

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

def test_multi_point_triplet(input_file):


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