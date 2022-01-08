#!/usr/bin/python
# -*- coding: utf-8 -*-
""" Module for testing structure factor"""


import pytest

from cdftpy.cdft1d.workflow import cdft1d_single_point, cdft1d_multi_solute

input_file_data = """
# Cl- solvation in s2 solvent
<solute>
# site   sigma        eps(kj/mol)    charge(e)   
Cl       4.83         0.05349244     -1.0        
<simulation>
solvent n2
method rism
tol 1.0E-9
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
    fe = cdft1d_single_point(input_file, method = method, solvent = solvent_model)
    # fe_ref = -296.7836464403522
    fe_ref = -301.56429559691776
    assert fe == pytest.approx(fe_ref, abs=1e-3)

# def test_single_point_rism(input_file):
#
#     solvent_model = "s2"
#     fe = cdft1d_single_point(input_file, solvent_model)
#     fe_ref = -296.7836464403522
#     assert fe == pytest.approx(fe_ref, abs=1e-3)


def test_single_point_adjust(input_file):
    method = "rsdft"
    solvent_model = "s2"
    adjust = [("charge", "0")]
    fe = cdft1d_single_point(input_file, method = method, solvent = solvent_model, adjust=adjust)
    # fe_ref = 71.16765356431705
    fe_ref = 70.55222234327249
    assert fe == pytest.approx(fe_ref, abs=1e-3)


def test_multi_point(input_file):

    method = "rsdft"
    solvent = "s2"
    var = "charge"

    fe = cdft1d_multi_solute(input_file, var, method=method, solvent=solvent,
                        values=[-0.5,0],
                        )

    # fe_ref = [-14.061937854150358,71.16765356431705]
    fe_ref = [-15.812939268973821, 70.55222234327249]
    assert fe == pytest.approx(fe_ref, abs=1e-3)

    fe = cdft1d_multi_solute(input_file, var, method=method, solvent=solvent,
                        values=[-0.5],
                        )

    fe_ref = [-14.061937854150358]
    fe_ref = [-15.812939268973821]
    assert fe == pytest.approx(fe_ref, abs=1e-3)

def test_multi_point_triplet(input_file):


    method = "rsdft"
    solvent = "s2"
    var = "charge"

    fe = cdft1d_multi_solute(input_file, var,
                             method=method, solvent=solvent,
                             start=-1,
                             stop=0,
                             nsteps=2
                             )

    # fe_ref = [-296.7836464403522, 71.16765356431705]
    fe_ref = [-301.56429559691776, 70.55222234327249]
    assert fe == pytest.approx(fe_ref, abs=1e-3)