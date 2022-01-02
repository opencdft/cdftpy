#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Task wrappers to CDFT calculations
"""
import copy
import math
import sys

import numpy as np
from prettytable import PrettyTable, PLAIN_COLUMNS

from cdftpy.cdft1d.io_utils import read_key_value, print_banner
from cdftpy.cdft1d.io_utils import read_solute
from cdftpy.cdft1d.simulation import SolvatedIon
from cdftpy.cdft1d.rdf import analyze_rdf_peaks_sim
from cdftpy.cdft1d.rdf import write_rdf_sim
from cdftpy.cdft1d.exceptions import ConvergenceError
from cdftpy.cdft1d.viz import single_point_viz, multi_solute_viz


def my_linspace(start, stop, nsteps):
    values, dv = np.linspace(start, stop, num=nsteps, retstep=True)

    if dv < 1:
        nr = int(-math.log10(dv)) + 1
    else:
        nr = 1
    return list(map(lambda x: x.round(nr), values))


def cdft1d_multi_solute(input_file, method, solvent_model, var,
                        values=None,
                        stop=None,
                        nsteps=11,
                        start=None,
                        dashboard=None
                        ):
    try:
        solute = read_solute(input_file)
    except FileNotFoundError:
        print(f"Cannot locate input file {input_file}")
        sys.exit(1)

    for k, v in solute.items():
        solute[k] = v[0]

    if values is None:
        if stop is None:
            print("Neither values or stop is provided")
            sys.exit(1)

        if start is None:
            start = solute[var]

        values = my_linspace(start, stop, nsteps)

    if len(values) == 1:
        adjust = [(var,values[0])]
        fe = cdft1d_single_point(input_file,
                            method, solvent_model,
                            dashboard=dashboard, adjust=adjust)
        return [fe]

    parameters = read_key_value(input_file, section="simulation")

    sim_array = []
    fe_array = []

    sim = SolvatedIon(solute=solute, solvent=solvent_model, params=parameters, method=method)

    print_level = {'header', 'parameters', 'solute', 'solvent'}
    print(sim.to_string(print_level=print_level))
    print(F"Performing calculation over the range of {var}s \n {values}")

    for v in values:
        sim.solute[var] = v
        try:
            print(F"{var}={v}")
            fe = sim.cdft()
            fe_array.append(fe)
            print(sim.log)

        except ConvergenceError:
            print(F"cannot converge {var}={v} point")
            print("skipping the rest of the cycle")
            break

        sim_array.append(copy.deepcopy(sim))

    tbl = PrettyTable()
    tbl.set_style(PLAIN_COLUMNS)
    tbl.field_names = [var.capitalize(), "Solvation Free Energy",
                       "Solvation Free Energy "]
    tbl.add_row([" ", "total (kj/mol)", "diff(kj/mol)"])
    fe_ref = fe_array[0]
    for fe,v in zip(fe_array,values):
        tbl.add_row([v, fe, fe - fe_ref])
    tbl.align = "r"
    tbl.float_format = ".3"

    print_banner(F"         Final results:          ")
    print(tbl)

    if dashboard is not None:
        multi_solute_viz(var, values[:len(sim_array)], sim_array, dashboard_dest=dashboard)

    return fe_array


def cdft1d_single_point(input_file, method, solvent_model, dashboard=None, adjust=None):

    try:
        solute = read_solute(input_file)
    except FileNotFoundError:
        print(f"Cannot locate input file {input_file}")
        sys.exit(1)

    for k, v in solute.items():
        solute[k] = v[0]

    if adjust is not None:
        for par_val in adjust:
            par, val = par_val
            solute[par] = float(val)

    parameters = read_key_value(input_file, section="simulation")

    sim = SolvatedIon(solute=solute, solvent=solvent_model, params=parameters, method=method)
    print(sim.to_string())

    fe = sim.cdft(quiet=False)

    analysis = read_key_value(input_file, section="analysis")
    if analysis is not None:
        if "rdf_peaks" in analysis:
            analyze_rdf_peaks_sim(sim)

    output = read_key_value(input_file, section="output")
    if output is not None:
        if "rdf" in output:
            write_rdf_sim(sim)

    if dashboard is not None:
        single_point_viz(sim, dashboard_dest=dashboard)

    return fe
