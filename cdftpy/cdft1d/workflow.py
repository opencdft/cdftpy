#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Task wrappers to CDFT calculations
"""

import math
import sys

import numpy as np
from prettytable import PrettyTable, PLAIN_COLUMNS

from cdftpy.cdft1d.io_utils import read_key_value, print_banner
from cdftpy.cdft1d.io_utils import read_solute
from cdftpy.cdft1d.rdf import analyze_rdf_peaks_sim
from cdftpy.cdft1d.rdf import write_rdf_sim
from cdftpy.cdft1d.rism import rism_1d
from cdftpy.cdft1d.rsdft import rsdft_1d
from cdftpy.cdft1d.solvent import solvent_model_locate, Solvent
from cdftpy.cdft1d.exceptions import ConvergenceError
from cdftpy.cdft1d.viz import single_point_viz, multi_solute_viz

_RUNNERS = dict(rism=rism_1d, rsdft=rsdft_1d)


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

    parameters = read_key_value(input_file, section="simulation")

    parameters["solvent"] = solvent_model
    solvent_name = parameters["solvent"]
    filename = solvent_model_locate(solvent_name)

    rism_patch = (method == "rism")
    solvent = Solvent.from_file(filename, rism_patch=rism_patch)

    runner = _RUNNERS[method]
    sim = []

    print_level = {'header', 'parameters','solute','solvent'}
    for v in values:
        solute[var] = v
        try:
            s = runner(solute, solvent, params=parameters,
                       print_level=print_level)
            print("\n")
        except ConvergenceError as e:
            print(F"cannot converge {var}={v} point")
            print("skipping the rest of the cycle")
            break
        print_level = {'solute'}
        sim.append(s)

    if len(sim) > 1:
        tbl = PrettyTable()

        tbl.set_style(PLAIN_COLUMNS)
        tbl.field_names = [var.capitalize(), "Solvation Free Energy",
                           "Solvation Free Energy "]
        tbl.add_row([" ", "total (kj/mol)", "diff(kj/mol)"])
        fe_ref = sim[0].fe_tot
        for v, s in zip(values, sim):
            fe_tot = s.fe_tot
            tbl.add_row([v,fe_tot,fe_tot-fe_ref])
        tbl.align = "r"
        tbl.float_format = ".3"

        print_banner(F"Final results:")

        print(tbl)

        if dashboard is not None:

            multi_solute_viz(var, values[:len(sim)], sim, dashboard_dest=dashboard)

    else:
        analysis = read_key_value(input_file, section="analysis")
        if analysis is not None:
            if "rdf_peaks" in analysis:
                analyze_rdf_peaks_sim(sim[0])

        output = read_key_value(input_file, section="output")
        if output is not None:
            if "rdf" in output:
                write_rdf_sim(sim[0])

        if dashboard is not None:
            single_point_viz(sim[0], dashboard_dest=dashboard)

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

    parameters["solvent"] = solvent_model
    solvent_name = parameters["solvent"]
    filename = solvent_model_locate(solvent_name)

    if method == "rism":
        solvent = Solvent.from_file(filename, rism_patch=True)
        sim = rism_1d(solute, solvent, params=parameters)
    elif method == "rsdft":
        solvent = Solvent.from_file(filename, rism_patch=False)
        sim = rsdft_1d(solute, solvent, params=parameters)
    else:
        print(f"Unknown method {theory}")
        sys.exit(1)

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

