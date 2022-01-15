#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Command line interface to 1D CDFT calculations
"""
import pathlib
import sys
from collections import defaultdict
from contextlib import suppress
from shutil import copy

import click
from click import BadParameter

from cdftpy.cdft1d._version import __version__
from cdftpy.cdft1d.globals import DATA_DIR, EXAMPLES_DIR
from cdftpy.cdft1d.workflow import cdft1d_single_point, cdft1d_multi_solute
from cdftpy.cdft1d.simulation import RUNNERS

KNOWN_METHODS = [k for k in RUNNERS.keys()]

def validate_adjust(ctx, param, value):
    for v in value:
        if v[0] not in ["charge","sigma","eps"]:
            raise BadParameter(F"{value[0]} format must be charge|sigma|eps <float> ")
    return value

def validate_range(ctx, param, value):

    if value is None:
        return value

    msg_par = F"incorrect parameter {value[0]} \n" \
              F" format must be charge|sigma|eps"
    msg_val = F"incorrect parameter {value[1]} \n" \
              F"format must be <float>,<float>,.. or " \
              F" [start:]stop:nsteps "
    if value[0] not in ["charge","sigma","eps"]:
            raise BadParameter(msg_par)

    seq = defaultdict(type(None))

    buffer = value[1]
    if ":" in buffer:
        buffer = buffer.strip(':')
        buffer = buffer.split(':')
        buffer = buffer[::-1]
        try:
            with suppress(IndexError):
                seq["nsteps"] = int(buffer[0])
                seq["stop"] = float(buffer[1])
                seq["start"] = float(buffer[2])
        except ValueError:
            raise BadParameter(msg_val)
    elif "," in buffer:
        buffer = buffer.split(',')
        try:
            float_list = [float(x) for x in buffer if x != '']
        except ValueError:
            raise BadParameter(msg_val)

        seq["values"] = float_list

    if len(seq) == 0:
        raise BadParameter(msg_val)
    return (value[0],seq)

def cdft1d_generate_sample_input():

    input_file = "cl-example"
    for i in range(10):
        input_file = F"example-{i}.dat"
        filepath = pathlib.Path.cwd() / input_file
        if not filepath.exists():
            break
        input_file = None

    if input_file is not None:
        copy(EXAMPLES_DIR / "cl.dat", input_file)

    return input_file

# noinspection PyTypeChecker
@click.command()
@click.argument("input_file", default="")
@click.option("-m", "--method", default=None,
              type=click.Choice(KNOWN_METHODS, case_sensitive=True),
              show_default=False,
              help="Calculation method (default:rsdft)")
@click.option("-s", "--solvent", "solvent_model", default=None, type=str,
              metavar="<solvent model>",
              show_default=True,
              help="solvent model")
@click.option("-o", "--output",  is_flag=True, help="generate solvent density output")
@click.option("-a","--adjust", multiple=True, type=(str, float), callback=validate_adjust,
              metavar='[charge|sigma|eps] <value>',
              help="adjust solute parameters")
@click.option("-r", "--range", "scan", default=None, type=(str, str),
              callback=validate_range,
              metavar='[charge|sigma|eps] "<values>" ',
              help="""Run calculation over the range of solute "charge","sigma","eps" values. Values could specified as
                    comma separated sequence (e.g. "0,0.5,1.0") or in triplets notation [start]:stop:nsteps. To avoid 
                    issues with blank spaces, it is recommended that values are enclosed in double quotes.
                   """)
@click.option("-d", "--dashboard", is_flag=False, flag_value="browser", default=None,
              metavar="[filename]",
              type=click.Path(),
              help="Generate dashboard for analysis. The dashboard will be saved as html file "
                   "under the name provided by optional argument. In the absence of the latter "
                   "dashboard will be open in browser")
@click.option("--version", is_flag=True, help="display version")
def cdft_cli(input_file, method, solvent_model, version, scan, dashboard, adjust, output):
    """
    Perform CDFT calculation

    Args:
        input_file
    """

    if version:
        print(__version__)
        sys.exit(0)

    if input_file == "":
        print('Input file is required to run the calculation\n')
        value = click.prompt('Should I generate one? ', default="N")
        if value == 'y':
            input_file = cdft1d_generate_sample_input()
            print(F"Generated sample input file {input_file}")
            value = click.prompt('Go ahead and run it', default="Y")
            if value != 'Y':
                sys.exit(0)
        else:
            sys.exit(1)

    if output:
        output="density"
    else:
        output=None

    if scan is not None:
        var = scan[0]
        seq = scan[1]
        cdft1d_multi_solute(input_file, var,
                            method=method,
                            solvent=solvent_model,
                            values=seq["values"],
                            stop=seq["stop"],
                            nsteps=seq["nsteps"],
                            start=seq["start"],
                            dashboard=dashboard,
                            output=output
                            )

    else:
        cdft1d_single_point(input_file,
                            method=method,
                            solvent=solvent_model,
                            dashboard=dashboard,
                            adjust=adjust,
                            output=output)







