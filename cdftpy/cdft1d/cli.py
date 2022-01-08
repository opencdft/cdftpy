#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Command line interface to 1D CDFT calculations
"""

import datetime
import json
import pathlib
import sys
from collections import defaultdict
from contextlib import suppress

import click
from click import BadParameter
from prompt_toolkit import prompt
from prompt_toolkit.history import InMemoryHistory
from prompt_toolkit.shortcuts import confirm
from prompt_toolkit.validation import Validator

from cdftpy.cdft1d._version import __version__
from cdftpy.cdft1d.config import DATA_DIR
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


def cdft1d_generate_input():
    def is_float(num):
        try:
            float(num)
            return True
        except ValueError:
            return False

    def is_empty(line):
        return line.strip() != ""

    def file_exist(filename):
        filename = filename + ".dat"
        one_word = filename.strip().find(" ") == -1
        not_blank = filename.strip() != ""
        return one_word and not_blank

    validator = Validator.from_callable(
        is_float, error_message="This input is invalid number", move_cursor_to_end=True
    )
    validator_file = Validator.from_callable(
        file_exist, error_message="Bad filename", move_cursor_to_end=True
    )

    validator_empty = Validator.from_callable(
        is_empty, error_message="input cannot be blank", move_cursor_to_end=True
    )

    history = InMemoryHistory()
    history.append_string("rism")
    history.append_string("rsdft")
    history.append_string("s2_rsdft")
    history.append_string("s2_rism")

    with open(DATA_DIR / "ion_data.json") as fp:
        ion_data = json.load(fp)

    while True:
        try:

            prompt(
                "Simulation type: ",
                accept_default=True,
                default=" Single ion solvation",
            )

            solvent = prompt("Solvent model: ", default="s2")
            print("Provide ion parameters")
            name = prompt("  name: ", validator=validator_empty)

            name = name.lower()
            if name in ion_data:
                default_charge = str(ion_data[name]["charge"])
                default_eps = str(ion_data[name]["eps"])
                default_sigma = str(ion_data[name]["sigma"])
            else:
                default_charge = ""
                default_eps = ""
                default_sigma = ""

            charge = prompt("  charge: ", validator=validator, default=default_charge)
            eps = prompt(
                "  \u03B5 (kj/mol): ", validator=validator, default=default_eps
            )
            sigma = prompt("  \u03C3 (Å): ", validator=validator, default=default_sigma)

            rdf_analysis = prompt("RDF Analysis y/n?:") == "y"
            rdf_output = prompt("RDF Output y/n?: ") == "y"

            while 1:
                input_file = prompt(
                    "Choose name for the input file: ",
                    default=f"{name}.dat",
                    validator=validator_file,
                )
                filepath = pathlib.Path.cwd() / input_file
                if filepath.exists():
                    answer = confirm(
                        f"File {input_file} exists. Overwrite? ", suffix="y/N: "
                    )
                    if answer:
                        break
                    continue
                break

            now = datetime.datetime.now()
            with open(input_file, "w") as fp:
                fp.write("# Input file for single ion solvation calculation\n")
                fp.write(f"# generated on {now.strftime('%Y-%m-%d %H:%M:%S')}\n")
                fp.write("<solute>\n")
                fp.write(f"#name  sigma  eps charge \n")
                fp.write(f" {name}  {sigma} {eps} {charge} \n")
                fp.write("<simulation>\n")
                fp.write(f"solvent {solvent}\n")
                if rdf_analysis:
                    fp.write(f"<analysis>\n")
                    fp.write(f"rdf_peaks\n")
                if rdf_output:
                    fp.write(f"<output>\n")
                    fp.write(f"rdf\n")

            print(f"Generated input file {input_file}")

        except EOFError:
            quit(1)
        except KeyboardInterrupt:
            quit(0)
        else:
            break

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
@click.option("-o", "--output",  is_flag=True, help="generate data output")
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
            input_file = cdft1d_generate_input()
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







