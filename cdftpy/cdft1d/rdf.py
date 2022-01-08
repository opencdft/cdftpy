#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
RDF analysis module
"""
import os
from scipy.signal import argrelmax
from scipy.signal import argrelmin
from cdftpy.cdft1d.io_utils import print_banner


def analyze_rdf_peaks_sim(sim):
    return analyze_rdf_peaks(sim.ifft, sim.solvent.aname, sim.h_r)


def write_rdf_sim(sim, dirpath="./"):
    write_rdf(sim.ifft, sim.name, sim.solvent.aname, sim.h_r, dirpath=dirpath)


def analyze_rdf_peaks(ifft, nametag, h_r):
    rgrid = ifft.rgrid

    peaks = {}
    print_banner("Solvent density structure analysis")
    for i, name in enumerate(nametag):
        imax_array = argrelmax(h_r[i, :], order=10)
        imax0 = imax_array[0][0]
        imax1 = imax_array[0][1]

        imin = argrelmin(h_r[i, :], order=10)[0]
        imin_filter = imin > imax0
        imin0 = imin[imin_filter][0]

        peaks[name] = dict(
            first_peak=(h_r[i, imax0] + 1, rgrid[imax0]),
            second_peak=(h_r[i, imax1] + 1, rgrid[imax1]),
            first_min=(h_r[i, imin0] + 1, rgrid[imin0])
        )
        print(
            f"{name} 1st peak position/height: {rgrid[imax0]:<.2f} {h_r[i, imax0] + 1:<.3f}  "
        )
        print(
            f"{name} 2nd peak position/height: {rgrid[imax1]:<.2f} {h_r[i, imax1] + 1:<.3f}  "
        )
        print(
            f"{name} 1st min position/height: {rgrid[imin0]:<.2f} {h_r[i, imin0] + 1:<.3f}  "
        )
        print("  ")

    return peaks


def write_rdf(ifft, solute_name, solvent_name, h_r, dirpath="./"):
    rgrid = ifft.rgrid
    print("\nGenerating solvent-solute RDF's")
    for i in range(len(solvent_name)):
        filename = os.path.join(dirpath, f"rdf_{solute_name}{solvent_name[i]}.dat")
        print(f"Generating {os.path.abspath(filename)}")
        with open(filename, "w") as fp:
            fp.write(f"{'# r':<10} {'g(r)':<10}\n")

            for j in range(len(rgrid)):
                r = rgrid[j]
                rdf = h_r[i, j] - h_r[i, 0]
                fp.write(f"{r:<10.4} {rdf:<10.4}\n")
