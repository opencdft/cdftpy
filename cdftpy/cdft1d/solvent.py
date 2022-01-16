#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
This module provides collection
of routines related to solvent
"""
import copy
import logging
import os
import pathlib
import sys
from dataclasses import dataclass, field
from io import StringIO
from itertools import combinations as comb

import numpy as np
from prettytable import PrettyTable, PLAIN_COLUMNS

from cdftpy.cdft1d.globals import DATA_DIR
from cdftpy.cdft1d.io_utils import read_array
from cdftpy.cdft1d.io_utils import read_key_value
from cdftpy.cdft1d.io_utils import read_molecule
from cdftpy.cdft1d.rad_fft import RadFFT, fft_rgrid_iv
from cdftpy.cdft1d.units import hartree_2_kjmol, bohr_2_ang

logging.basicConfig(level=logging.INFO)


@dataclass
class Molecule1:
    x: list[float] = field(default_factory=list)
    y: list[float] = field(default_factory=list)
    z: list[float] = field(default_factory=list)

    sigma: list[float] = field(default_factory=list)
    eps: list[float] = field(default_factory=list)
    charge: list[float] = field(default_factory=list)

    aname: list[str] = field(default_factory=list)
    atype: list[str] = field(default_factory=list)

    dipole: float = field(init=False)
    nv: int = field(init=False)

    def __post_init__(self):
        nv = len(self.x)
        if nv != 0:
            self.nv = nv
            q = self.charge
            mu2 = np.dot(self.x, q) ** 2 + np.dot(self.y, q) ** 2 + np.dot(self.z, q) ** 2
            self.dipole = np.sqrt(mu2)

    @classmethod
    def from_file(cls, filename):
        solvent = read_molecule(filename)
        solvent["nv"] = len(solvent["aname"])

        return cls(**solvent)

    def sm_k(self, kgrid):
        return compute_sm_rigid_bond(self.x, self.y, self.z, kgrid)


@dataclass
class Molecule:
    nv: int
    x: list[float]
    y: list[float]
    z: list[float]

    sigma: list[float]
    eps: list[float]
    charge: list[float]

    aname: list[str]
    atype: list[str]

    dipole: float = field(init=False)

    def __post_init__(self):
        q = self.charge
        mu2 = np.dot(self.x, q) ** 2 + np.dot(self.y, q) ** 2 + np.dot(self.z, q) ** 2
        self.dipole = np.sqrt(mu2)

    @classmethod
    def from_file(cls, filename):
        solvent = read_system(filename)
        solvent["nv"] = len(solvent["aname"])

        return cls(**solvent)

    def sm_k(self, kgrid):
        return compute_sm_rigid_bond(self.x, self.y, self.z, kgrid)


@dataclass
class Solvent(Molecule):
    model: str
    filename: str
    file_location: str

    density: float
    temp: float

    dielectric: float = None

    ifft: object = None
    s_k: np.ndarray = None
    hbar_k: np.ndarray = None


    @classmethod
    def from_file(cls, filename, rism_patch=False):
        solvent = read_molecule(filename, rism_patch=rism_patch)
        nv = len(solvent["aname"])
        solvent["nv"] = nv
        solvent["filename"] = os.path.abspath(filename)
        solvent["file_location"] = os.path.split(solvent["filename"][0])
        solvent["model"] = os.path.basename(filename).split(".")[0]

        state = read_key_value(filename, "state")

        # structure factors
        hbar, kgrid = read_array(filename, "hbar_k")
        s_k, kgrid = read_array(filename, "s_k")
        solvent["s_k"] = s_k
        solvent["hbar_k"] = hbar
        if kgrid is not None:
            solvent["ifft"] = RadFFT.from_kgrid(kgrid)
        return cls(**solvent, **state)

    @classmethod
    def from_name(cls,name,rism_patch=False):
        filepath = solvent_model_locate(name)
        return cls.from_file(filepath,rism_patch=rism_patch)

    @property
    def rgrid(self):
        return self.ifft.rgrid

    def hbar(self,method):
        s_k = self.s_k
        sm = self.sm_k(self.kgrid)
        if method == "rsdft":
            inv_sm = np.linalg.inv(sm.transpose(2, 0, 1))
            delta_s_k = s_k - sm
            hbar = np.einsum('imk,kmj->ijk', np.einsum('kij,jmk->imk',
                                                       inv_sm, delta_s_k), inv_sm)
        elif method == "rism":
            delta = np.diagflat(np.ones(s_k.shape[0]))
            hbar = s_k - np.expand_dims(delta, axis=2)
        else:
            print(F"incorrect calculation method {method}")
            sys.exit(0)

        return hbar

    @property
    def kgrid(self):
        return self.ifft.kgrid

    def zeta(self, beta):
        mu2 = self.dipole ** 2
        if mu2 == 0:
            return 0
        rho = self.density
        y3 = (hartree_2_kjmol / bohr_2_ang) * beta * 4 * np.pi * mu2 * rho / 3
        y3m = 1.0 - 1.0 / self.dielectric
        zeta = y3m / y3
        return zeta

    def to_smdl_file(self, filename, structure_factor):
        with open(filename, "w") as fp:
            fp.write("<geometry>\n")
            fp.write("# site   type   x     y     z\n")
            for i in range(self.nv):
                fp.write(F"{self.aname[i]}   {self.atype[i]}    {self.x[i]}     {self.y[i]}     {self.z[i]}\n")

            fp.write("<interactions>\n")
            fp.write(f"#{'type' : <10}{'sigma' : ^10}{'eps' : ^10}{'charge' : >10}\n")
            for i in range(self.nv):
                fp.write(F" {self.atype[i]: <10}{self.sigma[i]: ^10}{self.eps[i]: ^10}{self.charge[i]: >10}\n")

            fp.write("<state>\n")
            fp.write(F"density {self.density}\n")
            fp.write(F"temp {self.temp}\n")
            fp.write(F"dielectric {self.dielectric}\n")

            s_k = structure_factor["s_k"]
            kgrid = structure_factor["kgrid"]
            tag = structure_factor["tag"]
            ngrid = s_k.shape[-1]
            nv = s_k.shape[0]

            fp.write(F"<{tag}>\n")
            fp.write(F"{nv}   {ngrid}\n")
            for k in range(ngrid):
                fp.write(F"{kgrid[k]}   ")
                for i in range(nv):
                    for j in range(i, nv):
                        fp.write(F"{s_k[i, j, k]}  ")
                fp.write("\n")

    def extend(self, rmax):

        ifft = self.ifft
        rgrid = ifft.rgrid

        dr = rgrid[1] - rgrid[0]
        nr_old = len(rgrid)
        nr_new = int(rmax / dr)

        if nr_new > nr_old:
            s_r = np.apply_along_axis(ifft.to_rspace, 2, self.s_k)
            dims = list(s_r.shape)
            dims[-1] = nr_new
            s_r_new = np.zeros(dims)
            s_r_new[..., :nr_old] = s_r[..., :nr_old]
            rgrid_new = fft_rgrid_iv(dr, nr_new)
            self.ifft = RadFFT.from_rgrid(rgrid_new)
            self.s_k = np.apply_along_axis(self.ifft.to_kspace, 2, s_r_new)

    def distance_matrix(self):

        xv = self.x
        yv = self.y
        zv = self.z

        d = (
                np.subtract.outer(xv, xv) ** 2
                + np.subtract.outer(yv, yv) ** 2
                + np.subtract.outer(zv, zv) ** 2
        )
        d = np.sqrt(d)

        return d

    def to_string(self):
        sp = StringIO()
        nv = self.nv
        print(f"Solvent:", file=sp)
        print(f"  model: {self.model}", file=sp)
        print(f"  file: {os.path.relpath(self.filename)}", file=sp)
        tbl = PrettyTable()
        tbl.set_style(PLAIN_COLUMNS)
        tbl.field_names = ["site", "sigma(Å)", "epsilon(kj/mol)", "charge"]
        for i in range(nv):
            tbl.add_row([self.aname[i], self.sigma[i], self.eps[i], self.charge[i]])
        tbl.align = "l"
        tbl.align["charge"] = "r"
        print("  geometry:", file=sp)

        for line in tbl.get_string().split("\n"):
            print(F"  {line}", file=sp)

        print("  bonds:", file=sp)
        aname = self.aname
        d = self.distance_matrix()
        for i, j in comb(range(nv), 2):
            print(F"  {aname[i]}-{aname[j]} {d[i, j]} Å", file=sp),
        print(F"  size {self.rmax:.2f} Å", file=sp)
        print(F"  kmax {self.kmax:.2f} 1/Å", file=sp)
        print(F"  reference density {self.density}", file=sp)
        print(F"  reference temp(K) {self.temp}", file=sp)
        print(F"  dielectric {self.dielectric}", file=sp)

        return sp.getvalue()

    @property
    def rmax(self):
        return self.rgrid[-1]

    @property
    def kmax(self):
        return self.kgrid[-1]

def solvent_model_locate(solvent_name):
    solvent = solvent_name + ".smdl"
    cwd = pathlib.Path.cwd()
    for path in [pathlib.Path.cwd(), DATA_DIR]:
        solvent_file = path / solvent
        if solvent_file.exists():
            break
    else:
        print(f"Cannot find {solvent=} in ", cwd, DATA_DIR)
        print("Searched in ", cwd, DATA_DIR)
        sys.exit(1)
    return solvent_file


def sik(x):
    if x < 1e-22:
        y = 1.0
    else:
        y = np.sin(x) / x
    return y


def compute_sm_rigid_bond(xv, yv, zv, kgrid):
    """
    Args:
        xv, yv, zv : x,y,z coordinates
        kgrid: reciprocal grid
    Returns:
        structure factor matrix (nv x nv x ngrid)
    """
    # form distance matrix
    d = (
            np.subtract.outer(xv, xv) ** 2
            + np.subtract.outer(yv, yv) ** 2
            + np.subtract.outer(zv, zv) ** 2
    )
    d = np.sqrt(d)

    d2 = np.multiply.outer(d, kgrid)
    sik_array = np.vectorize(sik)
    d2 = sik_array(d2)

    return d2


def extend(solvent0, rmax):


    ifft = solvent0.ifft
    rgrid = ifft.rgrid

    dr = rgrid[1] - rgrid[0]
    nr_old = len(rgrid)
    nr_new = int(rmax / dr)

    if nr_new > nr_old:

        s_r = np.apply_along_axis(ifft.to_rspace, 2, solvent0.s_k)
        dims = list(s_r.shape)
        dims[-1] = nr_new
        s_r_new = np.zeros(dims)
        s_r_new[..., :nr_old] = s_r[..., :nr_old]
        rgrid_new = fft_rgrid_iv(dr, nr_new)
        solvent = copy.deepcopy(solvent0)
        solvent.ifft = RadFFT.from_rgrid(rgrid_new)
        solvent.s_k = np.apply_along_axis(solvent.ifft.to_kspace, 2, s_r_new)
        return solvent

    return solvent0
