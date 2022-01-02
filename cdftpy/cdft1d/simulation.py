import inspect
import sys
import traceback
from io import StringIO
from types import SimpleNamespace

import numpy as np

from cdftpy.cdft1d.coulomb import compute_long_range_coul_pot_kspace, compute_long_range_coul_pot_rspace, \
    compute_short_range_coul_pot_rspace, compute_coulomb_potential
from cdftpy.cdft1d.exceptions import ConvergenceError
from cdftpy.cdft1d.potential import compute_lj_potential_mod, compute_lj_potential
from cdftpy.utils.units import kb
import cdftpy.cdft1d.rism as rism
import cdftpy.cdft1d.rsdft as rsdft
from cdftpy.cdft1d.solvent import Solvent, solvent_model_locate
from cdftpy import __version__

DEFAULT_PRINT_LEVEL=frozenset(['parameters','solvent','solute','header'])
DEFAULT_PARAMS = dict(diis_iterations=2, tol=1.0e-9, output_rate=10, max_iter=200, rcoul=1.25, method="rsdft")

HEADER = F"""
==================================
1D CDFT PROGRAM

version {__version__}
Marat Valiev and Gennady Chuev
==================================
"""

class SolvatedIon():

    def __init__(self,
                 method="rsdft",
                 params=None,
                 solute=None,
                 solvent=None):

        self._solvent = None
        self.method = method

        if params is None:
            params = DEFAULT_PARAMS
        params = {**DEFAULT_PARAMS, **params}

        self.ndiis = int(params["diis_iterations"])
        self.rcoul = float(params["rcoul"])
        self.tol = float(params["tol"])
        self.max_iter = int(params["max_iter"])
        self.output_rate = int(params["output_rate"])
        if "temp" in params:
            self.temp = float(params["temp"])
            self.beta = 1.0 / (kb * temp)
        else:
            self.temp = None
            self.beta = None

        if "rmax" in params:
            self.rmax = float(params["rmax"])
        else:
            self.rmax = None

        self.solute = solute

        if solvent is not None:
            if isinstance(solvent, str):
                rism_patch = (method == "rism")
                self.solvent = Solvent.from_name(solvent, rism_patch=rism_patch)
            else:
                self.solvent = solvent

        self.log = None
        self.h_r = None
    @property
    def solvent(self):
        return self._solvent

    @solvent.setter
    def solvent(self, slv):
        if slv.nv > 2 and self.method == "rsdft":
            print(F"{self.method.upper()} is currently restricted to diatomic liquids")
            print("As an alternative you can try RISM")
            sys.exit(1)
        if self.temp is None:
            self.temp = slv.temp
            self.beta = 1.0 / (kb * self.temp)

        if self.rmax is not None:
            slv.extend(self.rmax)
        self.rmax = slv.ifft.rgrid[-1]

        self.hbar = slv.hbar(self.method)
        self.sm = slv.sm_k(slv.kgrid)
        self._solvent = slv

    @property
    def vl_k(self):
        vl_k = compute_long_range_coul_pot_kspace(self.solute["charge"], self.qv,
                                                  self.kgrid, r_s=self.rcoul)
        return self.beta * vl_k

    @property
    def vl_r(self):
        vl_r = compute_long_range_coul_pot_rspace(self.solute["charge"], self.qv,
                                                  self.rgrid, r_s=self.rcoul)
        return self.beta * vl_r

    @property
    def vs_r(self):
        vcs_r = compute_short_range_coul_pot_rspace(self.solute["charge"], self.qv,
                                                    self.rgrid, r_s=self.rcoul)
        if self.method == "rsdft":
            vlj_r = compute_lj_potential_mod(self.solute['sigma'], self.solute['eps'], self.sig_v, self.eps_v, self.rgrid)
        else:
            vlj_r = compute_lj_potential(self.solute['sigma'], self.solute['eps'], self.sig_v, self.eps_v, self.rgrid)
        vs_r = vcs_r + vlj_r

        return self.beta * vs_r


    @property
    def ifft(self):
        return self.solvent.ifft

    @property
    def rgrid(self):
        return self.solvent.rgrid

    @property
    def kgrid(self):
        return self.solvent.kgrid

    @property
    def s_k(self):
        return self.solvent.s_k

    @property
    def zeta(self):
        return self.solvent.zeta(self.beta)

    @property
    def rho_0(self):
        return self.solvent.density

    @property
    def sig_v(self):
        return self.solvent.sigma

    @property
    def eps_v(self):
        return self.solvent.eps

    @property
    def qv(self):
        return self.solvent.charge

    def to_string(self, print_level=None):

        sp = StringIO()
        if print_level is None:
            print_level = DEFAULT_PRINT_LEVEL
        if 'header' in print_level:
            print(HEADER)
        if 'solute' in print_level:
            print(F"Solute:  {self.solute['name']} charge={self.solute['charge']} "
                  F"sigma={self.solute['sigma']} "
                  F"epsilon={self.solute['eps']}", file=sp)
        if 'parameters' in print_level:
            buffer = F"Parameters: \n" \
                     F"  method: {self.method}\n" \
                     F"  coulomb cutoff: {self.rcoul} Å\n" \
                     F"  gamma tolerance: {self.tol} \n" \
                     F"  max iter: {self.max_iter} "
            print(buffer, file=sp)
        if 'solvent' in print_level:
            if self.solvent is None:
                print("Solvent:None", file=sp)
            else:
                print(self.solvent.to_string(), file=sp)
        return sp.getvalue()

    def cdft(self, quiet=True):
        try:
            if self.method == "rism":
                fe = rism.rism_1d(self,quiet=quiet)
                return fe
            elif self.method == "rsdft":
                fe = rsdft.rsdft_1d(self,quiet=quiet)
                return fe
        except ConvergenceError as e:
            print(e)
            sys.exit(1)
        except (AttributeError, TypeError) as e:
            traceback.print_exc()
            print("Simulation has not been properly initialized")
            self.integrity_check()
            sys.exit(1)
        except BaseException as e:
            raise e
    @property
    def epot_r(self):
        ifft = self.solvent.ifft
        h_k = np.apply_along_axis(ifft.to_kspace, 1, self.h_r)
        epot_r, epot_k = compute_coulomb_potential(self.rho_0,
                                                   self.solvent.charge, ifft, h_k)
        return epot_r

    def integrity_check(self):
        pre_message = "It seems that"
        if self.solvent is None:
            print(f"{pre_message} solvent has not been defined")
        if self.solute is None:
            print(f"{pre_message} solute has not been defined")

if __name__ == '__main__':
    solvent_name = "s2"
    filename = solvent_model_locate(solvent_name)
    solvent = Solvent.from_file(filename)
    solute = dict(name="Cl", charge=-1.0, sigma=4.83, eps=0.05349244)
    # solute = dict(name="Cl", charge=-1.0, sigma=60, eps=0.05349244)
    params = dict(diis_iterations=2, tol=1.0e-7, max_iter=500)

    sim = SolvatedIon(solvent=solvent, params=params, method="rism")

    # sim.solvent = solvent
    # sim = SolvatedIon(solute=solute, solvent="s2", params=params, method="rism")
    fe=sim.cdft()
    print(sim.to_string())
    # print(F"Free energy of solvation {fe:.3f} kj/mol")
    #
