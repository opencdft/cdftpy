import os
import sys
import traceback
from io import StringIO

import numpy as np

from cdftpy.cdft1d.coulomb import compute_long_range_coul_pot_kspace, compute_long_range_coul_pot_rspace, \
    compute_short_range_coul_pot_rspace, compute_coulomb_potential
from cdftpy.cdft1d.exceptions import ConvergenceError
from cdftpy.cdft1d.io_utils import read_solute, read_key_value
from cdftpy.cdft1d.potential import compute_lj_potential_mod, compute_lj_potential
from cdftpy.cdft1d.units import kb
import cdftpy.cdft1d.rism as rism
import cdftpy.cdft1d.rsdft as rsdft
from cdftpy.cdft1d.solvent import Solvent
from cdftpy.cdft1d._version import __version__

DEFAULT_PRINT_LEVEL = frozenset(['parameters', 'solvent', 'solute', 'header'])
DEFAULT_PARAMS = dict(ndiis=2, tol=1.0e-7, output_rate=10, max_iter=200, rcoul=1.25,
                      rmax=None, solvent='s2', method='rsdft')
RUNNERS = dict(rsdft=rsdft.rsdft_1d, rism=rism.rism_1d)
HEADER = F"""
==================================
1D CDFT PROGRAM

version {__version__}
Marat Valiev and Gennady Chuev
==================================
"""


class IonSolvation():

    def __init__(self,
                 name="ion",
                 charge=None,
                 sigma=None,
                 eps=None,
                 **params):

        self.name = name
        if charge is None:
            print("Charge has to be provided")
            raise ValueError
        else:
            self.charge = charge
        if sigma is None:
            print("Sigma has to be provided")
            raise ValueError
        else:
            self.sigma = sigma
        if eps is None:
            print("Sigma has to be provided")
            raise ValueError
        else:
            self.eps = eps

        self._solvent = None

        if params == {}:
            params = DEFAULT_PARAMS
        params = {**DEFAULT_PARAMS, **params}

        self.method = params['method']
        self.ndiis = int(params["ndiis"])
        self.rcoul = float(params["rcoul"])
        self.tol = float(params["tol"])
        self.max_iter = int(params["max_iter"])
        self.output_rate = int(params["output_rate"])
        self.rmax = params['rmax']

        if self.rmax is not None:
            self.rmax = float(self.rmax)

        self.temp = None
        self.beta = None

        solvent = params['solvent']
        if solvent is not None:
            if isinstance(solvent, str):
                rism_patch = (self.method == "rism")
                self.solvent = Solvent.from_name(solvent, rism_patch=rism_patch)
            else:
                self.solvent = solvent

        self.log = None
        self.h_r = None

        self.bridge = False
    @classmethod
    def from_input_file(cls, input_file, solvent=None, method=None):

        try:
            solute = read_solute(input_file)
        except FileNotFoundError:
            print(f"Cannot locate input file {input_file}")
            sys.exit(1)

        for k, v in solute.items():
            solute[k] = v[0]

        parameters = read_key_value(input_file, section="simulation")
        if parameters is None:
            parameters = DEFAULT_PARAMS
        parameters = {**DEFAULT_PARAMS, **parameters}

        if solvent is not None:
            parameters['solvent'] = solvent
        if method is not None:
            parameters['method'] = method

        return cls(**solute, **parameters)

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
        vl_k = compute_long_range_coul_pot_kspace(self.charge, self.qv,
                                                  self.kgrid, r_s=self.rcoul)
        return self.beta * vl_k

    @property
    def vl_r(self):
        vl_r = compute_long_range_coul_pot_rspace(self.charge, self.qv,
                                                  self.rgrid, r_s=self.rcoul)
        return self.beta * vl_r

    @property
    def vs_r(self):
        vcs_r = compute_short_range_coul_pot_rspace(self.charge, self.qv,
                                                    self.rgrid, r_s=self.rcoul)
        if self.bridge:
            vlj_r = compute_lj_potential_mod(self.sigma, self.eps, self.sig_v, self.eps_v, self.rgrid)
        else:
            vlj_r = compute_lj_potential(self.sigma, self.eps, self.sig_v, self.eps_v, self.rgrid)
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
            print(F"Solute:  {self.name} charge={self.charge} "
                  F"sigma={self.sigma} "
                  F"epsilon={self.eps}", file=sp)
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

        fe = None
        try:
            fe = RUNNERS[self.method](self, quiet=quiet)
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
        return fe

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

    def write_density(self,postfix="",dirpath="./"):

        filename = os.path.join(dirpath, f"density{postfix}.dat")
        with open(filename, "w") as fp:
            fp.write("# site density data\n")
            fp.write(f"# {self.name} charge={self.charge} "
                     F"sigma={self.sigma}  "
                     F"epsilon={self.eps}  "
                     F"solvent={self.solvent.model}\n")
            fp.write(f"{'# rgrid':<10} ")
            for tag in self.solvent.aname:
                fp.write(f"{self.name}-{tag}   ")
            fp.write("\n")
            for j in range(len(self.rgrid)):
                fp.write(f"{self.rgrid[j]}  ")
                for i in range(len(self.solvent.aname)):
                    fp.write(f"{self.h_r[i,j]+1.0}  ")
                fp.write("\n")

