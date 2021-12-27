#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Provides various conversion units according to 2018 CODATA recommended values
see https://physics.nist.gov/cuu/Constants/index.html
"""

# Hartree energy 4.3597447222071e-21 kJ
from decimal import Decimal

hartree = 4.3597447222071e-21
# Avogadro number 6.02214076e+23 mol-1
nav = 6.02214076e+23
# Molar gas constant 0.008314462618 kJ/(mol*K)
# R = 0.008314462618
R = 0.00831446  # Gennady's value

# Bohr radius 0.529177210903 Å
bohr = 0.529177210903
# The units used in the code are
#   energy kj/mol
#   length Å
#   charge e (elementary charge)

# hartree_2_kjmol = hartree*nav
hartree_2_kjmol = 2625.4956964241155 # gennady's value

bohr_2_ang = 1.0/bohr
# epot_2_kjmol = hartree_2_kjmol/bohr_2_ang
# epot_2_kjmol = 1389.352489871543
epot_2_kjmol = 2625.4956964241155
# print(epot_2_kjmol,hartree_2_kjmol)

# OLD VALUES (to be removed)
# N_AVO = 6.02214179e+23
# E_CHARGE = 4.803e-10
# E_CHARGE_2 = 23.068809e-20
# # CK = N_AVO * E_CHARGE ** 2 * 0.01
# CK = 60.2214179*23.068809

