#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Provides various conversion units according to 2018 CODATA recommended values
see https://physics.nist.gov/cuu/Constants/index.html
"""
# The units used in the code are
#   energy kj/mol
#   length Å
#   charge e (elementary charge)

# Hartree energy 4.3597447222071e-21 kJ
hartree = 4.3597447222071e-21

# Avogadro number 6.02214076e+23 mol-1
nav = 6.02214076e+23

# Molar gas constant 0.008314462618 kJ/(mol*K)
kb = 0.008314462618

# Bohr radius 0.529177210903 Å
bohr = 0.529177210903
bohr_2_ang = 1.0/bohr

# hartree_2_kjmol = hartree*nav
hartree_2_kjmol = 2625.4996394798254 # gennady's value




