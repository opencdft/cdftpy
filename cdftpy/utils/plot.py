#!/usr/bin/python
# -*- coding: utf-8 -*-
"""This is the module."""

import numpy as np
import matplotlib.pyplot as plt
import math


def simple_plot(r, f, xrange=None, yrange=None):
    plt.plot(r,f)
    if xrange is not None:
        plt.xlim(xrange)
    if yrange is not None:
        plt.ylim(yrange)

    # plt.axhline(0)
    plt.grid()
    plt.show()

def simple_2_plot(r, f1, f2, xrange=None, yrange=None, label=("f1","f2")):
    # plt.plot(r,f1,marker='o',markevery=50,linestyle='solid')
    plt.plot(r,f1,linestyle='solid', label=label[0])
    plt.plot(r,f2,linestyle='dashed',label=label[1])
    plt.legend()
    if xrange is not None:
        plt.xlim(xrange)
    if yrange is not None:
        plt.ylim(yrange)

    # plt.axhline(0)
    plt.grid()
    plt.show()

def simple_3_plot(r, f1, f2, f3, xrange=None, yrange=None):
    plt.plot(r,f1,marker='o',markevery=50,linestyle='solid')
    plt.plot(r,f2,marker='x',markevery=50,linestyle='solid')
    plt.plot(r,f3)
    if xrange is not None:
        plt.xlim(xrange)
    if yrange is not None:
        plt.ylim(yrange)

    plt.axhline(0)
    plt.grid()
    plt.show()


if __name__ == '__main__':

    pi = math.pi
    rmax = 3.0
    ntot = 30
    dr = rmax / ntot
    rg = np.linspace(dr, rmax, ntot, endpoint=True)
    a = 1.0
    fr = rg * np.exp(-a * rg ** 2)
    # fk =
    simple_plot(rg, fr)
