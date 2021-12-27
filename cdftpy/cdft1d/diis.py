#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
This module provides collection
of routines related to diis procedure
"""
import logging
from collections import deque


import numpy as np

logging.basicConfig(level=logging.DEBUG)

_initialized = False


class Diis:

    def __init__(self):
        self.f = deque([])
        self.df = deque([])

    def size(self):
        return len(self.f)

    def accumulate(self, g, dg):
        self.f.append(g)
        self.df.append(dg)

    def solve(self):
        nd = self.size()
        ngrid = self.f[0].shape[-1]
        a = -np.ones([nd + 1, nd + 1])
        a[nd, nd] = 0
        for i in range(nd):
            for j in range(i, nd):
                dij = self.df[i] * self.df[j]
                a[i][j] = np.sum(dij) / ngrid
                a[j][i] = a[i][j]

        b = np.zeros([nd + 1])
        b[nd] = -1
        x = np.linalg.solve(a, b)

        g = np.zeros(self.f[0].shape)
        for k in range(nd):
            g = g + x[k] * self.f[k]

        self.f.popleft()
        self.df.popleft()
        return g


def diis_session():
    f = deque([])
    df = deque([])

    def diis(nd, g, dg):
        f.append(g)
        df.append(dg)
        if len(f) == nd:
            ngrid = f[0].shape[-1]
            a = -np.ones([nd + 1, nd + 1])
            a[nd, nd] = 0
            for i in range(nd):
                for j in range(i, nd):
                    dij = df[i] * df[j]
                    a[i][j] = np.sum(dij) / ngrid
                    a[j][i] = a[i][j]
            b = np.zeros([nd + 1])
            b[nd] = -1
            x = np.linalg.solve(a, b)
            g = np.zeros(f[0].shape)
            for i in range(nd):
                g = g + x[i] * f[i]

            f.popleft()
            df.popleft()

        return g

    return diis


