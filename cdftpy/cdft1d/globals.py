#!/usr/bin/python
# -*- coding: utf-8 -*-
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
PROJ_DIR = THIS_DIR.parent
PROJ_DIR_PARENT = PROJ_DIR.parent

DATA_DIR = PROJ_DIR / "data"

