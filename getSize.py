#!/usr/bin/env python
import os
import sys
import random
import time
from random import seed, randint
import argparse
import platform
from datetime import datetime
import imp
import numpy as np
import fileinput
from itertools import product
import pandas as pd
from scipy.interpolate import griddata
from scipy.interpolate import interp2d
import seaborn as sns
from os import listdir

import matplotlib.pyplot as plt
import seaborn as sns
from scipy.interpolate import griddata
import matplotlib as mpl
# sys.path.insert(0,'..')
# from notebookFunctions import *
# from .. import notebookFunctions

from Bio.PDB.PDBParser import PDBParser
from pyCodeLib import *



parser = argparse.ArgumentParser(
    description="This is a python3 script to\
    automatic copy the template file, \
    run simulation and analysis")

parser.add_argument("protein", help="The name of the protein")
parser.add_argument("-v", "--verbose", action="store_true", default=False)
# parser.add_argument("-l", "--label", type=str, default="/Users/weilu/opt/ff_contact/1r69/")
# parser.add_argument("-f", "--familyFold", action="store_true", default=False)
args = parser.parse_args()


parser = PDBParser(PERMISSIVE=1,QUIET=True)
structure = parser.get_structure("X", args.protein)
count = 0
for res in structure.get_residues():
    if res.get_id()[0] != " ":
        continue   # skip
    try:
        res["CA"].get_vector()
    except:
        print(structure, res.get_full_id())
        continue
    count += 1
    if args.verbose:
        print(count, res.get_full_id()[2], res.get_id()[1], res.get_id()[1]-count)
size = count
print(args.protein, size)
# print("-----")
