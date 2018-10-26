#!/usr/bin/env python3
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
import scipy
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.interpolate import griddata
import matplotlib as mpl
# sys.path.insert(0,'..')
# from notebookFunctions import *
# from .. import notebookFunctions

parser = argparse.ArgumentParser(description="plot the energy file")
# parser.add_argument("protein", help="the name of protein")
parser.add_argument("file", help="the name of data file")
args = parser.parse_args()

data = pd.read_csv(args.file)
data.columns = data.columns.str.strip()
x = "Steps"
y = "Membrane"
data.plot(x, y)
plt.show()