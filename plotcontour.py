#!/usr/bin/env python3
'''
Generates contour plot from 3 columns of data.
'''

import sys
import argparse
import os
import numpy as np
from numpy.random import uniform

#from matplotlib import rc
#rc('text', usetex=True)

import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from scipy.interpolate import interp2d

parser = argparse.ArgumentParser(description='Plots pmf data.')
parser.add_argument("filename", nargs='?', help="input filename", default="pmf-350.dat")
parser.add_argument("outname", nargs='?', help="output filename", default="test.png")
parser.add_argument("-dpi", default=150, type=int, help="figure dpi")
parser.add_argument("-x", default=1, type=int, help="x column number in f")
parser.add_argument("-xmin", default=0, type=float, help="x axis lower limit")
parser.add_argument("-xmax", default=1, type=float, help="x axis upper limit")
parser.add_argument("-y", default=2, type=int, help="y column number in f")
parser.add_argument("-ymin", default=0, type=float, help="y axis lower limit")
parser.add_argument("-ymax", default=1, type=float, help="y axis upper limit")
parser.add_argument("-z", default=3, type=int, help="z column number in f")
parser.add_argument("-zmin", default=0, type=float, help="z axis lower limit")
parser.add_argument("-zmax", default=30, type=float, help="z axis upper limit")
parser.add_argument("-title", default='', help="title")
parser.add_argument("-xlabel", default='', help="xlabel")
parser.add_argument("-ylabel", default='', help="ylabel")
parser.add_argument("-axisfontsize", default=18, type=float, help="font size of xlabel, ylabel")
parser.add_argument("-titlefontsize", default=28, type=float, help="font size of title")
args = parser.parse_args()

mpl.rcParams.update({'font.size': args.axisfontsize})

data = np.loadtxt(args.filename)
data = data[~np.isnan(data).any(axis=1)] # remove rows with nan
data = data[~(data[:,args.z] > args.zmax)] # remove rows of data for z not in [zmin zmax]
data = data[~(data[:,args.z] < args.zmin)]

xi = np.linspace(min(data[:,args.x]), max(data[:,args.x]), 50)
yi = np.linspace(min(data[:,args.y]), max(data[:,args.y]), 50)
zi = griddata((data[:,args.x], data[:,args.y]), data[:,args.z], (xi[None,:], yi[:,None]), method='linear')
#plt.contour(xi, yi, zi, 50, linewidths=0.25,colors='k')
jet = cm = plt.get_cmap('jet')
print(jet)
plt.contourf(xi, yi, zi, 50, cmap='rainbow')

plt.xlim(args.xmin, args.xmax)
plt.ylim(args.ymin, args.ymax)
plt.clim(args.zmin, args.zmax)
plt.colorbar()

plt.xlabel(args.xlabel)
plt.ylabel(args.ylabel)
plt.title(args.title, y=1.02, fontsize = args.titlefontsize)
#plt.tight_layout()
#plt.axis('equal')
#plt.axes().set_aspect('equal')
#plt.axes().set_aspect('scaled')
plt.savefig(args.outname, dpi=args.dpi, bbox_inches='tight')
os.system("open " + args.outname)
