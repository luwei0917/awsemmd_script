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
from myPersonalFunctions import *
import glob
import numpy


parser = argparse.ArgumentParser(description="This code read the dump.lammpstrj file.\
                                output imformation as requested")

parser.add_argument("-t", "--test", help="test ", action="store_true", default=False)
parser.add_argument("-d", "--debug", action="store_true", default=False)
args = parser.parse_args()


if(args.debug):
    do = print
    cd = print

    def check(a):
        print(a)
else:
    do = os.system
    cd = os.chdir

    def check(a):
        pass

file_name = "dump.lammpstrj"
# file_name = "test"
with open("distance", "w") as out:
    with open(file_name) as f:
        while True:
            try:
                assert next(f).strip() == "ITEM: TIMESTEP"
            except:
                break
            timesteps = int(next(f))
            check(timesteps)
            assert next(f).strip() == "ITEM: NUMBER OF ATOMS"
            num = int(next(f).strip())
            check(num)
            assert next(f).strip() == "ITEM: BOX BOUNDS ss ss ss"
            xStart, xEnd = next(f).split()
            xStart, xEnd = float(xStart), float(xEnd)
            yStart, yEnd = next(f).split()
            yStart, yEnd = float(yStart), float(yEnd)
            zStart, zEnd = next(f).split()
            zStart, zEnd = float(zStart), float(zEnd)
            assert next(f).strip() == "ITEM: ATOMS id type xs ys zs"
            xi = yi = zi = xj = yj = zj = 0
            for i in range(num):
                n, _, x, y, z = next(f).split()
                n = int(n)

                if n == 1:
                    xi,yi,zi = float(x), float(y), float(z)
                if n == num-2:
                    xj,yj,zj = float(x), float(y), float(z)
            dis = (((xi-xj)*(xEnd-xStart))**2 + ((yi-yj)*(yEnd-yStart))**2 + ((zi-zj)*(zEnd-zStart))**2)**0.5
            out.write(str(dis)+"\n")
            # print(xi, yi, zi, xj, yj, zj)













print("Hello World")
