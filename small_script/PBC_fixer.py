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
import pandas as pd
import numpy as np
import glob



parser = argparse.ArgumentParser(description="This code read the dump file.\
                                output a new dump file with PBC fixed.")

parser.add_argument("-s", "--source", help="from", default="dump.unwraped")
parser.add_argument("-d", "--destination", help="to", default="new.dump")
args = parser.parse_args()

chain_num = 12
# def check(a):
#     print(a)

def check(a):
    pass

from_file = args.source
to_file = args.destination
# file = "dump.small"

with open(to_file, "w") as out:
    with open(from_file, "r") as f:
        while True:
            try:
                line = next(f)
                out.write(line)
                assert line.strip() == "ITEM: TIMESTEP"
            except:
                break
            line = next(f)
            out.write(line)
            timesteps = int(line)
            if timesteps % 1e4 == 0:
                print(timesteps)
            check(timesteps)

            line = next(f)
            out.write(line)
            assert line.strip() == "ITEM: NUMBER OF ATOMS"
            line = next(f)
            num = int(line.strip())
            check(num)
            out.write(line)

            atom_in_one_chain = int(num/chain_num)
            check(atom_in_one_chain)
    #         assert next(f).strip() == "ITEM: BOX BOUNDS ss ss ss"
            line = next(f)
            out.write(line)

            line = next(f)
            out.write(line)
            xStart, xEnd = line.split()
            xStart, xEnd = float(xStart), float(xEnd)
            line = next(f)
            out.write(line)
            yStart, yEnd = line.split()
            yStart, yEnd = float(yStart), float(yEnd)
            line = next(f)
            out.write(line)
            zStart, zEnd = line.split()
            zStart, zEnd = float(zStart), float(zEnd)
    #         assert next(f).strip() == "ITEM: ATOMS id type xs ys zs"
            line = next(f)
            out.write(line)
            for chain in range(chain_num):
                data = []
                for i in range(atom_in_one_chain):
                    n, atomType, x, y, z = next(f).split()
                    n = int(n)
                    atomType = int(atomType)
                    x,y,z = float(x), float(y), float(z)
                    tmp = [n, atomType, x, y, z]
                    data.append(tmp)
                data = pd.DataFrame(data, columns=["n", "atomType", "x", "y", "z"])

                for direction in ["x", "y", "z"]:
                    if (data[direction] > 1.1).all():
                        data[direction] = data[direction] - 1
                    if (data[direction] < -0.1).all():
                        data[direction] = data[direction] + 1
                for i, row in data.iterrows():
                    out.write(f'{int(row["n"])} {int(row["atomType"])} {row["x"]} {row["y"]} {row["z"]} \n')
