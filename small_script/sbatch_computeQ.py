#!/usr/bin/env python3
import os
import argparse

parser = argparse.ArgumentParser(description="This is used to compute the contact")
parser.add_argument("-i", help="Test run", type=int, default=-1)
parser.add_argument("-t", type=int, default=-1)
parser.add_argument("-n", type=int, default=-1)
args = parser.parse_args()


def computeMutualQ(i, thousand, n):
    try:
        os.remove(f"rm line_{i}_{thousand}")
    except OSError:
        pass
    for jj in range(1000):
        j = thousand*1000 + jj
        if j >= n:
            break
        os.system(f"python2 ~/opt/small_script/CalcQValueFromTwoPdb_2.py lowTstructure{i}.pdb lowTstructure{j}.pdb | tail -n 1 >> line_{i}_{thousand}")

computeMutualQ(args.i, args.t, args.n)
