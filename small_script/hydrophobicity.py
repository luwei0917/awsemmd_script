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
import subprocess

table = {"A":0.5, "R":1.81, "N":0.85, "D":3.64, "B":0.43, "C":-0.02, "Q":0.77, "E":3.63, "Z":0.11,
         "G":1.15, "H":0.11, "I":-1.12, "L":-1.25, "K":2.8, "M":-0.67, "F":-1.71, "P":0.14, "S":0.46,
         "T":0.25, "W":-2.09, "Y":-0.71, "V":-0.46}

# print("A", table["A"])
with open("2xov.seq", "r") as f:
    for line in f:
        seq = line.strip()
        for i in seq:
            print(table[i])
            # print(i, table[i])
