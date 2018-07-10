#!/usr/bin/env python3
import os
import argparse
import sys
from time import sleep
import subprocess
import fileinput
from myFunctions import duplicate_pdb
parser = argparse.ArgumentParser(
    description="The goal of this python3 code is to automatically create \
    the project template as fast as possible. Written by Wei Lu."
)
parser.add_argument("From", help="Duplicate from")
parser.add_argument("To", help="Duplicate to")

args = parser.parse_args()


duplicate_pdb(args.From, args.To, offset_x=40.0, new_chain="B")