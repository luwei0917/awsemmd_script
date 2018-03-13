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
import glob
from time import sleep
import fileinput
import numpy as np
from small_script.variable_test import variable_test
from small_script.variable_test2 import variable_test2
import subprocess
from small_script.myFunctions import compute_theta_for_each_helix
from small_script.myFunctions import *

# temp_list = ["350", "400", "450"]
# temp_list = ["350", "400", "450", "500", "550"]
temp_list = ["400"]
make_metadata_2(temps_list=temp_list,k=0.02)
