#!/usr/bin/env python3
import subprocess
import os
import sys
# my_env = os.environ.copy()
# my_env["PATH"] = "/home/wl45/python/bin:" + my_env["PATH"]
# print(my_env)
# subprocess.Popen("vmd", env=my_env)
# os.system("vmd")
# subprocess.run(["ls", "-l"])
# subprocess.run(["vmd", "-l"], shell=True)
# hello world
config = open('config.py', 'w')
n = 1
simulation_steps = 20
warm_up_steps = 20
config.write("number_of_run = %d\nsimulation_steps = %d\n\
warm_up_steps = %d\n" % (n, simulation_steps, warm_up_steps))
config.close()
sys.exit(0)
