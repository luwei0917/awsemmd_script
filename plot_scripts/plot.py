#!/usr/bin/env python3
import os
os.system("echo hello world")
os.system("sort -r 4cpv_HA/analysis/list_of_max_q > ha_q")
os.system("sort -r 4cpv_time_step_3/analysis/list_of_max_q > he_q")
os.system("gnuplot q_values.plt")
