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
from time import sleep

# from run_parameter import *
parser = argparse.ArgumentParser(
    description="This is a python3 script to\
    automatic copy the template file, \
    run simulation and analysis")

parser.add_argument("template", help="the name of template file")
args = parser.parse_args()
protein_name = args.template.strip('./')
energy_cmd = r"""
variable    step equal step

# f_4 = Total energy computed in the fix backbone
# f_4[i] = {eChain eShake eChi eRama eExcluded eDSSP eP_AP eWater eBurial eHelix eAMH-Go eFrag_Mem eSSB}
variable E_chain equal f_4[1]
variable E_shake equal f_4[2]
variable E_chi equal f_4[3]
variable E_rama equal f_4[4]
variable E_excluded equal f_4[5]
variable E_dssp equal f_4[6]
variable E_pap equal f_4[7]
variable E_water equal f_4[8]
variable E_burial equal f_4[9]
variable E_helix equal f_4[10]
variable E_amhgo equal f_4[11]
variable E_fmem equal f_4[12]
variable E_ssb equal f_4[13]
variable E_vecFM equal f_4[15]
variable E_membrane equal f_4[14]
variable E_bond equal emol
variable E_excl equal epair

variable E_P equal v_E_chain+v_E_shake+v_E_chi+v_E_rama+v_E_excluded+v_E_dssp+v_E_pap+v_E_water+v_E_burial+v_E_helix+v_E_amhgo+v_E_fmem+v_E_vecFM+v_E_membrane+v_E_ssb+v_E_bond+v_E_excl
variable E_P2 equal v_E_chain+v_E_shake+v_E_chi+v_E_rama+v_E_excluded+v_E_dssp+v_E_pap+v_E_water+v_E_burial+v_E_helix+v_E_amhgo+v_E_membrane+v_E_ssb+v_E_bond+v_E_excl
variable E_K     equal ke
variable E_total equal v_E_P+v_E_K



fix energy all print 1000 "${step} ${E_chain} ${E_shake} ${E_chi} ${E_rama} ${E_excluded} ${E_dssp} ${E_pap} ${E_water} ${E_burial} ${E_helix} ${E_amhgo} ${E_fmem} ${E_vecFM} ${E_membrane} ${E_ssb}  ${E_P} ${E_bond} ${E_excl}" file energy.dat screen no title "Step    Chain           Shake           Chi             Rama            Excluded        DSSP            P_AP            Water           Burial          Helix           AMH-Go       Frag_Mem        Vec_FM          Membrane        SSB          VTotal      Ebond     Epair"

fix         	wham all print 1000 "${step}  ${E_P}" file wham.dat screen no title "# timestep   energy"
"""
with open("energy.in", "w") as out:
    with open(protein_name+".in") as f:
        for line in f:
            # print(line.split())
            data = line.split()
            if(len(data) > 0):
                if(data[0] == "run"):
                    out.write("rerun dump.lammpstrj start 1000 dump x y z\n")
                elif(data[0] == "velocity"):
                    pass
                elif(data[0] == "dump"):
                    pass
                elif(data[0] == "dump_modify"):
                    pass
                elif(data[0] == "reset_timestep"):
                    pass
                elif(data[0] == "restart"):
                    pass
                else:
                    out.write(line)
            else:
                out.write(line)
            match = "# calculate and output collective variables during the simulation #\n"
            if line == match:
                out.write(energy_cmd)
