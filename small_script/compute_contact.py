#!/usr/bin/env python3
import os
import argparse
import sys
from time import sleep
import subprocess
from Bio.PDB import *
import numpy as np



parser = argparse.ArgumentParser(description="This is used to compute the contact")
parser.add_argument("-t", "--test", help="Test run", action="store_true", default=False)
parser.add_argument("-d", "--debug", action="store_true", default=False)
parser.add_argument("protein", help="the name of protein file")
args = parser.parse_args()

if(args.debug):
    do = print
    cd = print
else:
    do = os.system
    cd = os.chdir

# protein_name = args.protein.split('.')[0]
parser = PDBParser()

target = args.protein
structure = parser.get_structure('target', target)

# Assume only one chain
model = structure[0]
# chain = model['A']
chain = model.child_list[0]
with open("contact.csv", "w") as out:
    out.write("i, j, Distance\n")
    for i_residue in chain:
        is_regular_res = i_residue.has_id('CA')
        # print(i_residue, is_regular_res)
        if i_residue.get_resname() == "GLY":
            res_i = i_residue['CA']
        else:
            res_i = i_residue['CB']
        for j_residue in chain:
            if j_residue.get_resname() == "GLY":
                res_j = j_residue['CA']
            else:
                res_j = j_residue['CB']
            distance = res_i-res_j
            i = i_residue.id[1]
            j = j_residue.id[1]
            out_line = str(i) + ", " + str(j) + ", " + str(distance) + "\n"
            out.write(out_line)
        # print(j_residue.get_resname())
    # if id[0] != ' ' or not is_regular_res:
    #     # print(residue.id[1], -2, is_regular_res, id[0])
    #     residue.id = (' ', -1, ' ')
    #     # chain.detach_child(id)
    # # chain.detach_child(id)
    # else:
    #     subject_start = id[1]
    #     # print(subject_start)
    #     inside = False
    #     for seq in subject_seq_list:
    #         if subject_start >= seq.start and subject_start <= seq.end:
    #             residue.id = (' ', subject_start - seq.start + seq.query_start, ' ')
    #             inside = True
    #             # print(residue.id)
    #             break
    #     if not inside:
    #         # print(residue.id[1], -1)
    #         residue.id = (' ', -1, ' ')
    #         # chain.detach_child(id)
# if len(chain) == 0:
#     model.detach_child(chain.id)
#
# xyz_CA1 = []
# xyz_CA2 = []
# xyz_CA3 = []
# xyz_CA4 = []
#
# s = p.get_structure(struct_id, filename)
# chains = s[0].get_list()
# for ch in chains:
#     sequance = []
#     bonds = []
#     angles = []
#     dihedrals = []
#     if output_fn!="":
#     pass
#     else:
#         print "Chain:", ch.get_id()
#     four_res = [None, None, None, None]
#     all_res = []
#     for res in ch:
# #        is_regular_res = res.has_id('N') and res.has_id('CA') and res.has_id('C')
#     is_regular_res = res.has_id('CA') and res.has_id('O')
#
#     res_id = res.get_id()[0]
#         if (res_id==' ' or res_id=='H_MSE' or res_id=='H_M3L' or res_id=='H_CAS') and is_regular_res:
#             all_res.append(res)
#             four_res.append(res)
#             p_res = four_res.pop(0)
#             sequance.append(res.get_resname())
#             if four_res[2] and four_res[3]:
#                 if four_res[2].get_id()[1]+1!=four_res[3].get_id()[1]:
#                     print "Error: Wrong residue order"
#                 xyz_CA3 = four_res[2]['CA'].get_coord()
#                 xyz_CA4 = four_res[3]['CA'].get_coord()
#                 r = calc_bond(xyz_CA3, xyz_CA4)
#                 bonds.append(r)
#             if four_res[1] and four_res[2] and four_res[3]:
#                 if four_res[1].get_id()[1]+1!=four_res[2].get_id()[1]:
#                     print "Error: Wrong residue order"
#                 xyz_CA2 = four_res[1]['CA'].get_coord()
#                 theta = calc_angle(xyz_CA2, xyz_CA3, xyz_CA4)
#                 angles.append(theta)
#             if four_res[0] and four_res[1] and four_res[2] and four_res[3]:
#                 if four_res[0].get_id()[1]+1!=four_res[1].get_id()[1]:
#                     print "Error: Wrong residue order"
#                 xyz_CA1 = four_res[0]['CA'].get_coord()
#                 phi = calc_dihedral_angle(xyz_CA1, xyz_CA2, xyz_CA3, xyz_CA4)
#                 dihedrals.append(phi)
#
#     isNative = []
#     sigma = []
#     for i in range( 0, len(all_res) ):
#     zero = []
#     zerob = []
#     for ii in range( 0, len(all_res) ):
#         zero.append(0)
#         zerob.append(False)
#     isNative.append(zerob)
#     sigma.append(zero)
#         for j in range( i+4, len(all_res) ):
#             ires = all_res[i]
#         jres = all_res[j]
#             isNative[i][j] = checkIfNative(ires, jres)
#         if isNative[i][j]:
#                 xyz_CAi = ires['CA'].get_coord()
#                 xyz_CAj = jres['CA'].get_coord()
#         v = vector(xyz_CAi, xyz_CAj)
#                 sigma[i][j] = vabs(v)
#         else:
#                 sigma[i][j] = sigma0
#
#     if output_fn!="":
#         if splite and len(sequance)==0: continue
#         if splite:
#         if len(chains)==1:
#         file_name = output_fn+".data"
#         else:
#         file_name = output_fn+"_"+ch.get_id()+".data"
#             out = open( file_name, 'w' )
#
