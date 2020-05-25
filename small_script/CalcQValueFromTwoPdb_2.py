#!/usr/bin/python

# ----------------------------------------------------------------------
# Copyright (2010) Aram Davtyan and Garegin Papoian

# Papoian's Group, University of Maryland at Collage Park
# http://papoian.chem.umd.edu/

# Last Update: 03/04/2011
# ----------------------------------------------------------------------

import sys
from VectorAlgebra import *

#from Bio.PDB.PDBParser import PDBParser


if len(sys.argv)!=4 and len(sys.argv)!=3:
    print("\nCalcQValue.py PDB_Id pdb2  [sigma_exp] [-i]\n")
    # print
    print("\t\t-i\tcalculate individual q values for each chain")
    # print
    exit()

splitq = False
for iarg in range(0, len(sys.argv)):
    if sys.argv[iarg]=="-i":
        splitq = True
        sys.argv.pop(iarg)

struct_id = sys.argv[1]
if struct_id[-4:].lower()==".pdb":
    pdb_file = struct_id
else:
    pdb_file = struct_id + ".pdb"
struct_id2 = sys.argv[2]
if struct_id2[-4:].lower()==".pdb":
    pdb_file2 = struct_id2
else:
    pdb_file2 = struct_id2 + ".pdb"


sigma_exp = 0.15
if len(sys.argv)==4:
    sigma_exp = float(sys.argv[3])

n_atoms = 0
i_atom = 0
item = ''
step = 0
ca_atoms_pdb = {}
pdb_chain_id = []
pdb_residue_id = {}
ca_atoms_pdb_2 = {}
pdb_chain_id_2 = []
pdb_residue_id_2 = {}
box = []
A = []
sigma = []
sigma_sq = []




from Bio.PDB.PDBParser import PDBParser

p = PDBParser(PERMISSIVE=1, QUIET=True)

def computeQ(ca_atoms_pdb, ca_atoms_pdb_2):
    if len(ca_atoms_pdb_2)!=len(ca_atoms_pdb):
        print("Notice. Length mismatch!")
        print("Pdb: ", len(ca_atoms_pdb), "trj: ", len(ca_atoms_pdb_2))
        # exit()
    # print("hi")
    max_length = max(len(ca_atoms_pdb_2), len(ca_atoms_pdb))
    # print(max_length)
    # ca_atoms = ca_atoms[:max_length]
    # ca_atoms_pdb = ca_atoms_pdb[:max_length]
    # print ca_atoms_pdb
    # print ca_atoms_pdb_2
    Q = 0
    norm = 0
    N = max_length
    for ia in range(0, N):
        for ja in range(ia+3, N):
            if (ia+1) in ca_atoms_pdb and (ja+1) in ca_atoms_pdb:
                if (ia+1) in ca_atoms_pdb_2 and (ja+1) in ca_atoms_pdb_2:
                        # print ia-pdbBegin+1, ja-pdbBegin+1
                    # print ca_atoms_pdb[ia+1]
                    rn = vabs(vector(ca_atoms_pdb[ia+1], ca_atoms_pdb[ja+1]))
                    r = vabs(vector(ca_atoms_pdb_2[ia+1], ca_atoms_pdb_2[ja+1]))
                    dr = r - rn

                    Q = Q + exp(-dr * dr / (2 * sigma_sq[ja - ia]))
                    norm = norm + 1
    Q = Q / norm
    return Q

s = p.get_structure(struct_id, pdb_file)
chains = s[0].get_list()
#chain = chains[0]
ichain = 0
for chain in chains:
    ichain = ichain + 1
    for res in chain:
        is_regular_res = res.has_id('CA') and res.has_id('O')
        res_id = res.get_id()[0]
        if (res_id==' ' or res_id=='H_MSE' or res_id=='H_M3L' or res_id=='H_CAS') and is_regular_res:
            residue_id = res.id[1]
            ca_atoms_pdb[residue_id] = res['CA'].get_coord()
            pdb_chain_id.append(ichain)
            pdb_residue_id[res.id[1]] = 1






s = p.get_structure(struct_id2, pdb_file2)
chains = s[0].get_list()
#chain = chains[0]
ichain = 0
for chain in chains:
    ichain = ichain + 1
    for res in chain:
        is_regular_res = res.has_id('CA') and res.has_id('O')
        res_id = res.get_id()[0]
        if (res_id==' ' or res_id=='H_MSE' or res_id=='H_M3L' or res_id=='H_CAS' ) and is_regular_res:
            # ca_atoms.append(res['CA'].get_coord())
            residue_id_2 = res.id[1]
            if residue_id_2 in pdb_residue_id:
                # print(residue_id_2)
                ca_atoms_pdb_2[residue_id_2] = res['CA'].get_coord()
                pdb_chain_id_2.append(ichain)
                pdb_residue_id_2[res.id[1]] = 1
            else:
                print("no: " + str(residue_id_2))
# print("------")
n = max(len(ca_atoms_pdb_2), len(ca_atoms_pdb))
for i in range(0, n+1):
    sigma.append( (1+i)**sigma_exp )
    sigma_sq.append(sigma[-1]*sigma[-1])

if len(ca_atoms_pdb_2)>0:
    q = computeQ(ca_atoms_pdb, ca_atoms_pdb_2)
    # print str(round(q,3)),
    sys.stdout.write(str(round(q,3)) + " ")
    n_atoms = len(ca_atoms_pdb_2)
