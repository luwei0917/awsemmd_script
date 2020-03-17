#!/usr/bin/python2
import sys
import os
from VectorAlgebra import *
from Bio.PDB.PDBParser import PDBParser

atom_type = {'1': 'C', '2': 'N', '3': 'O', '4': 'C', '5': 'H', '6': 'C'}
atom_desc = {'1': 'C-Alpha', '2': 'N', '3': 'O',
             '4': 'C-Beta', '5': 'H-Beta', '6': 'C-Prime'}
PDB_type = {'1': 'CA', '2': 'N', '3': 'O', '4': 'CB', '5': 'HB', '6': 'C'}


# if len(sys.argv) != 6 and len(sys.argv) != 5:
#     print "\nCalcQValue.py PDB_Id Input_file Output_file qonuchic_flag(1 for q_o, 0 for q_w) [-i]\n"
#     print
#     print "\t\t-i\tcalculate individual q values for each chain"
#     print
#     sys.exit()
if len(sys.argv) != 6 and len(sys.argv) != 5:
    print("\nCalcQValue.py PDB_Id Input_file Output_file qonuchic_flag(1 for q_o, 0 for q_w) [-i]\n")
    print("\t\t-i\tcalculate individual q values for each chain")

    sys.exit()
cutoff = 9.5
splitq = False
for iarg in range(0, len(sys.argv)):
    if sys.argv[iarg] == "-i":
        splitq = True
        sys.argv.pop(iarg)

struct_id = sys.argv[1]
if struct_id[-4:].lower() == ".pdb":
    pdb_file = struct_id
else:
    pdb_file = struct_id + ".pdb"
lammps_file = sys.argv[2]

output_file = ""
if len(sys.argv) > 3:
    output_file = sys.argv[3]

sigma_exp = 0.15
qo_flag = int(sys.argv[4])

n_atoms = 0
i_atom = 0
item = ''
step = 0
ca_atoms_pdb = {}
pdb_chain_id = []
pdb_residue_id = {}
ca_atoms = []
box = []
A = []
sigma = []
sigma_sq = []


out = open(output_file, 'w')


p = PDBParser(PERMISSIVE=1)

def computeQ():
    # print ca_atoms_pdb
    # if len(ca_atoms)!=len(ca_atoms_pdb):
        # print "Notice."
    # pdbBegin = int(
    #     os.popen("grep \"REMARK 220 TARGET BEGIN:\" " + pdb_file).read().split(":")[1])
    # pdbEnd = int(os.popen("grep \"REMARK 220 TARGET END:\" " +
    #                       pdb_file).read().split(":")[1])
    # print "Pdb: ", len(ca_atoms_pdb), "trj: ", len(ca_atoms), "pdb Begin: ", pdbBegin, "pdb End:", pdbEnd

    Q = 0
    norm = 0
    N = len(ca_atoms)
    min_sep = 3
    if qo_flag == 1:
        min_sep = 4
    for ia in range(N):
        for ja in range(ia + min_sep, N):
            if (ia+1) in ca_atoms_pdb and (ja+1) in ca_atoms_pdb:
                    # print ia-pdbBegin+1, ja-pdbBegin+1
                # print ca_atoms_pdb[ia+1]
                rn = vabs(vector(ca_atoms_pdb[ia+1], ca_atoms_pdb[ja+1]))
                if qo_flag == 1 and rn >= cutoff:
                    continue
                r = vabs(vector(ca_atoms[ia], ca_atoms[ja]))
                dr = r - rn
                Q = Q + exp(-dr * dr / (2 * sigma_sq[ja - ia]))
                norm = norm + 1

    Q = Q / norm
    return Q

s = p.get_structure(struct_id, pdb_file)
chains = s[0].get_list()
# chain = chains[0]
ichain = 0
for chain in chains:
    ichain = ichain + 1
    for res in chain:
        is_regular_res = res.has_id('CA') and res.has_id('O')
        res_id = res.get_id()[0]
        if (res_id == ' ' or res_id == 'H_MSE' or res_id == 'H_M3L' or res_id == 'H_CAS') and is_regular_res:
            residue_id = res.id[1]
            ca_atoms_pdb[residue_id] = res['CA'].get_coord()
            # ca_atoms_pdb.append(res['CA'].get_coord())
            pdb_chain_id.append(ichain)
            pdb_residue_id[res.id[1]] = 1
            # print res.id[1]
for i in range(0, len(ca_atoms_pdb) + 1):
    sigma.append((1 + i)**sigma_exp)
    sigma_sq.append(sigma[-1] * sigma[-1])

lfile = open(lammps_file)
for l in lfile:
    l = l.strip()
    if l[:5] == "ITEM:":
        item = l[6:]
    else:
        if item == "TIMESTEP":
            if len(ca_atoms) > 0:
                q = computeQ()
                out.write(str(round(q, 3)))
                out.write(' ')
                out.write('\n')
                n_atoms = len(ca_atoms)
            step = int(l)
            ca_atoms = []
            box = []
            A = []
        elif item == "NUMBER OF ATOMS":
            n_atoms = int(l)
        elif item[:10] == "BOX BOUNDS":
            box.append(l)
            l = l.split()
            A.append([float(l[0]), float(l[1])])
        elif item[:5] == "ATOMS":
            l = l.split()
            i_atom = l[0]
            x = float(l[2])
            y = float(l[3])
            z = float(l[4])
            x = (A[0][1] - A[0][0]) * x + A[0][0]
            y = (A[1][1] - A[1][0]) * y + A[1][0]
            z = (A[2][1] - A[2][0]) * z + A[2][0]
            desc = atom_desc[l[1]]
            if desc == 'C-Alpha':
                #                atom = Atom(i_atom, atom_type[l[1]], l[1], x, y, z, desc)
                residue_id = int((int(i_atom) + 2) / 3)
                # print residue_id
                atom = [x, y, z, residue_id]
                ca_atoms.append(atom)
lfile.close()

if len(ca_atoms) > 0:
    q = computeQ()
    out.write(str(round(q, 3)))
    out.write(' ')
    out.write('\n')
    n_atoms = len(ca_atoms)

out.close()
