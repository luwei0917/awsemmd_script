#!/usr/bin/python

# ----------------------------------------------------------------------
# Copyright (2010) Aram Davtyan and Garegin Papoian

# Papoian's Group, University of Maryland at Collage Park
# http://papoian.chem.umd.edu/

# Last Update: 03/04/2011
# ----------------------------------------------------------------------

import sys, os
import numpy
from Bio.SVDSuperimposer import SVDSuperimposer
from VectorAlgebra import *

#from Bio.PDB.PDBParser import PDBParser

atom_type = {'1' : 'C', '2' : 'N', '3' : 'O', '4' : 'C', '5' : 'H', '6' : 'C'}
atom_desc = {'1' : 'C-Alpha', '2' : 'N', '3' : 'O', '4' : 'C-Beta', '5' : 'H-Beta', '6' : 'C-Prime'}
PDB_type = {'1' : 'CA', '2' : 'N', '3' : 'O', '4' : 'CB', '5' : 'HB', '6' : 'C' }

class PDB_Atom:
	no = 0
	ty = ''
	res = 'UNK'
	res_no = 0
	x = 0.0
	y = 0.0
	z = 0.0
	atm = 'C'

	def __init__(self, no, ty, res, res_no, x, y, z, atm):
		self.no = no
		self.ty = ty
		self.res = res
		self.res_no = res_no
		self.x = x
		self.y = y
		self.z = z
		self.atm = atm

	def write_(self, f):
		f.write('ATOM')
		f.write(('       '+str(self.no))[-7:])
		f.write('  ')
		f.write((self.ty+'    ')[:4])
		f.write(self.res)
		f.write(' ')
		f.write('T')
		f.write(('    '+str(self.res_no))[-4:])
		f.write(('            '+str(round(self.x,3)))[-12:])
		f.write(('        '+str(round(self.y,3)))[-8:])
		f.write(('        '+str(round(self.z,3)))[-8:])
		f.write('  1.00')
		f.write('  0.00')
		f.write(('            '+self.atm)[-12:]+'  ')
		f.write('\n')

class Atom:
	No = 0
	ty = ''
	x = 0.0
	y = 0.0
	z = 0.0
	desc = ''

	def __init__(self, No, ty, No_m, x, y, z, desc=''):
		self.No = No
		self.ty = ty
		self.No_m = No_m
		self.x = x
		self.y = y
		self.z = z
		self.desc = desc

	def write_(self, f):
		f.write(str(self.No))
		f.write(' ')
		f.write(PDB_type[self.No_m])
		f.write(' ')
		f.write(str(round(self.x,8)))
		f.write(' ')
		f.write(str(round(self.y,8)))
		f.write(' ')
		f.write(str(round(self.z,8)))
		f.write(' ')
		f.write(self.desc)
		f.write('\n')

if len(sys.argv)!=7 and len(sys.argv)!=8 :
	print "# of arguments", len(sys.argv)
	print "\nCalcQValue.py PDB_Id Frag_PDB_path FragPDB+chain i_start j_start length (1 for rmsd calculation) \n"
	exit()

struct_id = sys.argv[1]
if struct_id[-4:].lower()==".pdb":
	pdb_file = struct_id
else:
	pdb_file = struct_id + ".pdb"

#myhome = os.environ.get("HOME")
#pdbDir = myhome + "/opt/script/PDBs/"
Frag_struct_id = sys.argv[3]

if len(Frag_struct_id) < 5:
  print "\nCalcQValue.py PDB_Id Frag_PDB_path FragPDB+chain i_start j_start length\n"
  print "format of FragPDB+chain, like 1r69a "
  exit()

pdbDir = sys.argv[2]
Frag_pdbID = Frag_struct_id[0:4].upper()
Frag_chain_name = Frag_struct_id[4].upper()
Frag_pdb_file = pdbDir + Frag_pdbID + ".pdb"

#output_file = ""
#if len(sys.argv)>3: output_file = sys.argv[3]

sigma_exp = 0.15
res_flag = 0

res_flag  = 1
res_Start = int(sys.argv[4])
res_len   = int(sys.argv[6])
rmsd_flag = 0
if len(sys.argv) == 8 :
  rmsd_flag = int(sys.argv[7])
  if rmsd_flag != 1 :
  	print "rmsd_flag has to be 1 !!"
	exit()

#res_End = res_Start + res_len - 1
Frag_Start = int(sys.argv[5])
#Frag_End   = res_len + Frag_Start - 1
#print "res_Start:", res_Start, "Frag_Start:", Frag_Start

n_atoms = 0
i_atom = 0
item = ''
step = 0
ca_atoms_pdb = []
ca_atoms = []
box = []
A = []
sigma = []
sigma_sq = []

#out = open(output_file, 'w')

from Bio.PDB.PDBParser import PDBParser

p = PDBParser(PERMISSIVE=1)

def compute_frag_RMSD(res_len):
        if len(ca_atoms)!=len(ca_atoms_pdb):
		print "Error. Length mismatch! target:frag", len(ca_atoms_pdb), len(ca_atoms)
		return 0
        l = len(ca_atoms)
	N = res_len
	if l != N :
		print "atom list length mismatches the fragment length!", str(l), str(N)
		return 0

        fixed_coord  = numpy.zeros((l, 3))
        moving_coord = numpy.zeros((l, 3))

        for i in range(0, l):
                fixed_coord[i]  = numpy.array ([ca_atoms_pdb[i][0], ca_atoms_pdb[i][1], ca_atoms_pdb[i][2]])
                moving_coord[i] = numpy.array ([ca_atoms[i][0], ca_atoms[i][1], ca_atoms[i][2]])
        sup = SVDSuperimposer()
        sup.set(fixed_coord, moving_coord)
        sup.run()
        rms = sup.get_rms()
        return rms

def computeQ1(res_len):
	#print "res_len:", res_len
	#print "target:frag", len(ca_atoms_pdb), len(ca_atoms)
	if len(ca_atoms)!=len(ca_atoms_pdb):
		print "Error. Length mismatch! target:frag", len(ca_atoms_pdb), len(ca_atoms)
		return 0
	Q = 0
	N = res_len

	for ia in range(0, res_len-2):
	  for ja in range(ia+2, res_len):
		  #print ia, ja
		  r = vabs(vector(ca_atoms[ia], ca_atoms[ja]))
		  rn = vabs(vector(ca_atoms_pdb[ia], ca_atoms_pdb[ja]))
		  dr = r - rn
		  Q = Q + exp(-dr*dr/(2*sigma_sq[ja-ia]));
	Q = 2*Q/((N-1)*(N-2))
	return Q

s = p.get_structure(struct_id, pdb_file)
chains = s[0].get_list()
chain = chains[0]
count=0
shift = 0
first_flag = 0
for res in chain:
	res_index = res.get_id()[1]
	#print "res_index: ", res_index, "first_flag: ", first_flag

	#pdb.fasta file is modified, shift is needed to get back pbd index
	if first_flag == 0 and res_index != 1:
	    shift = res_index - 1
	#   print "shift:",shift

	first_flag = 1
	res_index = res_index - shift
	if (res_index < res_Start):
		continue
	if (count >= res_len ):
		break

	is_regular_res = res.has_id('CA') and res.has_id('O')
	res_id = res.get_id()[0]
        if (res_id ==' ' or res_id =='H_MSE' or res_id == 'H_M3L' or res_id == 'H_CAS') and is_regular_res:
		#print "Add res_index:", res_index + shift, "resname:", res.get_resname()
		ca_atoms_pdb.append(res['CA'].get_coord())
		#print 'res: ', res.get_resname(), res.get_id()[1]
		count += 1
	else:
		print "res_id:", res_id, "res_index: ", res_index, "is_regular_res: ", is_regular_res
		#exit()

s = p.get_structure(Frag_pdbID, Frag_pdb_file)
chains = s[0].get_list()
for chain in chains:
  if chain.get_id() == Frag_chain_name:
    count_frag=0
    for res in chain:
	res_index = res.get_id()[1]
	if (res_index < Frag_Start):
		continue
	if (count_frag >= res_len ):
		break

	count_frag += 1
	is_regular_res = res.has_id('CA')
	res_id = res.get_id()[0]
        if (res_id ==' ' or res_id =='H_MSE' or res_id == 'H_M3L' or res_id == 'H_CAS') and is_regular_res:
		ca_atoms.append(res['CA'].get_coord())
	else:
		print "res_id:", res_id, "res_index: ", res_index, "is_regular_res: ", is_regular_res

for i in range(0, len(ca_atoms_pdb)+1):
	#sigma.append( (1+i)**sigma_exp )
	#sigma_sq.append(sigma[-1]*sigma[-1])
	sigma_tmp = (1+i)**sigma_exp
	sigma_sq.append(sigma_tmp*sigma_tmp)

if rmsd_flag == 1:
	rmsd = compute_frag_RMSD(res_len)
	print round(rmsd,3)
else :
	Q = computeQ1(res_len)
	strs=str(round(Q,3))
	sys.stdout.write(strs+"\n")
	#print round(Q,3)
