#!/usr/bin/python
# ----------------------------------------------------------------------
# Copyright (2010) Aram Davtyan and Garegin Papoian
# Papoian's Group, University of Maryland at Collage Park
# http://papoian.chem.umd.edu/
# Modified by Weihua Zheng 03/2014 
# Last Update: 03/04/2011
# ----------------------------------------------------------------------

import sys
import Zheng_func
from VectorAlgebra import *

def calc_dihedral_angle(p1, p2, p3, p4):
    v1 = vector(p1, p2)
    v2 = vector(p2, p3)
    v3 = vector(p3, p4)
    return 180*dihedral_angle(v1, v2, v3)/3.14159265358979

if len(sys.argv)!=3 :
	print "\n.py dump_file output\n"
	print "Take dump_file, output psi for all residues in each snapshot.\n"
	print
	sys.exit()

dump_file = sys.argv[1]
out_file  = sys.argv[2]

an = 0.4831806                                                                                                        
bn = 0.7032820                                                                                                        
cn = -0.1864262                                                                                                       
ap = 0.4436538                                                                                                        
bp = 0.2352006                                                                                                        
cp = 0.3211455

#Get n_snapshot
file_len       = Zheng_func.get_file_len(dump_file)
nline_snapshot = Zheng_func.get_nline_snapshot(dump_file)
n_snapshot     = file_len / nline_snapshot
print "n_snapshot of the dump file is = ", n_snapshot

out=open(out_file,'w')
for i in range(n_snapshot):
	## i_snapshot of the dump_file
	i_dump  = Zheng_func.get_dump_i(dump_file, i) 
	## get atom coordinates
	ca_atoms, cb_atoms, o_atoms= Zheng_func.get_atoms_dump(i_dump) 

	## calculate N, C-prime coordinates
	Nres = len(ca_atoms); c_atoms = []; n_atoms= []
	n0 = [0,0,0]; n_atoms.append(n0) ## first N atom can not be calculated
	for ires in range(Nres-1):
		ca_xi = ca_atoms[ires][0];   ca_yi  = ca_atoms[ires][1];   ca_zi  = ca_atoms[ires][2]
		ca_xi1= ca_atoms[ires+1][0]; ca_yi1 = ca_atoms[ires+1][1]; ca_zi1 = ca_atoms[ires+1][2]
		o_xi = o_atoms[ires][0];   o_yi  = o_atoms[ires][1];   o_zi  = o_atoms[ires][2]

		c_xi = ap*ca_xi + bp*ca_xi1 + cp*o_xi
		c_yi = ap*ca_yi + bp*ca_yi1 + cp*o_yi
		c_zi = ap*ca_zi + bp*ca_zi1 + cp*o_zi
		atom_cp = [c_xi, c_yi, c_zi]
		c_atoms.append(atom_cp)

		n_xi1 =an*ca_xi + bn*ca_xi1 + cn*o_xi 
		n_yi1 =an*ca_yi + bn*ca_yi1 + cn*o_yi 
		n_zi1 =an*ca_zi + bn*ca_zi1 + cn*o_zi 
		atom_n1 =[n_xi1, n_yi1, n_zi1]
		n_atoms.append(atom_n1)
	c_Nres = [0,0,0]; c_atoms.append(c_Nres) ## last Cp atom can not be calculated
		
	## calculate dihedral angles
	psi=[0]*(Nres-1)
	for ia in range(Nres-1): # psi: [0,Nres-2]; phi: [1,Nres-1]
		psi[ia] =calc_dihedral_angle(n_atoms[ia], ca_atoms[ia], c_atoms[ia], n_atoms[ia+1])
		#phi =calc_dihedral_angle(c_atoms[ia], n_atoms[ia+1], ca_atoms[ia+1], c_atoms[ia+1])
	## output
	psi_round = [ round(elem,1) for elem in psi ]
	str_out = ' '.join(map(str,psi_round))
	#print str_out
	#print "haha"
	out.write(str_out+'\n')

out.close()
