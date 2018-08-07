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

#def NewDump(dump_lines,ChainNo,ChainLength,boxInfo):
#	atoms_dump = [];

#	for l in dump_lines:
#		l=l.strip()

#	for i in range(ChainNo):
#		for j in range(ChainLength):


def calc_dihedral_angle(p1, p2, p3, p4):
    v1 = vector(p1, p2)
    v2 = vector(p2, p3)
    v3 = vector(p3, p4)
    return 180*dihedral_angle(v1, v2, v3)/3.14159265358979

if len(sys.argv)!=5 :
    print "\n.py dump_file output ChainNo ChainLength\n"
    print "Take dump_file, output psi for all residues in each snapshot.\n"
    print
    sys.exit()

dump_file = sys.argv[1]
out_file  = sys.argv[2]
ChainNo  = int (sys.argv[3])
ChainLength =int (sys.argv[4])

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


o=open(out_file,'w')
for i in range(n_snapshot):
    ## i_snapshot of the dump_file
    i_dump  = Zheng_func.get_dump_i(dump_file, i)
    ## get atom coordinates
    #ca_atoms, cb_atoms, o_atoms= Zheng_func.get_atoms_dump(i_dump)
    #print i_dump[0];
    o.write(i_dump[0]);
    o.write(i_dump[1]);
    o.write(i_dump[2]);
    o.write(i_dump[3]);
    o.write(i_dump[4]);
    o.write(i_dump[5]);
    o.write(i_dump[6]);
    o.write(i_dump[7]);
    o.write(i_dump[8]);
    ChainCenterCoord=[ 0.0 for i in range(36) ];
    for j in range(ChainNo):
        l = i_dump[ChainLength*3*j+9];l = l.split();
        #imagexref=float(l[5]);imageyref=float(l[6]);imagezref=float(l[7]);
        #ChainCenterCoord=[ 0.0 for i in range(18) ];

        for k in range(ChainLength*3):

            #ChainCenterCoord=[0.0,0.0,0.0];
            #print  ChainCenterCoord;
            l = i_dump[ChainLength*3*j+k+9];l = l.split();
            i_atom = l[0]; type =l[1];
            x = float(l[2]); y = float(l[3]); z = float(l[4])
            xnew=x; ynew=y; znew=z;labelx=0;labely=0;labelz=0;
            if k== 0 and j==0:
                xold = xnew; yold = ynew; zold = znew;
                o.write(str(i_atom)+" "+str(type)+" "+str(xold)+" "+str(yold)+" "+str(zold)+"\n");
            if k== 0 and j!=0:
                print str(xnew)+" "+str(ynew)+" "+str(znew)+"\n";
                NearChainDist=10000000.0;NearChainDistNative=100000;
                NearChainLabel=100;
                dist=1000.0;
                #print "this is j:"+ str(j);
                for label in range(j):
                #	#print "this is label:"+str(label);
                    dist=(xnew-ChainCenterCoord[3*label])*(xnew-ChainCenterCoord[3*label])+(ynew-ChainCenterCoord[3*label+1])*(ynew-ChainCenterCoord[3*label+1])+(znew-ChainCenterCoord[3*label+2])*(znew-ChainCenterCoord[3*label+2]);
                #	#print dist;
                    if dist < NearChainDistNative:
                        NearChainDistNative=dist;
                #print NearChainLabel
                dist=900;
                #xold = ChainCenterCoord[0+3*NearChainLabel];yold = ChainCenterCoord[1+3*NearChainLabel];zold = ChainCenterCoord[2+3*NearChainLabel];
                for x0 in range(-1,2):
                    for y0 in range(-1,2):
                        for z0 in range(-1,2):
                            for label in range(j):
                                xold=xnew+x0;yold=ynew+y0;zold=znew+z0;
                                dist=(xold-ChainCenterCoord[3*label])*(xold-ChainCenterCoord[3*label])+(yold-ChainCenterCoord[3*label+1])*(yold-ChainCenterCoord[3*label+1])+(zold-ChainCenterCoord[3*label+2])*(zold-ChainCenterCoord[3*label+2]);
                                if dist < NearChainDist:
                                                            NearChainDist=dist;labelx=x0;labely=y0;labelz=z0;
                                #print dist
                if abs(NearChainDist-NearChainDistNative)<0.1:
                    xold = xnew; yold = ynew; zold = znew;
                else:
                    xold=xnew+labelx;yold=ynew+labely;zold=znew+labelz;
                print ChainCenterCoord
                print str(labelx)+" "+str(labely)+" "+str(labelz);
                o.write(str(i_atom)+" "+str(type)+" "+str(xold)+" "+str(yold)+" "+str(zold)+"\n")
                #print str(xnew)+" "+str(xnew)+" "+str(xnew)+"\n"
                #ChainCenterCoord[0+3*j]=xnew/(ChainLength*3)+ChainCenterCoord[0];
                                #ChainCenterCoord[1+3*j]=ynew/(ChainLength*3)+ChainCenterCoord[1];
                                #ChainCenterCoord[2+3*j]=znew/(ChainLength*3)+ChainCenterCoord[2];
                print str(xold)+" "+str(yold)+" "+str(zold)+"\n"
            if k > 0:
                if xnew > xold and abs(xnew-xold)>0.5:
                    xnew=xnew-1.0;
                if xnew < xold and abs(xnew-xold)>0.5:
                                        xnew=xnew+1.0;
                if ynew > yold and abs(ynew-yold)>0.5:
                                        ynew=ynew-1.0;
                if ynew < yold and abs(ynew-yold)>0.5:
                                        ynew=ynew+1.0;
                if znew > zold and abs(znew-zold)>0.5:
                                        znew=znew-1.0;
                if znew < zold and abs(znew-zold)>0.5:
                                        znew=znew+1.0;
                o.write(str(i_atom)+" "+str(type)+" "+str(xnew)+" "+str(ynew)+" "+str(znew)+"\n")
                xold = xnew; yold = ynew; zold = znew;
            ChainCenterCoord[0+3*j]=xold/(ChainLength*3)+ChainCenterCoord[0+3*j];
            ChainCenterCoord[1+3*j]=yold/(ChainLength*3)+ChainCenterCoord[1+3*j];
            ChainCenterCoord[2+3*j]=zold/(ChainLength*3)+ChainCenterCoord[2+3*j];


    ## calculate N, C-prime coordinates
    #Nres = len(ca_atoms); c_atoms = []; n_atoms= []
    #n0 = [0,0,0]; n_atoms.append(n0) ## first N atom can not be calculated
    #for ires in range(Nres-1):
        #ca_xi = ca_atoms[ires][0];   ca_yi  = ca_atoms[ires][1];   ca_zi  = ca_atoms[ires][2]
        #ca_xi1= ca_atoms[ires+1][0]; ca_yi1 = ca_atoms[ires+1][1]; ca_zi1 = ca_atoms[ires+1][2]
        #o_xi = o_atoms[ires][0];   o_yi  = o_atoms[ires][1];   o_zi  = o_atoms[ires][2]

        #c_xi = ap*ca_xi + bp*ca_xi1 + cp*o_xi
        #c_yi = ap*ca_yi + bp*ca_yi1 + cp*o_yi
        #c_zi = ap*ca_zi + bp*ca_zi1 + cp*o_zi
        #atom_cp = [c_xi, c_yi, c_zi]
        #c_atoms.append(atom_cp)

        #n_xi1 =an*ca_xi + bn*ca_xi1 + cn*o_xi
        #n_yi1 =an*ca_yi + bn*ca_yi1 + cn*o_yi
        #n_zi1 =an*ca_zi + bn*ca_zi1 + cn*o_zi
        #atom_n1 =[n_xi1, n_yi1, n_zi1]
        #n_atoms.append(atom_n1)
    #c_Nres = [0,0,0]; c_atoms.append(c_Nres) ## last Cp atom can not be calculated

    ## calculate dihedral angles
    #psi=[0]*(Nres-1)
    #for ia in range(Nres-1): # psi: [0,Nres-2]; phi: [1,Nres-1]
    #	psi[ia] =calc_dihedral_angle(n_atoms[ia], ca_atoms[ia], c_atoms[ia], n_atoms[ia+1])
        #phi =calc_dihedral_angle(c_atoms[ia], n_atoms[ia+1], ca_atoms[ia+1], c_atoms[ia+1])
    ## output
    #psi_round = [ round(elem,1) for elem in psi ]
    #str_out = ' '.join(map(str,psi_round))
    #print str_out
    #print "haha"
    #out.write(str_out+'\n')

o.close()
