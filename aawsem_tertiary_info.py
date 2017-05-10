#!/usr/bin/env python3
import os
import argparse
import sys
from time import sleep
import subprocess
from Bio.PDB import *
import numpy as np

parser = argparse.ArgumentParser(description="This is my aawsem project headquarter")
parser.add_argument("-t", "--test", help="Test run", action="store_true", default=False)
parser.add_argument("--fix", action="store_true", default=False)
parser.add_argument("--move", action="store_true", default=False)
parser.add_argument("-d", "--debug", action="store_true", default=False)
parser.add_argument("-e", "--extract", action="store_true", default=False)
# parser.add_argument("protein", help="the name of protein file")
args = parser.parse_args()

# try:
#     protein_name,_ = args.protein.split('.')
# except:
#     try:
#         protein_name,_ = args.protein.split('/')
#     except:
#         protein_name = args.protein

if(args.debug):
    do = print
    cd = print
else:
    do = os.system
    cd = os.chdir

def download(pdbfull):
    # pdbfull = "4OUQA"
    pdbID = pdbfull[0:4].lower()
    pdbIDsecond = pdbfull[1:2].lower()
    pdbIDthird = pdbfull[2:3].lower()
    chainID = pdbfull[4:5].lower()

    exeline = "wget ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/"
    exeline += pdbIDsecond + pdbIDthird + "/pdb" + pdbID + ".ent.gz"
    os.system(exeline)
    os.system("nice gunzip pdb" + pdbID + ".ent.gz; mv pdb" +
              pdbID + ".ent " + pdbID.upper() + ".pdb")
    print("nice gunzip pdb" + pdbID + ".ent.gz; mv pdb" +
              pdbID + ".ent " + pdbID.upper() + ".pdb")
    if not os.path.isfile(pdbID.upper() + ".pdb"):
        print(":::Cannot build PDB for PDB ID, failed to download:" + pdbID.upper())


class Lammps(Select):
    # def accept_atom(self, atom):
    #     if atom.get_name() == 'CA':
    #         return 1
    #     else:
    #         return 0

    def accept_residue(self, residue):
        # print(residue.get_id())[1]
        # return 0
        if residue.get_id()[1] == -1:
            return 0
        else:
            return 1

database = "~/opt/database/cullpdb_pc95_res3.0_R0.3_d170428_chains35871"
# database = "cullpdb_pc80_res3.0_R1.0_d160504_chains29712"
# # database = "cullpdb_pc90_res3.0_R1.0_d170427_chains34250"
#
#
exeline = "psiblast -num_iterations 5 -comp_based_stats 0 -word_size 2 -evalue 10000"
exeline += " -outfmt '6 sseqid qstart qend sstart send qseq sseq length gaps bitscore evalue' -matrix BLOSUM62 -threshold 9 -window_size 0"
exeline += " -db " + database + " -query T0766.fasta"



# exeline = "psiblast -num_iterations 1 -word_size 3 -evalue 0.005"
# exeline = "psiblast -num_iterations 1 -word_size 3 -evalue 0.5"
# # exeline += " -matrix BLOSUM62 -db " + \
# #     database + " -query fragment.fasta"
# # exeline += " -outfmt '6 sseqid slen bitscore score evalue pident' -matrix BLOSUM62 -db " + \
# #     database + " -query fragment.fasta"
#
# exeline += " -outfmt '6 sseqid qstart qend sstart send bitscore score evalue length pident qseq sseq' -matrix BLOSUM62 -db " + \
#     database + " -query T0803.fasta"
# do(exeline)
class sequence:
    def __init__(self, start, end, query_start):
        self.start = start
        self.end = end
        self.query_start = query_start

    def __str__(self):
        return "{}, {}, {}".format(self.start, self.end, self.query_start)

if(args.extract):
    # download PDBs if not exist    ##from script 'pdbget' (original author
    # unknown)
    subject_seq_list = []
    homolog = "4OUQA".lower()
    with open("log.mem", "r") as f:
        for line in f:
            tmp = line.split()
            pdbfull = tmp[4][:5]
            query_start = int(tmp[5])
            subject_start = int(tmp[6])
            segment_length = int(tmp[7])
            rank = tmp[1]
            if rank != "1":
                continue
            if pdbfull != homolog:
                continue
            # print(pdbfull)
            added = False
            for seq in subject_seq_list:
                if subject_start >= seq.start and subject_start <= seq.end:
                    seq.end = max(subject_start+segment_length-1, seq.end)
                    added = True
                    break
            if not added:
                subject_seq_list.append(sequence(subject_start, subject_start+segment_length-1, query_start))
            # print(tmp)

        # print(subject_seq_list)
    for i in subject_seq_list:
        print(i)

    # download("1x3ka")
    # pdbDir = os.path.expanduser("~") + "/opt/script/PDBs/"
    # pdbfull = "4OUQA"
    pdbID = homolog[:4]
    if not os.path.isfile(pdbID.upper() + ".pdb"):
        download(homolog)
    parser = PDBParser()
    structure = parser.get_structure('test', homolog[:4].upper() + ".pdb")
    # # print(structure)
    for model in structure:
        for chain in model:
            for residue in chain:
                is_regular_res = residue.has_id('CA')
                id = residue.id
                # print(residue, is_regular_res)
                if id[0] != ' ' or not is_regular_res:
                    residue.id = (' ', -1, ' ')
                    # chain.detach_child(id)
                # chain.detach_child(id)
                else:
                    subject_start = id[1]
                    # print(subject_start)
                    inside = False
                    for seq in subject_seq_list:
                        if subject_start >= seq.start and subject_start <= seq.end:
                            residue.id = (' ', subject_start - seq.start + seq.query_start, ' ')
                            inside = True
                            # print(residue.id)
                            break
                    if not inside:
                        residue.id = (' ', -1, ' ')
                        # chain.detach_child(id)
            # if len(chain) == 0:
            #     model.detach_child(chain.id)
    #
    io = PDBIO()
    io.set_structure(structure)
    io.save('ca_only.pdb', Lammps())
    # io.save('ca_only.pdb')

if(args.test):
    # folder_list = ["T0766", "T0778", "T0782", "T0792", "T0815", "T0833", "T0844", "T0842", "T0846", "T0803"]
    # folder_list = ["T0792", "T0778", "T0782", "T0833", "T0844"]
    # folder_list = ["T0792"]
    # folder_list = ["T0766", "T0778", "T0782", "T0792"]
    folder_list = ["T0766"]
    file_name = "/Users/weilu/opt/AAWSEM/HO_Mem/T0766/homologues.txt"
    with open(file_name, "r") as f:
        for line in f:
            name, *middle, pident = line.split()
            print(name, pident)
            # For now, we choose the one with highest Percentage of identical matches but not higher than 95%.

    # for protein_name in folder_list:
    #     os.chdir(protein_name)
    #     n = 20
    #     os.chdir("simulation")
    #     for i in range(n):
    #         os.chdir(str(i))
    #         os.system("cp ~/opt/AAWSEM/aawsem.slurm .")
    #         os.system(
    #             "sed -i.bak 's/PROTEIN/'" +
    #             protein_name +
    #             "'/g' aawsem.slurm")
    #         os.system("sbatch aawsem.slurm")
    #         os.chdir("..")
    #     # os.system("run.py -n 20 -o "+protein_name+"\/")
    #     os.chdir("../..")
