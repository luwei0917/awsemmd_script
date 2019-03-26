from pyCodeLib import *
import warnings
import glob
import re
import argparse
warnings.filterwarnings('ignore')

parser = argparse.ArgumentParser(
    description="compute phi for protein list")

parser.add_argument("proteins", help="The name of the protein list")
parser.add_argument("-m", "--mode", type=int, default=0)
args = parser.parse_args()


# transfer database
def transferPDB(file=None, source="database/dompdb_cleaned/"):
    if file is None:
        p_list = glob.glob(source + "*")
    else:
        with open(file) as f:
            content = f.readlines()
        # you may also want to remove whitespace characters like `\n` at the end of each line
        p_list = [source + x.strip() for x in content]
    for p in p_list:
        new_p = p.replace(source, "database/dompdb/")
        os.system(f"cp {p} {new_p}.pdb")
# add location
def addLocation(source, target):
    with open(target, "w") as out:
        with open(source, "r") as f:
            for l in f:
                nl = "database/dompdb/" + l
                out.write(nl)
def convertTo4cName():
    with open("proteins_4c_name_list.txt", "w") as out:
        with open("proteins_name_list.txt", "r") as f:
            for l in f:
                nl = l[:4]
                out.write(nl+"\n")
# extract seq info from fa file
def extractSeqFromFa(size=7, source="database/cath-dataset-nonredundant-S20Clean.atom.fa"):
    os.system("mkdir -p database/S20_seq")
    with open(source, "r") as f:
        count = 0
        for l in f:
            if count % 2 == 0:
                #extract protein id
                assert(l[0] == ">")
    #             print(l)
                if size == 7:
                    name = re.search('>cath\|(.*)\|(\w{7})\/(.*)', l).group(2)
                if size == 4:
                    name = re.search('>cath\|(.*)\|(\w{4})(.*)\/(.*)', l).group(2)
    #             name = "test"
    #             print(name)
            else:
                assert(l[0] != ">")
    #             print(l)
                with open(f"database/S20_seq/{name}.seq", "w") as out:
                    out.write(l)
            count += 1


# convertTo4cName()
# extractSeqFromFa()

def computePhis(proteins, multiSeq=False, sampleK=1000):
    # if addGylcines:
    proteins_location = "".join(proteins.split("/")[:-1]) + "/location_" + proteins.split("/")[-1]
    # transferPDB(proteins)
    addLocation(proteins, proteins_location)
    # os.chdir('database')
    add_virtual_glycines_list(proteins_location)
    # generate_decoy_sequences(proteins, methods=['shuffle'], databaseLocation="../../")
    # evaluate_phis_over_training_set_for_native_structures_Wei(proteins, "phi_list.txt", decoy_method='shuffle', max_decoys=1e+10, tm_only=False, num_processors=1)
    if multiSeq:
        evaluate_phis_over_training_set_for_native_structures_Wei(proteins, "phi_list.txt", decoy_method='multiShuffle', max_decoys=1e+10, tm_only=False, num_processors=1, multi_seq=True, sampleK=sampleK)
# extractSeqFromFa()

def computePhisForDecoys(proteins, **kwargs):
    proteins_location = "".join(proteins.split("/")[:-1]) + "/location_" + proteins.split("/")[-1]
    # transferPDB(proteins)
    # addLocation(proteins, proteins_location)
    # # os.chdir('database')
    # add_virtual_glycines_list(proteins_location)
    # generate_decoy_structures(proteins, methods=['lammps'], num_decoys=[10], databaseLocation="../../")
    evaluate_phis_over_training_set_for_decoy_structures_Wei(proteins, "phi_list.txt", decoy_method='lammps', max_decoys=1e+5, tm_only=False, num_processors=1, **kwargs)

if args.mode == 0:
    computePhis(args.proteins)
elif args.mode == 1:
    computePhisForDecoys(args.proteins)
elif args.mode == 2:
    computePhisForDecoys(args.proteins, withBiased=True)
elif args.mode == 3:
    computePhisForDecoys(args.proteins, withBiased=True, mode=1)
if args.mode == 4:
    computePhis(args.proteins, multiSeq=True, sampleK=1000)
