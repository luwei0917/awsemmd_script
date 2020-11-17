from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBIO import Select
import argparse
parser = argparse.ArgumentParser(description="python extract_membrane_part.py 1su4.pdb 1su4_membrane.pdb ")
parser.add_argument("source_pdb", help="pdb file")
parser.add_argument("to_pdb", help="pdb file")
parser.add_argument("-c", "--cutoff", default=15, help="default is 15A")
args = parser.parse_args()


def extractTransmembrane(toLocation, location, cutoff=15):
    x = PDBParser(QUIET=True).get_structure("x", location)

    class Transmembrane(Select):
        def accept_residue(self, residue):
            try:
                if residue.get_id()[0] == ' ' and abs(float(residue["CA"].get_vector()[-1])) < cutoff:
                    return 1
                else:
                    return 0
            except Exception as e:
                print(e)
                print(residue["CA"].get_vector())
                print(residue)
                return 0
    io = PDBIO()
    io.set_structure(x)
    io.save(toLocation, Transmembrane())

# pdb_file = "1su4.pdb"
# pdb_file_2 = "1su4_membrane.pdb"
# extractTransmembrane(pdb_file_2, pdb_file)
extractTransmembrane(args.to_pdb, args.source_pdb, cutoff=float(args.cutoff))