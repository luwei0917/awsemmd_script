from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBIO import Select
import argparse
parser = argparse.ArgumentParser(description="Compute Rg of pdb, python extract_membrane_part.py 1su4.pdb 1su4_membrane.pdb ")
parser.add_argument("source_pdb", help="pdb file")
parser.add_argument("to_pdb", help="pdb file")
args = parser.parse_args()


def extractTransmembrane(toLocation, location, cutoff=15):
    x = PDBParser().get_structure("x", location)

    class Transmembrane(Select):
        def accept_residue(self, residue):
            if residue.get_id()[0] == ' ' and abs(residue["CA"].get_vector()[-1]) < cutoff:
                return 1
            else:
                return 0

    io = PDBIO()
    io.set_structure(x)
    io.save(toLocation, Transmembrane())

# pdb_file = "1su4.pdb"
# pdb_file_2 = "1su4_membrane.pdb"
# extractTransmembrane(pdb_file_2, pdb_file)
extractTransmembrane(args.to_pdb, args.source_pdb)