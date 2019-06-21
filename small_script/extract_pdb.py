from Bio.PDB import *
import os
class MyChainSelect(Select):
    def __init__(self, chain, start, end):
        super(MyChainSelect, self).__init__()
        self.chain=chain
        self.start=start
        self.end=end

    def accept_chain(self, chain):
        if chain.get_id()==self.chain:
            return True
        else:
            return False

    def accept_residue(self, residue):
        if residue.get_id()[1] >= self.start and residue.get_id()[1] <= self.end:
            return True
        else:
            return False

def extract_pdb(pwd, protein, chain, residue_start, residue_end):
    parser = PDBParser(PERMISSIVE=1,QUIET=True)
    structure_id = protein
    fileName = pwd + 'original_pdbs/' + structure_id + ".pdb"
    outFileName = pwd + "extracted/" + structure_id + ".pdb"
    structure = parser.get_structure(structure_id, fileName)
    io = PDBIO()
    io.set_structure(structure)
    io.save(pwd + "tmp.pdb", MyChainSelect(chain, residue_start, residue_end))
    os.system("~/opt/small_script/pdb_reres.py -1 " + pwd + "tmp.pdb > " + outFileName)
