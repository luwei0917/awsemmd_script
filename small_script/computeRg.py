from Bio.PDB.PDBParser import PDBParser
import argparse
parser = argparse.ArgumentParser(description="Compute Rg of pdb")
parser.add_argument("pdb", help="pdb file")
args = parser.parse_args()

def computeRg(pdb_file, chain="A"):
    # compute Radius of gyration
    # pdb_file = f"/Users/weilu/Research/server/feb_2019/iterative_optimization_new_temp_range/all_simulations/{p}/{p}/crystal_structure.pdb"
    chain_name = chain
    parser = PDBParser()
    structure = parser.get_structure('X', pdb_file)
    chain = list(structure[0][chain_name])
    n = len(chain)
    regular_res_list = [res for res in chain  if res.get_id()[0] == ' ']
    cutoff = 15
    for residue in regular_res_list:
        if residue.get_id()[0] == ' ' and abs(residue["CA"].get_vector()[-1]) < cutoff:
            print(residue.get_id()[1])
    rg = 0.0
    for i, residue_i in enumerate(regular_res_list):
        for j, residue_j in enumerate(regular_res_list[i+1:]):
            try:
                r = residue_i["CA"] - residue_j["CA"]
            except:
                print(residue_i, residue_j)
            rg += r**2
    return (rg/(n**2))**0.5

rg = computeRg(args.pdb)
print(rg)


def cylindrical_rg_bias_term(oa, k_rg=4.184, rg0=0, atomGroup=-1, forceGroup=27):
    nres, ca = oa.nres, oa.ca
    if atomGroup == -1:
        group = list(range(nres))
    else:
        group = atomGroup     # atomGroup = [0, 1, 10, 12]  means include residue 1, 2, 11, 13.
    n = len(group)
    rg_square = CustomBondForce("1/normalization*(x^2+y^2)")
    # rg = CustomBondForce("1")
    rg_square.addGlobalParameter("normalization", n*n)
    for i in group:
        for j in group:
            if j <= i:
                continue
            rg_square.addBond(ca[i], ca[j], [])
    rg = CustomCVForce(f"{k_rg}*(rg_square^0.5-{rg0})^2")
    rg.addCollectiveVariable("rg_square", rg_square)
    rg.setForceGroup(forceGroup)
    return rg
