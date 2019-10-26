#!/usr/bin/env python3
import pandas as pd
import numpy as np
import os
import argparse
import sys

parser = argparse.ArgumentParser(description="convert to simulation format")

parser.add_argument("source", help="The file you want to convert")
parser.add_argument("folder", help="The name of folder where you want to save the formated gamma.")
# parser.add_argument("-l", "--label", type=str, default="label")
# parser.add_argument("-t", "--test", action="store_true", default=False)
args = parser.parse_args()

# if args.test:
#     do = print
# else:
#     do = os.system
with open('cmd_gg_server.txt', 'a') as f:
    f.write(' '.join(sys.argv))
    f.write('\n')


do = os.system
def gamma_format_convertion_iteration_to_simulation(iteration_gamma, gamma_for_simulation, burial_gamma_for_simulation=None, debye=None):
    from Bio.PDB.Polypeptide import one_to_three, three_to_one
    res_type_map = {
        'A': 0,
        'C': 4,
        'D': 3,
        'E': 6,
        'F': 13,
        'G': 7,
        'H': 8,
        'I': 9,
        'K': 11,
        'L': 10,
        'M': 12,
        'N': 2,
        'P': 14,
        'Q': 5,
        'R': 1,
        'S': 15,
        'T': 16,
        'V': 19,
        'W': 17,
        'Y': 18
    }
    # res_type_map = gamma_se_map_1_letter = {   'A': 0,  'R': 1,  'N': 2,  'D': 3,  'C': 4,
    #                             'Q': 5,  'E': 6,  'G': 7,  'H': 8,  'I': 9,
    #                             'L': 10, 'K': 11, 'M': 12, 'F': 13, 'P': 14,
    #                             'S': 15, 'T': 16, 'W': 17, 'Y': 18, 'V': 19}
    res_type_map_letters = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G',
                            'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

    inverse_res_type_map = dict(list(zip(list(range(20)), res_type_map_letters)))

    # gamma_location = "/Users/weilu/Research/server_backup/jan_2019/optimization/gammas_dec30/cath-dataset-nonredundant-S20Clean_phi_pairwise_contact_well4.5_6.5_5.0_10phi_density_mediated_contact_well6.5_9.5_5.0_10_2.6_7.0_gamma"
    # gamma_for_simulation = "/Users/weilu/Research/server_backup/jan_2019/optimization/iteration_gamma.dat"
    gamma = iteration_gamma
    gamma = -gamma  # caused by tradition.
    # convert gamma to gamma used by simulation
    with open(gamma_for_simulation, "w") as out:
        c = 0
        for i in range(20):
            for j in range(i, 20):
                out.write(f"{gamma[c]:<.5f} {gamma[c]:10.5f}\n")
                c += 1
        out.write("\n")
        for i in range(20):
            for j in range(i, 20):
                # protein, water
                out.write(f"{gamma[c]:<.5f} {gamma[c+210]:10.5f}\n")
                c += 1
    if burial_gamma_for_simulation:
        rhoGamma = pd.DataFrame(gamma[630:690].reshape(3,20).T, columns=["rho1", "rho2", "rho3"]).reset_index()
        rhoGamma["oneLetter"] = rhoGamma["index"].apply(lambda x: inverse_res_type_map[x])
        rhoGamma["Residue"] = rhoGamma["index"].apply(lambda x: one_to_three(inverse_res_type_map[x]))
        rhoGamma = rhoGamma[["Residue", "rho1", "rho2", "rho3", "index", "oneLetter"]]
        g = rhoGamma[["rho1", "rho2", "rho3"]].values
        np.savetxt(burial_gamma_for_simulation, g, fmt='%7.4f')
    if debye:
        # np.savetxt(debye, [iteration_gamma[690]], fmt='%7.4f')
        try:
            np.savetxt(debye, [iteration_gamma[690]], fmt='%7.4f')
        except:
            # print("no debye")
            pass


# iteration = "trial_2"
pre = args.folder + "/"
do(f"mkdir -p {pre}")
# iter_gamma = np.loadtxt("trial_2_mixed_original_and_cutoff400_impose_Aprime_constraint_90")
iter_gamma = np.loadtxt(args.source)
gamma_for_simulation = pre + f"gamma.dat"
burial_gamma_for_simulation = pre + f"burial_gamma.dat"
debye = pre + f"k_debye.dat"
# gamma_format_convertion_iteration_to_simulation(iter_gamma, gamma_for_simulation, burial_gamma_for_simulation=burial_gamma_for_simulation)
gamma_format_convertion_iteration_to_simulation(iter_gamma, gamma_for_simulation, burial_gamma_for_simulation=burial_gamma_for_simulation, debye=debye)
# do("mv original_*.dat for_simulation/")
