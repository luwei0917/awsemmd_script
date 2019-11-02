#!/usr/bin/env python
from pyCodeLib import *
import warnings
import glob
import re
import numpy as np
import pandas as pd
from Bio.PDB.Polypeptide import one_to_three
import argparse
import os

def plot_contact_well(gammas, ax, fig, invert_sign=True, fix_colorbar=True, inferBound=False,
                        vmin=-0.3, vmax=0.3, fix_confidence_colorbar=True, confidence_vmin=0,
                        confidence_vmax=1.0, plot_confidence=False, confidence_lower=None, confidence_upper=None):
    size = 20
    interaction_matrix = np.zeros((size, size))
    i_content = 0
    for i in range(size):
        for j in range(i, size):
            index1 = hydrophobicity_map[inverse_res_type_map[i]]
            index2 = hydrophobicity_map[inverse_res_type_map[j]]
            interaction_matrix[index1][index2] = gammas[i_content]
            interaction_matrix[index2][index1] = gammas[i_content]
            i_content += 1


    # The minus sign is here to be consistent with the way AWSEM thinks about gammas
    if invert_sign:
        interaction_matrix *= -1

    if inferBound:
        vmin = np.min(interaction_matrix)
        vmax = np.max(interaction_matrix)

    if fix_colorbar:
        cax = ax.pcolor(interaction_matrix, vmin=vmin,
                        vmax=vmax, cmap="bwr")
        # cax = ax.pcolor(interaction_matrix, vmin=vmin,
        #                 vmax=vmax, cmap="coolwarm")
        # cax = ax.pcolor(interaction_matrix, vmin=vmin,
        #                 vmax=vmax, cmap="jet")
    else:
        cax = ax.pcolor(interaction_matrix, cmap="RdBu_r")
    fig.colorbar(cax)

    # put the major ticks at the middle of each cell
    ax.set_yticks(np.arange(interaction_matrix.shape[0]) + 0.5, minor=False)
    ax.set_xticks(np.arange(interaction_matrix.shape[1]) + 0.5, minor=False)

    ax.set_xticklabels(hydrophobicity_letters)
    ax.set_yticklabels(hydrophobicity_letters)

    if plot_confidence:
        confidence_interval_size = confidence_upper - confidence_lower
        confidence_matrix = np.zeros((size, size))
        i_content = 0
        for i in range(size):
            for j in range(i, size):
                index1 = hydrophobicity_map[inverse_res_type_map[i]]
                index2 = hydrophobicity_map[inverse_res_type_map[j]]
                confidence_matrix[index1][index2] = confidence_interval_size[i_content]
                confidence_matrix[index2][index1] = confidence_interval_size[i_content]
                i_content += 1

        fig = plt.figure()
        ax = fig.add_subplot(111)
        if fix_confidence_colorbar:
            cax = ax.pcolor(confidence_matrix, vmin=confidence_vmin,
                            vmax=confidence_vmax, cmap="RdBu_r")
        else:
            cax = ax.pcolor(confidence_matrix, cmap="RdBu_r")
        fig.colorbar(cax)

        # put the major ticks at the middle of each cell
        ax.set_yticks(np.arange(confidence_matrix.shape[0]) + 0.5, minor=False)
        ax.set_xticks(np.arange(confidence_matrix.shape[1]) + 0.5, minor=False)

        ax.set_xticklabels(hydrophobicity_letters)
        ax.set_yticklabels(hydrophobicity_letters)


def plot_contact_well_all(gamma_complete, title=None, **kwargs):

    fig = plt.figure()
    if title is None:
        title = os.path.dirname(args.location) + "\n" + os.path.basename(args.location)
    st = fig.suptitle(title, fontsize="x-large")
    # st = fig.suptitle(args.location, fontsize="x-large")
    ax = fig.add_subplot(131)
    gammas = gamma_complete[:210]
    _ = ax.set_title("Direct")
    plot_contact_well(gammas, ax, fig, **kwargs)

    ax = fig.add_subplot(132)
    gammas = gamma_complete[210:420]
    _ = ax.set_title("Protein")
    plot_contact_well(gammas, ax, fig, **kwargs)

    ax = fig.add_subplot(133)
    gammas = gamma_complete[420:630]
    _ = ax.set_title("Water")
    plot_contact_well(gammas, ax, fig, **kwargs)

    fig.tight_layout()
    # shift subplots down:
    st.set_y(0.95)
    fig.subplots_adjust(top=0.8)
#     plt.savefig('direct_contact.pdf')
    # plt.show()


parser = argparse.ArgumentParser(description="This is my playground for current project")
# parser.add_argument("protein", help="the name of protein")
# parser.add_argument("template", help="the name of template file")
parser.add_argument("location", help="gamma source ")
parser.add_argument("-o", "--output", help="name of output figure ", default="gamma.png")
parser.add_argument("-t", "--title", help="title of output figure ", default=None)
args = parser.parse_args()

with open('commandline_args.txt', 'a') as f:
    f.write(' '.join(sys.argv))
    f.write('\n')

scale = 0.4
plt.rcParams['figure.figsize'] = 0.9*np.array([scale*16.18033*3, scale*10/0.8])
# pre = "/Users/weilu/Research/server/april_2019/optimization_with_frag_iter2/gammas/"
# location=pre + "proteins_name_list_phi_pairwise_contact_well4.5_6.5_5.0_10phi_density_mediated_contact_well6.5_9.5_5.0_10_2.6_7.0phi_burial_well4.0_gamma_filtered"
# gamma = np.loadtxt(location)
gamma = np.loadtxt(args.location)
plot_contact_well_all(gamma, inferBound=True, invert_sign=True, vmin=-2, vmax=2,title=args.title)
# plot_contact_well_all(gamma, inferBound=False, invert_sign=True, vmin=-1, vmax=1,title=args.title)
# plot_contact_well_all(gamma, inferBound=False, invert_sign=True, vmin=-3, vmax=3,title=args.title)
# plot_contact_well_all(gamma, inferBound=False, invert_sign=True, vmin=-0.2, vmax=0.2,title=args.title)
plt.savefig(args.output)
plt.show()