import subprocess
import argparse
import os
from pathlib import Path
import numpy as np
# from run_parameter import *
parser = argparse.ArgumentParser(
    description="Generate hybrid memory form HE.mem and HO.mem")
parser.add_argument("start", help="starting res id")
parser.add_argument("end", help="ending res id")
# parser.add_argument("-m", "--mode", help="choose mode", type=int, default=1)
args = parser.parse_args()
# sed -i.bak 's/\/work\/pw8\/mc70\/script/\/home\/wl45/g' frag_HE.mem

do = os.system

def file_len(fname):
    p = subprocess.Popen(['wc', '-l', fname], stdout=subprocess.PIPE,
                                              stderr=subprocess.PIPE)
    result, err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)
    return int(result.strip().split()[0])


def file_width(fname):
    p = subprocess.Popen(['wc', '-c', fname], stdout=subprocess.PIPE,
                                              stderr=subprocess.PIPE)
    result, err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)
    return int(result.strip().split()[0])-1

protein_size = 0
end = int(args.end)
start = int(args.start)
def backup():
    my_file = Path("origin_ssweight")
    if my_file.is_file():
        print("backup exist, won't overide")
    else:
        # print("hi")
        do("cp ssweight origin_ssweight")
        do("cp Hybrid.mem origin_Hybrid.mem")
        do("cp go_rnativeCACA.dat origin_go_rnativeCACA.dat")
        do("cp go_rnativeCACB.dat origin_go_rnativeCACB.dat")
        do("cp go_rnativeCBCB.dat origin_go_rnativeCBCB.dat")
        do("cp casp.seq origin_casp.seq")

def ssweight():
    protein_size = file_len("origin_ssweight")
    print(protein_size)

    cmd = f"head -n {end} origin_ssweight | sed 1,{start}d > ssweight"
    do(cmd)
    print(cmd)


def fragMemory():
    with open('Hybrid.mem', 'w') as w:
        with open('origin_Hybrid.mem', 'r') as f:
            l1 = next(f)
            l2 = next(f)
            l3 = next(f)
            l4 = next(f)
            w.write(l1)
            w.write(l2)
            w.write(l3)
            w.write(l4)
            lines = f.readlines()
            if end < 8:
                print("end is less than 8, won't have any frag memeory")
            else:
                truncated_lines = lines[(int(start)-1)*20:(int(end)-8)*20]
                print("origin", len(lines)//20, "now", len(truncated_lines)//20)
                for line in truncated_lines:
                    path, loc, start_i, fragLens, weight = line.split()
                    # w.write(line)
                    w.write(" ".join([path, str(int(loc)-start+1), start_i, fragLens, weight])+"\n")


def er():
    data = np.loadtxt("origin_go_rnativeCACA.dat")
    np.savetxt("go_rnativeCACA.dat", data[start-1:end,start-1:end], fmt='%.3f')
    data = np.loadtxt("origin_go_rnativeCACB.dat")
    np.savetxt("go_rnativeCACB.dat", data[start-1:end,start-1:end], fmt='%.3f')
    data = np.loadtxt("origin_go_rnativeCBCB.dat")
    np.savetxt("go_rnativeCBCB.dat", data[start-1:end,start-1:end], fmt='%.3f')

def seq():
    with open("casp.seq", "w") as w:
        with open("origin_casp.seq", "r") as f:
            for line in f:
                truncated = line.strip()[start-1:end]
                print(truncated)
                print(len(truncated))
                w.write(truncated+"\n")
backup()
ssweight()
fragMemory()
er()
seq()