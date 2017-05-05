import subprocess
import argparse
# from run_parameter import *
parser = argparse.ArgumentParser(
    description="Generate hybrid memory form HE.mem and HO.mem")
parser.add_argument("protein", help="the name of protein")
parser.add_argument("-m", "--mode", help="choose mode", type=int, default=1)
args = parser.parse_args()

try:
    protein_name, _ = args.protein.split('.')
except:
    try:
        protein_name, _ = args.protein.split('/')
    except:
        protein_name = args.protein
        print("ATTENSION! protein name {}\n Correct?".format(protein_name))


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
    return int(result.strip().split()[0])


def hybrid(protein_name):
    protein_size = file_width(protein_name+".seq")
    print(protein_size)
    n = 0
    nn = 0
    if(args.mode == 1):
        with open('Hybrid.mem', 'w') as w:
            with open('frag_HE.mem', 'r') as f:
                l1 = next(f)
                l2 = next(f)
                l3 = next(f)
                l4 = next(f)
                w.write(l1)
                w.write(l2)
                w.write(l3)
                w.write(l4)
                he_lines = f.readlines()
            # print(ha_lines)
            with open('frag_HO.mem', 'r') as f:
                next(f)
                next(f)
                next(f)
                next(f)
                ho_lines = f.readlines()
            for i in range(1, protein_size-7):
                fragGroup = []
                for line in ho_lines:
                    # _, loc, _, fragLens, _ = line.split(" ")
                    # window, _, score, evalue, name, loc, j_start, fragLens, *weight = line.split()
                    path, loc, start, fragLens, weight = line.split()
                    locEnd = int(loc) + int(fragLens)
                    # print(loc, fragLens, locEnd)
                    if(i <= int(loc) and i+9 >= locEnd):
                        fragGroup += [line]
                groupLens = len(fragGroup)
                if(groupLens != 0):
                    nn += 1
                    for j in range(20):
                        w.write(fragGroup[j % groupLens])
                else:
                    n += 1
                    for j in range(20):
                        try:
                            w.write(he_lines[(i-1)*20+j])
                        except:
                            print(i*20+j)

    if(args.mode == 2):
        with open('Hybrid.mem', 'w') as w:
            with open('HE.mem', 'r') as f:
                he_lines = f.readlines()
            # print(ha_lines)
            with open('HO.mem', 'r') as f:
                ho_lines = f.readlines()
            for i in range(1, protein_size-7):
                fragGroup = []
                for line in ho_lines:
                    # _, loc, _, fragLens, _ = line.split(" ")
                    window, _, score, evalue, name, loc, j_start, fragLens, *weight = line.split()
                    locEnd = int(loc) + int(fragLens)
                    # print(loc, fragLens, locEnd)
                    if(i <= int(loc) and i+9 >= locEnd):
                        fragGroup += [line]
                groupLens = len(fragGroup)
                if(groupLens != 0):
                    nn += 1
                    for j in range(20):
                        w.write(fragGroup[j % groupLens])
                else:
                    n += 1
                    for j in range(20):
                        try:
                            w.write(he_lines[(i-1)*20+j])
                        except:
                            print(i*20+j)

hybrid(protein_name)
