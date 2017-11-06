import subprocess


def length_from_fasta(fasta):
    length = -1
    with open(fasta) as input_data:
        data = ""
        for line in input_data:
            if(line[0] == ">"):
                pass
            elif(line == "\n"):
                pass
            else:
                data += line.strip("\n")
        length = len(data)
    return length



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
