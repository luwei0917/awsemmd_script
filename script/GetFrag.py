import sys
import os
import glob
import subprocess

# subprocess could only be used for python 2.7 or above.


def number_assign(strength, number):

    # strength is the list of number of memories in all, 20 memories needs to
    # be produced.
    totalnumber = sum(strength)
    print(totalnumber)
    for i in range(len(strength)):
        if strength[i] <= totalnumber / (number + 3):
            strength[i] = 0
    totalnumber = sum(strength)
    MemInd = []
    for i in range(len(strength)):
        aaa = (strength[i] * number * 1.0) / totalnumber - \
            int((strength[i] * number) / totalnumber)
        print(aaa)
        if aaa >= 0.5:
            MemInd.append(int((strength[i] * number) / totalnumber) + 1)
        if aaa < 0.5:
            MemInd.append(int((strength[i] * number) / totalnumber))
    # print MemInd;
    return MemInd


def memory_assign(number, wd):
    # number: is the number of memories to be assigned, in our case should be
    # 20
    memfile = open("frag.mem", 'w')
    memfile.write("[Target]" + "\n")
    memfile.write("query" + "\n" + "\n")

    memfile.write("[Memories]" + "\n")
    # T089
    # start1=[1,22,42,56];
    # start2=[86,107,127,141]
    # length=[41,34,28,26];

    # TOP7
    # start1=[1,9,22,40];
    # start2=[1,9,22,40];
    # length=[18,27,24,37];

    # T120
    start1 = [1, 17, 31, 41, 48, 62, 82, 92]
    start2 = [1, 17, 31, 41, 48, 62, 82, 92]
    length = [25, 22, 18, 20, 28, 29, 21, 24]

    # T251
    # start1=[1,17,40,53,71,77];
    # start2=[1,17,40,53,71,77];
    # length=[33,29,27,24,13,23];

    # 1N2X
    # start1=[1,22,31,48,64,77];
    # start2=[116,137,146,163,179,192];
    # length=[30,23,30,29,21,25];

    # TOP7
    # start1=[1,13,24,44,54,74];
    # start2=[1,13,24,44,54,74];
    # length=[22,30,29,30,30,19];

    # 1R69
    # start1=[1,17,27,44];
    # start2=[1,17,27,44];
    # length=[23,21,26,20];

    # GAGB
    # start1=[1,22];
    #     start2=[1,22];
    #     length=[39,35];

    fragnum = len(start1)
    for ii in range(1, fragnum + 1):
        # wd = "/dascratch/mc70/opt/script/1R69_AWSEM/fragment_final";
        fragwd = wd + "/frag" + str(ii)
        # qqq = subprocess.Popen(['cd',fragwd], stderr=subprocess.STDOUT, stdout=subprocess.PIPE)
        qqq = subprocess.check_output("cd %s" % (fragwd), shell=True)
        os.chdir(fragwd)
        print(fragwd)
        filelist = glob.glob("clustersize*.pdb")
        print(filelist)
        strength = []
        for i in range(len(filelist)):
            words = filelist[i].split(".")
            strength.append(int(words[1]))
        print(strength)
        MemInd = number_assign(strength, number)
        print(MemInd)
        for i in range(len(filelist)):
            if MemInd[i] != 0:
                namegro = "frag_" + str(i) + '.gro'
                # q = subprocess.Popen(['python','~/opt/script/Pdb2Gro.py',filelist[i],namegro], stderr=subprocess.STDOUT, stdout=subprocess.PIPE);
                os.chdir(fragwd)
                q = subprocess.check_output(
                    "python2 ~/opt/script/Pdb2Gro.py %s %s" % (filelist[i][0:-4], namegro), shell=True)
                for j in range(MemInd[i]):
                    memfile.write(fragwd + "/" + namegro + " " + str(start1[ii - 1]) + " " + str(
                        start2[ii - 1]) + " " + str(length[ii - 1]) + " " + "1" + "\n")
                # qq = subprocess.Popen(['cd','..'], stderr=subprocess.STDOUT, stdout=subprocess.PIPE);
        # qq = subprocess.check_output("cd ..",shell=True)

# wd = sys.argv[1]
wd = "/Users/weilu/Documents/Research/July-18/t120_aa"
memory_assign(20, wd)
