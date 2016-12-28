#!/usr/bin/env python3
import glob
import os
os.system("mkdir -p _data")
os.chdir("_output")
folder_list = glob.glob("*.jnet")
# print(folder_list, len(folder_list))
# folder_list = ['T0759.jnet']
ssdata = "../_data/"

for protein in folder_list:
    with open(protein) as input_data:
        protein_name = protein.split(".")[0]
        print(protein_name)
        out = open(ssdata+protein_name, 'w')
        for line in input_data:
            # print(line)
            top, *data = line.split(",")
            if(top == "jnetpred:-"):
                # print("hello")
                # print(line)
                print(data)
                a = 0.0
                b = 0.0
                out.write(str(round(a, 1)))
                out.write(' ')
                out.write(str(round(b, 1)))
                out.write('\n')
                for i in data[:-1]:     # last one is \n
                    a = 0.0
                    b = 0.0
                    if(i == "E"):
                        b = 1.0
                    elif(i == "H"):
                        a = 1.0
                    out.write(str(round(a, 1)))
                    out.write(' ')
                    out.write(str(round(b, 1)))
                    out.write('\n')

# if len(sys.argv)!=3:
# 	print "\n" + sys.argv[0] + " inpute_file output_file\n"
# 	exit()
#
# input_file = sys.argv[1]
# output_file = sys.argv[2]
#
# inp = open(input_file, 'r')
# st = inp.read().strip()
# inp.close()
#
# out =  open(output_file, 'w')
# for s in st:
# 	a = 0.0
# 	b = 0.0
# 	if s=='H':
# 		a = 1.0
# 	elif s=='E':
# 		b = 1.0
# 	out.write(str(round(a,1)))
# 	out.write(' ')
# 	out.write(str(round(b,1)))
# 	out.write('\n')
# out.close()
