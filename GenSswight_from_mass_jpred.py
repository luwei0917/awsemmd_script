#!/usr/bin/env python3
import sys

if len(sys.argv) != 3:
    print("\n" + sys.argv[0] + " inpute_file output_file\n")
    exit()

input_file = sys.argv[1]
output_file = sys.argv[2]

with open(output_file, "w") as out:
	with open(input_file, "r") as f:
		next(f)
		for line in f:
			name, data = line.rstrip().split(":")
			if name == "jnetpred":
				data = data.split(",")
				print(name, data)
				for s in data[:-1]:
					a = 0.0
					b = 0.0
					if s == 'H':
						a = 1.0
					elif s == 'E':
						b = 1.0
					out.write(str(round(a, 1)))
					out.write(' ')
					out.write(str(round(b, 1)))
					out.write('\n')
