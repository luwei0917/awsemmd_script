import sys

if len(sys.argv)!=4:
    print "\n"+sys.argv[0]+" inputpdb.pdb outputpdb.pdb frame_number\n"
    exit()

inputfile=sys.argv[1]
outputfile=sys.argv[2]
fo=open(inputfile,'r')
table=open(outputfile,'w');
time=0;
a=int(sys.argv[3]);
for line in fo:
    strs=line.split()
    if strs[0]=='END':
        time=time+1;
    if time==a and strs[0]=='ATOM':
        table.write(line);
    if time==a+1:
	break
table.write("END\n");
table.close()
fo.close()


