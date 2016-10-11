bg_color white

#cd ~/Dropbox/IAAWSEM/Figure/Figure3/TOP7/best_2nd
python
for i in range(1,21):
	cmd.load(str(i)+".pdb")
python end

python
import sys

for i in range(1,21):
  for j in range(1,21):
    result = cmd.cealign(str(i), str(j))
    #print result.keys()
    #print
    sys.stdout.write("%.3f " % float(result['RMSD']))
  print ""
    #print "%.3f %s " % (float(result['RMSD']),result['alignment_length'])

python end

quit
