#!/bin/csh

# Author: Alexey Drozdetskiy 
# Email: a.drozdetskiy@dundee.ac.uk
# ---
# Dr. Alexey Drozdetskiy
# Senior Postdoctoral Research Assistant
# The Barton Group
# Division of Computational Biology
# College of Life Sciences
# University of Dundee, DD1 5EH, Dundee, Scotland, UK.
# Tel:+44 1382 88731
# www.compbio.dundee.ac.uk
# The University of Dundee is registered Scottish charity: No.SC015096

# Please contact me with suggestions to improve/fix client/server part of the JPred API and the mass submission/scheduling scripts. 
# Do not assume it's OK to modify the API and to try your version on live JPred installation. 
# Users abusing JPred servers effectively preventing others from using the service will be banned.

# For documentation: run jpredapi without arguments like so: 'perl jpredapi'
# JPred API tutorial is available at: http://www.compbio.dundee.ac.uk/jpred4/api.shtml

set BLACK='\033[0;30m'
set RED='\033[0;31m'
set BLUE='\033[0;34m'
set NC='\033[0m' # No Color
set DG='\033[1;30m'

echo "${BLUE}"
echo "\n\nJPred REST API: prepare input FASTA files script (v.1)"
echo "Author: Dr. Alexey Drozdetskiy"
echo "Email: a.drozdetskiy@dundee.ac.uk"
echo "The Barton Group, Division of Computational Biology"
echo "College of Life Sciences, University of Dundee"
echo "www.compbio.dundee.ac.uk"
echo "${DG}"
echo "Please contact me with suggestions to improve/fix client/server\npart of the JPred API and the mass submission/scheduling scripts."
echo "\nUsage:"
echo "${BLACK}"
echo "\t./prepareInputs.csh filename"
echo "${DG}"
echo "\nwhere filename - is the name of the local file with multiple\nFASTA format entries you would want to submit as individual sequences."
echo ""
echo "NOTE: multiple controlled individual sequence submission through\nthe JPred REST API is a preferred way over single batch submission.\nFor details - please see tutorials on mass submission (and monitoring) at:\nhttp://www.compbio.dundee.ac.uk/jpred4/help.shtml\n\n"
echo "${BLACK}"

cat $1 | sed 's%>.*%|%g' | sed ':a;N;$\!ba;s/\n//g' | sed 's%|%\n%g' | sed '1d' > $1.seqs
echo "" >> $1.seqs
cat $1 | grep ">" | awk '{print $1}' > $1.names

set inputdir = "$1_dir"
rm -rf $inputdir
mkdir $inputdir
echo "Creating individual inputs in the directory: "$inputdir"...................."
set counter=0
foreach line (`cat $1.seqs`)
    @ counter = $counter + 1
    echo $line > $inputdir/$1_${counter}.seq
end
set counter=0
foreach line (`cat $1.names`)
    @ counter = $counter + 1
    set name = `echo $line | sed 's%>%%'`
    echo $line > $inputdir/${name}_${counter}.fasta
    cat $inputdir/$1_${counter}.seq >> $inputdir/${name}_${counter}.fasta
    rm $inputdir/$1_${counter}.seq
end
rm -f $1.seqs $1.names


echo "\nCOMPLETED.\n\n"
echo "${NC}"
