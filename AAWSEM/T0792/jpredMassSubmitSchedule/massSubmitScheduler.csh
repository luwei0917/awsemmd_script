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
echo "\n\nJPred REST API: mass submission and scheduling script (v.1)"
echo "Author: Dr. Alexey Drozdetskiy"
echo "Email: a.drozdetskiy@dundee.ac.uk"
echo "The Barton Group, Division of Computational Biology"
echo "College of Life Sciences, University of Dundee"
echo "www.compbio.dundee.ac.uk"
echo "${DG}"
echo "Please contact me with suggestions to improve/fix client/server\npart of the JPred API and the mass submission/scheduling scripts."
echo "\nUsage:"
echo "${BLACK}"
echo "\t./massSubmitScheduler.csh directoryname"
echo "${DG}"
echo "\nwhere directoryname - is the name of the local directory with single sequence\nFASTA format files you would want to submit as individual jobs."
echo ""
echo "NOTE: multiple controlled individual sequence submission through\nthe JPred REST API is a preferred way over single batch submission.\nFor details - please see tutorials on mass submission (and monitoring) at:\nhttp://www.compbio.dundee.ac.uk/jpred4/help.shtml\n\n"
echo "${BLACK}"
echo "${NC}"

set maxNJobs=3 # IMPORTANT: for tests - use maxNJobs between 3 and
	       # 10. This would allow you to check all the
	       # functionality, you will see if there are any errors,
	       # etc. Once tested, use maxNJobs between 10 and 30 for
	       # controlled mass submission. DO NOT go above 40!!
	       # Script/JPred would allow that, BUT it would mean that
	       # other users would experience longer waiting
	       # times... Note that this script allows you fully
	       # automatic submission, status check, results
	       # retrieval. So, just allow it to run in the background
	       # for a bit longer, but in a more mindful towards other
	       # users way. :)

set pause=120 # this is interval in seconds between job status
	      # checks. JPred job execution time is about 5 min
	      # (median) according to our estimations in Fall
	      # 2014. So, checking much more often than every 120 sec
	      # (2 min) will not help much.

# --------------------------- IMPORTANT: ---------------------------

# Only the above two parameters (actually, even the only one, the
# maxNJobs parameter) worth changing. Change the rest of the script
# below only if you know what you are doing. Code is delibirately
# written in as simple and transparent fashion as possible. However,
# if you never did shell scripting and don't care to learn it now BUT
# need some modifications, to help cater to your specific case -
# just drop me a line and I would be happy to help.


set dir=$1
rm -rf ${dir}_output ${dir}_error
mkdir ${dir}_output
mkdir ${dir}_error

set zero = 0
set counter=0
set nJobsRunning = 0
foreach line (`ls $dir/*.fasta`)
    @ counter = $counter + 1
    sleep 2

    set lengthCheck = `cat $line | sed '1d' | wc -m` 
    if ($lengthCheck<20 || $lengthCheck>800) then
        mv $line ${dir}_error/.
	echo "Length of sequence is outside of the allowed interval [20,800]: " $lengthCheck > $line.log_running_job_id
	mv $line.* ${dir}_error/.
    else
        while ($nJobsRunning >= $maxNJobs)
	    #---------------------------------------------------------
	    sleep $pause
	    foreach lineid (`ls $dir/*.log_running_job_id`)
		set id = `cat $lineid | awk '{print $6}'`
		set name = `echo $lineid | sed 's%/% %' | awk '{print $2}' | sed 's%.fasta.log_running_job_id%%g' | awk '{print $1}'`
		echo "Going to check status of jobid: " $id "name: " $name
		perl jpredapi status jobid=$id getResults=no checkEvery=once silent > $lineid.log_status
		set mystatus2 = `grep "ERROR" $lineid.log_status | wc -m`
		if ($mystatus2 > $zero) then 
		    echo "Jod $id error-ed. Movig to error directory."
		    mv $dir/$name.fasta ${dir}_error/.
		    mv -f $lineid* ${dir}_error/.
		else
		    set mystatus = `grep "finished. Results available at the following" $lineid.log_status | wc -m`
		    if ($mystatus > $zero) then
			echo "Jod $id finished. Getting results, creating job directory and moving there all the relevant files."
			set myurl = `grep "www.compbio.dundee.ac.uk" $lineid.log_status`
			set myurlJnet = `echo $myurl | sed 's%results.html%jnet%'`
			wget -o ${dir}_output/${name}_download.log $myurlJnet
			mv $id.jnet ${dir}_output/$name.jnet
			mv $dir/$name.fasta ${dir}_output/.
			mv $lineid* ${dir}_output/.
		    endif
		endif
	    end
	    #---------------------------------------------------------
	    set nJobsRunning = `ls $dir/*log_running_job_id | wc -l`
	    echo "n jobs running now: " $nJobsRunning
	end
	echo "Input file: " $line "submitting..."
	perl jpredapi submit mode=single format=fasta silent longtime=on file=$line > $line.log_running_job_id
    endif
    set nJobsRunning = `ls $dir/*log_running_job_id | wc -l`
    echo "n jobs running now: " $nJobsRunning
end

set moreLogs = 1
while ($moreLogs > 0) 
	    #---------------------------------------------------------
	    sleep $pause
	    foreach lineid (`ls $dir/*.log_running_job_id`)
		set id = `cat $lineid | awk '{print $6}'`
		echo "Going to check status of jobid: " $id
		set name = `echo $lineid | sed 's%/% %' | awk '{print $2}' | sed 's%.fasta.log_running_job_id%%g' | awk '{print $1}'`
		perl jpredapi status jobid=$id getResults=no checkEvery=once silent > $lineid.log_status
		set mystatus2 = `grep "ERROR" $lineid.log_status | wc -m`
		if ($mystatus2 > $zero) then 
		    echo "Jod $id error-ed. Movig to error directory."
		    mv $dir/$name.fasta ${dir}_error/.
		    mv -f $lineid* ${dir}_error/.
		else
		    set mystatus = `grep "finished. Results available at the following" $lineid.log_status | wc -m`
		    if ($mystatus > $zero) then
			echo "Jod $id finished. Getting results, creating job directory and moving there all the relevant files."
			set myurl = `grep "www.compbio.dundee.ac.uk" $lineid.log_status`
			set myurlJnet = `echo $myurl | sed 's%results.html%jnet%'`
			wget -o ${dir}_output/${name}_download.log $myurlJnet
			mv $id.jnet ${dir}_output/$name.jnet
			mv $dir/$name.fasta ${dir}_output/.
			mv $lineid* ${dir}_output/.
		    endif
		endif
	    end
	    #---------------------------------------------------------
    set moreLogs = `ls $dir/*.log_running_job_id | wc -l`
end

