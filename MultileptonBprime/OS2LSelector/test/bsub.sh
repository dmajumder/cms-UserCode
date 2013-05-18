#!/bin/bash

SUBDIR=${PWD}

for fin in OSOF_*.sh
do
	jobfile="${SUBDIR}/${fin}"
	outfile=`echo ${fin} | sed 's:\.sh:\.out:g'`
	echo " Submitting job ${jobfile}"
	bsub -o ${outfile} -q 8nh ${jobfile} 
done
