#!/bin/bash

RUNDIR=/afs/cern.ch/work/d/devdatta/Analysis/FourthGen/CMSSW_5_3_6_patch1/src/MultileptonBprime/OS2LSelector/test/NTuples/SkimmedBkg/
BATCHDIR=${PWD}

clean_up () {
	#try to recover log files and root files
	echo try to recover log files and root files ...
	cp -p *.root $RUNDIR
	cp -p *.log  $RUNDIR
	cp -p *.out  $RUNDIR
	exit
}
#LSF signals according to http://batch.web.cern.ch/batch/lsf-return-codes.html
trap clean_up HUP INT TERM SEGV USR2 XCPU XFSZ IO

cd ${RUNDIR}
echo Setting up ${PWD} as CMSSW environment. 
eval `scram runtime -sh`
cp ${RUNDIR}/INFILE.in $BATCHDIR 

cd $BATCHDIR
echo The running directory is ${PWD}.
time OS2LSel < INFILE.in > INFILE.out 
cp -pu *.root *.out ${RUNDIR}

