#!/bin/bash

# Setup CMSSW environment
CMSSWPATH="/afs/cern.ch/work/d/devdatta/CMSREL/CMSSW_5_3_2_patch4/src/"
cd ${CMSSWPATH}
eval `scramv1 runtime -csh` 
cd - 

SIZE=8

for (( i=0; i < $SIZE; i++ ))
do
        newmass=`echo "450.0+${i}*50.0" | bc -l | cut -c -5`
        newwidth=`echo "${newmass}/100" | bc -l | cut -c -5`
        newbwshapelimit=`echo "${newmass}/10" | bc -l| cut -c -5`

        masspoint=${newmass/.*}

	lheid=`echo "6145.0+${i}*2.0" | bc -l | cut -c -5`
        newlheid=${lheid/.*} 

echo 
echo $newmass $newwidth $newbwshapelimit $masspoint $newlheid

newmass="$newmass"
newwidth="$newwidth"
newbwshapelimit="$newbwshapelimit"
masspoint="$masspoint"
newlheid="$newlheid"

echo $newmass $newwidth $newbwshapelimit $masspoint $newlheid 

#FASTSIM
cmsDriver.py Configuration/GenProduction/python/EightTeV/Bprime/BprimeBprimeToBHBH_HToBB_M_${masspoint}_TuneZ2star_8TeV_madgraph_cff.py --filein lhe:${newlheid} --step GEN,FASTSIM,HLT:7E33v2 --conditions START53_V7A::All --pileup 2012_Startup_inTimeOnly --datamix NODATAMIXER --beamspot Realistic8TeVCollision --eventcontent AODSIM --datatier AODSIM --no_exec -n -1

#--pileup 2012_Startup_inTimeOnly 
#--pileup NoPileUp 

sed 's/2012_Startup_inTimeOnly/2012_Summer_inTimeOnly/g' BprimeBprimeToBHBH_HToBB_M_${masspoint}_TuneZ2star_8TeV_madgraph_cff_py_GEN_FASTSIM_HLT_PU.py > BprimeBprimeToBHBH_HToBB_M_${masspoint}_TuneZ2star_8TeV_madgraph_cff_py_GEN_FASTSIM_HLT_PUSummer12.py
#
#sed 's/MASS/'"${newmass}"'D0/g;s/WIDTH/'"${newwidth}"'D0/g;s/BWSHAPELIMIT/'"${newbwshapelimit}"'D0/g;' BprimeBprimeToBHBZinc_HToBB_M_Mass_TuneZ2star_8TeV_madgraph_cff.py > BprimeBprimeToBHBZinc_HToBB_M_${masspoint}_TuneZ2star_8TeV_madgraph_cff.py
#
#sed 's/MASS/'"${newmass}"'D0/g;s/WIDTH/'"${newwidth}"'D0/g;s/BWSHAPELIMIT/'"${newbwshapelimit}"'D0/g;' BprimeBprimeToBHTWinc_HToBB_M_Mass_TuneZ2star_8TeV_madgraph_cff.py > BprimeBprimeToBHTWinc_HToBB_M_${masspoint}_TuneZ2star_8TeV_madgraph_cff.py
#
#sed 's/MASS/'"${newmass}"'D0/g;s/WIDTH/'"${newwidth}"'D0/g;s/BWSHAPELIMIT/'"${newbwshapelimit}"'D0/g;' BprimeBprimeToBHBZTWinc_M_Mass_TuneZ2star_8TeV_madgraph_cff.py > BprimeBprimeToBHBZTWinc_M_${masspoint}_TuneZ2star_8TeV_madgraph_cff.py
#
#sed 's/MASS/'"${newmass}"'D0/g;s/WIDTH/'"${newwidth}"'D0/g;s/BWSHAPELIMIT/'"${newbwshapelimit}"'D0/g;' TprimeTprimeToTHTZBWinc_M_Mass_TuneZ2star_8TeV_madgraph_cff.py > TprimeTprimeToTHTZBWinc_M_${masspoint}_TuneZ2star_8TeV_madgraph_cff.py

done

