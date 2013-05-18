#!/bin/bash

for fin in *.in 
do
	echo "${fin}" 
	fout=`echo "${fin}" | sed 's:\(.*\)\.in:\1:g'`
	batch=`echo "${fin}" | sed 's:\(.*\)\.in:\1\.sh:g'`
	sed 's/INFILE/'"${fout}"'/g;' batch0.sh > "${batch}" 
	chmod +x ${batch} 
done

