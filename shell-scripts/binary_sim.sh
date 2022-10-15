#!/bin/bash

exe_dir="/mnt/c/Users/weiga/Dropbox/nov-search/python-scripts"
out_file=$(basename $1 .tsv)
for alg in {1..3}; do
	for trial in {1..20}; do
		python ${exe_dir}/binary_search.py -t alg$alg-trial$trial -land $land_file -gen 100 -alg $alg -pop 100 
	done
done
