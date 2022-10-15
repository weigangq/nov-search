#!/bin/bash

exe_dir="/mnt/c/Users/weiga/Dropbox/nov-search/python-scripts"
land_file=$1

for alg in {1..3}; do
	for trial in {1..50}; do
		python ${exe_dir}/aa_search.py -t alg$alg-trial$trial -land $land_file -gen 100 -alg $alg
	done
done
