#!/bin/bash

exe_dir="/mnt/c/Users/weiga/Dropbox/nov-search/python-scripts"
land_file=$1

for alg in {1..3}; do
	for trial in {1..30}; do
		python ${exe_dir}/cov_search.py -t alg$alg-trial$trial -land $land_file -alg $alg
	done
done
