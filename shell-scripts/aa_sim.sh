#!/bin/bash

#exe_dir="$HOMEDropbox/nov-search/python-scripts"
exe_dir=$1
land_file=$2
behavior=$3
for alg in {1..3}; do
	for trial in {1..30}; do
		python ${exe_dir}/aa_search.py -t alg$alg-trial$trial -land $land_file -alg $alg -beh $behavior -gen 200 -pop 200
	done
done
