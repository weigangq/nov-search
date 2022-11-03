#!/bin/bash

# get landscape files:
# for k in {0..14}; do python ../../python-scripts/nk_landscape.py -K $k -N 15 > nk15-${k}.tsv; done

# remove .search & .log files to run
exe_dir="/mnt/c/Users/weiga/Dropbox/nov-search/python-scripts"
out_file=$(basename $1 .tsv)
land_file=$1
alg=2
#for alg in {1}; do
	for trial in {1..50}; do
		start=$(date +%s.%N)
		#touch ${out_file}.log
		echo -ne "$land_file\t$alg\t$trial\t..." 
		python ${exe_dir}/binary_search.py -t alg$alg-trial$trial -land $land_file -alg $alg 2> /dev/null >> ${out_file}.${alg}.search
		echo "done"
#		end=$(date +%s.%N)
#		runtime=$( echo "$end - $start" | bc -l )
#		echo $runtime >> ${out_file}.log
#	done
done
