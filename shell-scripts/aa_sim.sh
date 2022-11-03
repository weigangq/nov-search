#!/bin/bash

#exe_dir="$HOMEDropbox/nov-search/python-scripts"
exe_dir=$1
land_file=$2
#behavior=$3
alg=2
out_file=$(basename $land_file .tsv)
for trial in {1..50}; do
    echo -en "$land_file\t$alg\t$trial\t..."
    python ${exe_dir}/aa_search.py -t alg$alg-trial$trial -land $land_file -alg $alg >> ${out_file}.${alg}.search
    echo "done"     
done

