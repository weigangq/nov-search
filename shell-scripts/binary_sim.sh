#!/bin/bash
for alg in {1..3}; do
	for trial in {1..50}; do
		python ../python-scripts/binary_search.py -t alg$alg-trial$trial -land land.tsv -gen 100 -alg $alg -pop 100
	done
done
