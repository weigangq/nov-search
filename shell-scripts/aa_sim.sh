#!/bin/bash
for alg in {1..3}; do
	for trial in {1..50}; do
		python ../python-scripts/aa_search.py -t alg$alg-trial$trial -land ../data/gb1_fitness.tsv -gen 100 -alg $alg
	done
done
