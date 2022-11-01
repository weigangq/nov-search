#!/bin/bash
for alg in {1..3}; do
	for trial in {1..20}; do
		python aa_run.py -alg $alg -t alg$alg-trial$trial -gen 100
	done
done
