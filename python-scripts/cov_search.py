#!/usr/bin/env python
from cov_sim_population import *
import argparse
import logging
import sys
import numpy as np
import pandas as pd
import re

parser = argparse.ArgumentParser(
    description='Simulated adaptive walks on a fitness landscape of peptides. Haplotypes may not be combinatorically complete. In the case of RBD mutation scan by Starr et al[2022], fitness values for only-single mutation haplotypes are available'
)

parser.add_argument('-t', '--tag', default='test',
                    help='Prefix for output files. Default: "test"')

parser.add_argument('-land', '--landscape_file', required=True,
                    help='landscape file, containing tab-sep 2 columns, named as "Variants" and "Fitness". Required')

parser.add_argument('-s', '--rng_seed', default=None,
                    help='RNG seed to generate reproducible rng results. Default = None (i.e., unpredictable rng.')

parser.add_argument('-pop', '--pop_size', type=int, default=100,
                    help='Number of individuals in the population. Default = 100.')

parser.add_argument('-gen', '--generation', type=int, default=100,
                    help='Number of generations to evolve. Default = 100.')

parser.add_argument('-alg', '--algorithm', choices=['1', '2', '3'], default='1',
                    help='Evolutionary search algorithms. 1. Objective search: Selects the most fit individuals. '
                         '2. Novelty search: Selects the most novel individuals. '
                         '3. Combo search: Combines fitness and novelty.')

parser.add_argument('-n', '--nearest_neighbors', type=int, default=10,
                    help='Number of nearest neighbors to use in novelty search. Default = 10.')

parser.add_argument('-arch', '--archive_method', choices=['1', '2'], default='1',
                    help='Method to use for selecting what to add to the archive. '
                         '1. Fixed size, random archive. 2. Adaptive threshold. Default = 1.')

parser.add_argument('-prob', '--prob_arch', type=float, default=0.1,
                    help='Probability for an individual to get randomly added to the archive during novelty '
                         'search method 1. Default = 0.1.')

parser.add_argument('-w', '--weight', type=float, default=0.5,
                    help='Weight for the combo search algorithm. Weight value is between 0 and 1. '
                         'Closer to 1 means more bias towards novelty. Default = 0.5.')

args = parser.parse_args()
tagRun = args.tag
logging.basicConfig(level=logging.DEBUG)

df = pd.read_csv(args.landscape_file, sep="\t")
wtIDs = df[df['Fitness']==0]['Variant'] # all wildtypes
mtIDs = df[df['Fitness']!=0]['Variant'] # all mutations
wt_seq = "".join([x[0] for x in wtIDs ])
pep_len = len(wt_seq)
pep_pos = [ int(x[1:]) for x in wtIDs ]
wt_AAs = [x for x in wtIDs ]
mt_AAs = [x for x in mtIDs ]
#print(p.elite)
#sys.exit()
# Find the variant with the highest fitness.
fitness_peak = df[df['Fitness'] == df['Fitness'].max()]
fitness_peak = fitness_peak.reset_index()
fitness_peak_value = fitness_peak.at[0, 'Fitness']
fitness_peak_hap = fitness_peak.at[0, 'Variant']
logging.info("fitness peak\t%s\t%s", fitness_peak_hap, fitness_peak_value)
#sys.exit()
df.set_index('Variant', inplace=True)

p = Population(pop_size=args.pop_size,
               wt = wt_AAs,
               mut = mt_AAs,
               land=df,
               rng_seed=args.rng_seed,
               peak = fitness_peak_hap
               )
#print(p.population)
#sys.exit()

# Fitness file (elite.tsv): each generation's 10 highest fitness strings and fitness.
#elite = open(f'{tagRun}-elite.tsv', 'w')
#elite.write('tag\tgen\telite_hap\telite_fitness\talgorithm\n')
print('tag\tgen\telite_hap\telite_fitness\tavg_fit\talgorithm')

for n in range(args.generation):
    avg_fit = round(np.average([ x[2] for x in p.elite ]),6)
    print(f"{tagRun}\t{p.generation}\t{p.elite1[1]}\t{p.elite1[2]}\t{avg_fit}\t{args.algorithm}")

    # elite.write(f"{tagRun}\t{p.generation}\t{p.elite1[1]}\t{p.elite1[2]}\t{args.algorithm}\n")

    # End the simulation when a sequence reaches the peak.
    if p.elite1[1] == fitness_peak_hap:
        break

    p.mutate()

    if args.algorithm == '1':
        p.objective_selection()

    if args.algorithm == '2':
        p.novelty_selection(nearest_neighbors=args.nearest_neighbors, archive_method=args.archive_method,
                            prob=args.prob_arch)

    if args.algorithm == '3':
        p.combo_selection(weight=args.weight, nearest_neighbors=args.nearest_neighbors,
                          archive_method=args.archive_method, prob=args.prob_arch)

logging.info("Done")
#elite.close()
sys.exit()
