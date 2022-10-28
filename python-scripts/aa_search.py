#!/usr/bin/env python
from aa_sim_population import *
import argparse
import logging
import sys
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(
    description='Simulated adaptive walks on a fitness landscape of peptides. Haplotypes should be combinatorially '
                'complete, although not necessarily for all 20 AAs. For example, 1st site could have 2 possible AAs, '
                '2nd site could have 5 possible AAs, then the landscape file should contain ALL 2 x 5 = 10 haplotypes. '
                'Error otherwise. Authors: Winston Koh (Qiu Lab)'
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

parser.add_argument('-mut', '--mutation_rate', type=float, default=1,
                    help='Number of mutations per generation per haplotype, randomly selected by poisson distribution. '
                         'Default = 1.')

parser.add_argument('-alg', '--algorithm', choices=['1', '2', '3'], default='1',
                    help='Evolutionary search algorithms. 1. Objective search: Selects the most fit individuals. '
                         '2. Novelty search: Selects the most novel individuals. '
                         '3. Combo search: Combines fitness and novelty.')

parser.add_argument('-beh', '--behavior', choices=['all', 'all_norm', 'blosum'], default='all',
                    help='Behavior measure of haplotypes: raw: raw distance (Euclidean pol, hydro, and iso) to peak '
                         'norm: normalized distance to peak (pol, hydro, and iso diffs equally weighted to 0-1)')

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

parser.add_argument('-out', action = 'store_true', help='Export landscape file')

args = parser.parse_args()
tagRun = args.tag
logging.basicConfig(level=logging.DEBUG)

df = pd.read_csv(args.landscape_file, sep="\t")
recs = df.to_dict('records')
#print(recs[0:10])
pep_len = len(df.at[0, 'Variants'])

# Find the variant with the highest fitness.
fitness_peak = df[df['Fitness'] == df['Fitness'].max()]
fitness_peak = fitness_peak.reset_index()
fitness_peak_value = fitness_peak.at[0, 'Fitness']
fitness_peak_hap = fitness_peak.at[0, 'Variants']
logging.info("fitness peak:\t%s\t%s\tbehavior:%s", fitness_peak_hap, round(fitness_peak_value,6), args.behavior)

# Get aa sets for each position
aa_states = {}
for pos in range(pep_len):
    aa_states[pos] = list(np.unique([x[pos] for x in df['Variants']]))

df.set_index('Variants', inplace=True)
#print(df.head)

# all or all_norm works both well for novelty
# blosum not good
def add_behavior_raw(seq):
    dist = euclidean_distance(seq, fitness_peak_hap)
    return round(dist['all'],4)

def add_behavior_norm(seq):
    dist = euclidean_distance(seq, fitness_peak_hap)
    return round(dist['all_norm'],4)

def add_behavior_pol(seq):
    dist = euclidean_distance(seq, fitness_peak_hap)
    return round(dist['pol'],4)

def add_behavior_pol_norm(seq):
    dist = euclidean_distance(seq, fitness_peak_hap)
    return round(dist['pol_norm'],4)

def add_behavior_blosum(seq):
    dist = euclidean_distance(seq, fitness_peak_hap)
    return round(dist['blosum'],4)

# write out behavior scores:
if args.out:
    for rec in recs:
        rec['beh_raw'] = add_behavior_raw(rec['Variants'])
        rec['beh_norm'] = add_behavior_norm(rec['Variants'])
        rec['beh_pol'] = add_behavior_pol(rec['Variants'])
        rec['beh_pol_norm'] = add_behavior_pol_norm(rec['Variants'])
        rec['beh_blosum'] = add_behavior_blosum(rec['Variants'])

    dfOut = pd.DataFrame.from_dict(recs)
    outfile = args.landscape_file + ".land"
    dfOut.to_csv(outfile, sep = "\t")
    sys.exit()

p = Population(pop_size=args.pop_size,
               pep_len=pep_len,
               aa_sets=aa_states,
               land=df,
               rng_seed=args.rng_seed,
               peak = fitness_peak_hap
               )

# Fitness file (elite.tsv): each generation's 10 highest fitness strings and fitness.
#elite = open(f'{tagRun}-elite.tsv', 'w')
#elite.write('tag\tgen\telite_hap\telite_fitness\talgorithm\n')
print('tag\tgen\telite_hap\telite_fitness\talgorithm')

for n in range(args.generation):
    print(f"{tagRun}\t{p.generation}\t{p.elite1[1]}\t{round(p.elite1[2],6)}\t{args.algorithm}")

    # elite.write(f"{tagRun}\t{p.generation}\t{p.elite1[1]}\t{p.elite1[2]}\t{args.algorithm}\n")

    # End the simulation when a sequence reaches the peak.
    #if p.elite1[2] == fitness_peak_value:
    if p.elite1[1] == fitness_peak_hap:
        break

    p.mutate(mut_rate=args.mutation_rate)

    if args.algorithm == '1':
        p.objective_selection()

    if args.algorithm == '2':
        p.novelty_selection(nearest_neighbors=args.nearest_neighbors, archive_method=args.archive_method,
                            prob=args.prob_arch, behave=args.behavior)

    if args.algorithm == '3':
        p.combo_selection(weight=args.weight, nearest_neighbors=args.nearest_neighbors,
                          archive_method=args.archive_method, prob=args.prob_arch, behave = args.behavior)

logging.info("Done")
#elite.close()
sys.exit()
