import argparse
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(
    description='Simulated adaptive walks on a fitness landscape of peptides. Haplotypes should be combinatorially '
                'complete, although not necessarily for all 20 AAs. For example, 1st site could have 2 possible AAs, '
                '2nd site could have 5 possible AAs, then the landscape file should contain ALL 2 x 5 = 10 haplotypes. '
                'Error otherwise. Authors: Winston Koh (Qiu Lab)'
)

parser.add_argument('-c', '--cython', action='store_true',
                    help='Use the Cython version to run novelty search faster. Requires a C compiler to be installed '
                         'separately.')

parser.add_argument('-t', '--tag', default='test',
                    help='Prefix for output files. Default: "test"')

parser.add_argument('-land', '--landscape_file', required=True,
                    help='landscape file, containing tab-sep 2 columns, named as "Variants" and "Fitness". Required')

parser.add_argument('-s', '--rng_seed', default=None, type=int,
                    help='RNG seed to generate reproducible rng results. Default = None (i.e., unpredictable rng.')

parser.add_argument('-pop', '--pop_size', type=int, default=100,
                    help='Number of individuals in the population. Default = 100.')

parser.add_argument('-gen', '--generation', type=int, default=100,
                    help='Number of generations to evolve. Default = 100.')

parser.add_argument('-mut', '--mutation_rate', type=float, default=1,
                    help='Average number of mutations per generation per haplotype, randomly selected by poisson '
                         'distribution. Default = 1.')

parser.add_argument('-alg', '--algorithm', choices=['1', '2', '3'], default='1',
                    help='Evolutionary search algorithms. 1. Objective search: Selects the most fit individuals. '
                         '2. Novelty search: Selects the most novel individuals. '
                         '3. Combo search: Combines fitness and novelty.')

parser.add_argument('-beh', '--behavior', default='all',
                    choices=['pol', 'hydro', 'iso', 'all', 'blosum', 'pol_norm','hydro_norm', 'iso_norm', 'all_norm'], 
                    help='Behavior measure of haplotypes to use for novelty search. '
                         'Options: Polarity, Hydropathy index, Isoelectric index, Euclidean distance between all 3,'
                         'Normalized (0 to 1) versions of polarity, hydropathy, isoelectric, and all, '
                         'and the blosum matrix.')

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

if args.cython:
    import pyximport; pyximport.install(language_level=3)
    from aa_sim_population_cython import *

else:
    from aa_sim_population import *

df = pd.read_csv(args.landscape_file, sep="\t")
pep_len = len(df.at[0, 'Variants'])

# Find the variant with the highest fitness.
fitness_peak = df[df['Fitness'] == df['Fitness'].max()]
fitness_peak = fitness_peak.reset_index()
fitness_peak_value = fitness_peak.at[0, 'Fitness']
fitness_peak_hap = fitness_peak.at[0, 'Variants']

# Get aa sets for each position
aa_states = {}
for pos in range(pep_len):
    aa_states[pos] = list(np.unique([x[pos] for x in df['Variants']]))

df.set_index('Variants', inplace=True)

p = Population(pop_size=args.pop_size,
               pep_len=pep_len,
               aa_sets=aa_states,
               land=df,
               rng_seed=args.rng_seed
               )

print('tag\tgen\telite_hap\telite_fitness\talgorithm')
for n in range(args.generation):
    print(f"{args.tag}\t{p.generation}\t{p.elite1[1]}\t{round(p.elite1[2],6)}\t{args.algorithm}")

    # End the simulation when a sequence reaches the peak.
    if p.elite1[1] == fitness_peak_hap:
        break

    p.mutate(mut_rate=args.mutation_rate)

    if args.algorithm == '1':
        p.objective_selection()

    elif args.algorithm == '2':
        p.novelty_selection(nearest_neighbors=args.nearest_neighbors, archive_method=args.archive_method,
                            prob=args.prob_arch, behave=args.behavior)

    elif args.algorithm == '3':
        p.combo_selection(weight=args.weight, nearest_neighbors=args.nearest_neighbors,
                          archive_method=args.archive_method, prob=args.prob_arch, behave=args.behavior)
