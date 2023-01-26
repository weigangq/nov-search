from binary_sim_population import *
import argparse
import pandas as pd

# Setup: arguments, parameters, and logging
parser = argparse.ArgumentParser(
    description="Search of fitness of random binary strings. Author: Winston Koh (Qiu Lab)")

parser.add_argument('-t', '--tag', default='test',
                    help='Prefix for output files. Default: "test"')

parser.add_argument('-land', '--landscape_file', required=True,
                    help='landscape file, output from nk_landscape.py. Required')

parser.add_argument('-s', '--rng_seed', type=int, default=None,
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

df = pd.read_csv(args.landscape_file, sep="\t", dtype={'haplotype': str})
seq_len = len(df.at[0, 'haplotype'])

# Find the variant with the highest fitness.
fitness_peak = df[df['global_peak'] == 1]
fitness_peak.reset_index(inplace=True)
df.set_index('haplotype', inplace=True)

# Print simulation parameters
print(f"Simulation parameters:\n\tLandscape file: {args.landscape_file}\n\tPopulation size: {args.pop_size}\n\t"
      f"Search algorithm: {args.algorithm}")


# Initialize a population
p = Population(pop_size=args.pop_size,
               seq_len=seq_len,
               land=df,
               rng_seed=args.rng_seed
               )

print(f"Starting haplotype on landscape: {arr_to_str(p.population[0])}\n")

# Search by evolution
print(f"Tag\tGen\ttop_elite\tfit\tN\tK\tsearch_algo")

for n in range(args.generation):
    print(f"{args.tag}\t{p.generation}\t{p.elite1[1]}\t{p.elite1[2]}\t{fitness_peak.at[0, 'N']}\t"
          f"{fitness_peak.at[0, 'K']}\t{args.algorithm}")

    # end when reaches the global peak
    if p.elite1[1] == fitness_peak.at[0, 'haplotype']:
        break

    p.mutate(mut_rate=args.mutation_rate)

    if args.algorithm == '1':
        p.objective_selection()

    if args.algorithm == '2':
        p.novelty_selection(nearest_neighbors=args.nearest_neighbors, archive_method=args.archive_method,
                            prob=args.prob_arch)

    if args.algorithm == '3':
        p.combo_selection(weight=args.weight, nearest_neighbors=args.nearest_neighbors,
                          archive_method=args.archive_method, prob=args.prob_arch)

