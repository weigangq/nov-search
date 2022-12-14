from aa_search import *
import argparse

parser = argparse.ArgumentParser(
    description='Simulate the evolution of Streptococcal bacteria protein G domain B1 using the fitness landscape from '
                'https://doi.org/10.7554/eLife.16965. Author: Winston Koh (Qiu Lab)'
)

parser.add_argument('-t', '--tag', default='test',
                    help='Prefix for output files. Default: "test"')

parser.add_argument('-r', '--rng_seed', default=None,
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
tagRun = args.tag

fitness_peak_value = fitness_peak.at[0, 'Fitness']
p = Population(pop_size=args.pop_size, rng_seed=args.rng_seed)

# Fitness file (elite.tsv): each generation's 10 highest fitness strings and fitness.
elite = open(f'results/{tagRun}-elite.tsv', 'w')
elite.write('tag\tgen\telite_hap\telite_fitness\talgorithm\n')
print('tag\tgen\telite_hap\telite_fitness\talgorithm')

for n in range(args.generation):
    print(f"{tagRun}\t{p.generation}\t{p.elite1[1]}\t{p.elite1[2]}\t{args.algorithm}")

    elite.write(f"{tagRun}\t{p.generation}\t{p.elite1[1]}\t{p.elite1[2]}\t{args.algorithm}\n")

    # End the simulation when a sequence reaches the peak.
    if p.elite1[2] == fitness_peak_value:
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

elite.close()
