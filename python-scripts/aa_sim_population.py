#!/usr/bin/env python
import numpy as np
import pandas as pd
import sys

################################################################
# AA distances
# List of all possible letters representing amino acids, sorted
amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

# Amino acid physicochemical properties
# Polarity (Grantham, 1974)
polarity = {'A': 8.1, 'C': 5.5, 'D': 13, 'E': 12.3, 'F': 5.2, 'G': 9, 'H': 10.4, 'I': 5.2, 'K': 11.3, 'L': 4.9,
            'M': 5.7, 'N': 11.6, 'P': 8, 'Q': 10.5, 'R': 10.5, 'S': 9.2, 'T': 8.6, 'V': 5.9, 'W': 5.4, 'Y': 6.2
            }

# Hydropathy index (Kyte-Doolittle, 1982)
hydropathy = {'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5, 'F': 2.8, 'G': -0.4, 'H': -3.2, 'I': 4.5, 'K': -3.9, 'L': 3.8,
              'M': 1.9, 'N': -3.5, 'P': -1.6, 'Q': -3.5, 'R': -4.5, 'S': -0.8, 'T': -0.7, 'V': 4.2, 'W': -0.9, 'Y': -1.3
              }

# Isoelectric index (Zimmerman, 1968)
iso = {'A': 6, 'C': 5.05, 'D': 2.77, 'E': 3.22, 'F': 5.48, 'G': 5.97, 'H': 7.59, 'I': 6.02, 'K': 9.74, 'L': 5.98,
              'M': 5.74, 'N': 5.41, 'P': 6.3, 'Q': 5.65, 'R': 10.76, 'S': 5.68, 'T': 5.66, 'V': 5.96, 'W': 5.89, 'Y': 5.66
              }

# weights for calculating novelty scores
# ACTY -> VCTY (pos: 1/4 * 0.7 + 2.1/max * 0.1 + 0, etc)
# ACTY -> GWSQ: low chem dist, but large hamming distance
dist_wts = { 'pos': 0.7, 'pol': 0.1, 'hydro': 0.1, 'iso': 0.1}

pol_dist = {}
hydro_dist = {}
iso_dist = {}
for i in range(len(amino_acids)): # sorted aa list
    for j in range(i, len(amino_acids)):
        pol_dist[(amino_acids[i], amino_acids[j])] = abs(polarity[amino_acids[i]] - polarity[amino_acids[j]])
        hydro_dist[(amino_acids[i], amino_acids[j])] = abs(hydropathy[amino_acids[i]] - hydropathy[amino_acids[j]])
        iso_dist[(amino_acids[i], amino_acids[j])] = abs(iso[amino_acids[i]] - iso[amino_acids[j]])

def normalize(arr, t_min, t_max):
    norm_arr = []
    diff = t_max - t_min
    diff_arr = max(arr) - min(arr)
    for i in arr:
        temp = (((i - min(arr))*diff)/diff_arr) + t_min
        norm_arr.append(temp)
    return norm_arr

def normalize_dict(dic, t_min, t_max):
    norm_dic = {}
    diff = t_max - t_min
    vals = list(dic.values())
    diff_arr = max(vals) - min(vals)
    for k in dic:
        temp = (((dic[k] - min(vals))*diff)/diff_arr) + t_min
        norm_dic[k] = {'val': round(dic[k],6), 'norm': round(temp,6)}
    return norm_dic

pol_norm = normalize_dict(pol_dist, 0, 1)
hyd_norm = normalize_dict(hydro_dist, 0, 1)
iso_norm = normalize_dict(iso_dist, 0, 1)


def get_distance(seq1: str, seq2: str) -> float:
    pep1 = list(seq1)
    pep2 = list(seq2)
    
    dist_pol = 0
    dist_hydro = 0
    dist_iso = 0
    hamming_dist = 0
    for i in range(len(seq1)):
        aa_pair = [pep1[i], pep2[i]]
        aa_pair.sort()
        aa1 = aa_pair[0]
        aa2 = aa_pair[1]
        key = (aa1, aa2)
        dist_pol += pol_norm[key]['norm'] 
        dist_hydro += hyd_norm[key]['norm']
        dist_iso += iso_norm[key]['norm']
        
        if pep1[i] != pep2[i]:
            hamming_dist += 1
    
    d = (dist_pol * dist_wts['pol'] + dist_hydro * dist_wts['hydro'] + dist_iso * dist_wts['iso'] + hamming_dist* dist_wts['pos'])/len(seq1)
    #print(d)
    return d

# defines behavior
def euclidean_distance(seq1: str, seq2: str) -> dict:
    pep1 = list(seq1)
    pep2 = list(seq2)
    
    dists = {}
    dists['pol'] = 0
    dists['hydro'] = 0
    dists['iso'] = 0
    dists['all'] = 0
    dists['pol_norm'] = 0
    dists['hydro_norm'] = 0
    dists['iso_norm'] = 0
    dists['all_norm'] = 0
    for i in range(len(seq1)):
        aa_pair = [pep1[i], pep2[i]]
        aa_pair.sort()
        aa1 = aa_pair[0]
        aa2 = aa_pair[1]
        key = (aa1, aa2)
        # use normalized diff to weight the three values equally
        dists['pol_norm'] +=  pol_norm[key]['norm'] 
        dists['hydro_norm'] += hyd_norm[key]['norm']
        dists['iso_norm'] += iso_norm[key]['norm']
        dists['all_norm'] += dists['pol_norm'] + dists['hydro_norm'] + dists['iso_norm'] # not good
        # not normalized:
        dists['pol'] += (polarity[aa1] - polarity[aa2]) ** 2 # good
        dists['hydro'] += (hydropathy[aa1] - hydropathy[aa2]) ** 2 
        dists['iso'] += (iso[aa1] - iso[aa2]) ** 2 
        dists['all'] += dists['pol'] + dists['hydro'] + dists['iso'] # this is better!!!!

    for key in dists:
        if key in ['pol', 'hydro', 'iso', 'all']:
            dists[key] = np.sqrt(dists[key])/len(pep1)
        else:
            dists[key] = dists[key]/len(pep1)            
    return dists
#######################################################
sequence_behavior = {}


def arr_to_str(nparray: np.ndarray) -> str:
    """
    Takes a 1D numpy array of strings and joins them into a single string.

    Args:
        nparray: numpy.ndarray with ndim = 1 and dtype = str.

    Returns:
        str
    """
    return ''.join(nparray)


def aa_fitness(aa_string: np.ndarray, land: pd.DataFrame) -> float:
    """
    Get the fitness of the amino acid string from the landscape data frame.

    Args:
        land: pandas.DataFrame
            Data frame of the fitness landscape. The landscape must contain all possible sequences and their
            corresponding fitness value.

        aa_string: numpy.ndarray with ndim = 1, dtype = str.
            Amino acid sequence represented by a 1D numpy ndarray, where each element is the 1 letter symbol of an
            amino acid.

    Returns:
        float
    """
    pep = arr_to_str(aa_string)
    try:
        fitness = land.at[pep, 'Fitness']
    except KeyError:
        print(f"Error: Haplotype '{pep}' not found in landscape. The landscape must be combinatorially complete.")
        sys.exit()

    return fitness


def population_fitness(pop: np.ndarray, land: pd.DataFrame) -> list:
    """
    Returns a list containing only the fitness values of the sequences in the population.
    Fitness values are obtained from the landscape data frame.
    Fitness values have the same index in this list as the index of the sequence in the population.

    Args:
        pop: numpy.ndarray with ndim = 2, dtype = str.
            2D numpy ndarray that contains all the amino acid sequences in a population.

        land: pandas.DataFrame
            Data frame of the fitness landscape. The landscape must contain all possible sequences and their
            corresponding fitness value.

    Returns:
        list
    """
    return [aa_fitness(pop[i], land) for i in range(pop.shape[0])]


def get_elites(pop: np.ndarray, population_metrics: list, land: pd.DataFrame) -> list:
    """
    Used by a Population object to find the top 10 sequences of the population based on fitness/novelty/combo without
    having to sort the entire list of values.

    Args:
        pop: numpy.ndarray
            Intended to be an array from a Population.population.

        population_metrics: list generated by population_fitness(), Population.population_novelty(), or
            Population.population_combo().
            This list contains only the fitness/novelty/combo scores of the sequences in the population.
            The scores have the same index in the list as the index of the corresponding sequence in the population.

        land: pandas.DataFrame
            Data frame of the fitness landscape. The landscape must contain all possible sequences and their
            corresponding fitness value.

    Returns:
        list of 10 * [Sequence's index in population, Sequence string, Fitness]
    """
    # Temporary elite list containing [Sequence's population index, Metric]
    # replace the smallest one, iteratively
    # metric is either fitness, sparsity, or a weighted sum of the two
    # regardless, get the top 10 highest
    elite = [[n, population_metrics[n]] for n in range(10)]
    elite_min = min(elite, key=lambda x: x[1])
    for index in range(10, pop.shape[0]):
        if population_metrics[index] > elite_min[1]:
            elite.remove(elite_min)
            elite.append([index, population_metrics[index]])
            elite_min = min(elite, key=lambda x: x[1])

    elite_fitness = [[e[0], arr_to_str(pop[e[0]]), aa_fitness(pop[e[0]], land)] for e in elite]
    return elite_fitness


class Population:
    def __init__(self, pop_size: int, pep_len: int, aa_sets: dict, land: pd.DataFrame, rng_seed: object, peak: str) -> None:
        """
        Creates Population object. Contains amino acid sequences, along with other information about the population.

        Args
            pop_size: int, default = 100.
                Number of sequences in the population. pop_size must be greater than 10. If pop_size is not a multiple
                of 10, it will change to become the largest multiple of 10 that is less than the given pop_size.
                After this, the population size remains constant in this simulation.

            pep_len: int
                Length of the amino acid sequences. This should be obtained from the landscape file.

            aa_sets: dict
                Dictionary containing the amino acids that have been seen at each site in the landscape.

            land: pandas.DataFrame
                Data frame of the fitness landscape. The landscape must contain all possible sequences and their
                corresponding fitness value.

            rng_seed: {None, int, array_like[ints], SeedSequence, BitGenerator, Generator}, default = None.
                To generate reproducible results, specify a seed to initialize numpy.random.default_rng().
                By default, the rng results will be unpredictable.
        """
        if pop_size < 10:
            raise ValueError('The population size should be at least 10.')
        self.pop_size = (pop_size // 10) * 10

        # Generation counter that counts the number of times the population has undergone mutations.
        self.generation = 0

        self.pep_len = pep_len

        self.land = land

        self.aa_sets = aa_sets
        self.global_peak = peak

        # Numpy rng
        self.rng = np.random.default_rng(seed=rng_seed)

        # Initialize a population of a single random amino acid sequence
        # Each amino acid has equal probability of being picked for a position
        pep_start = [self.rng.choice(aa_sets[pos]) for pos in aa_sets]
        self.population = np.broadcast_to(pep_start, (pop_size, pep_len)).copy()

        # Elite list contains each generation's top 10 sequences based on fitness, novelty, or a combination.
        # Will be updated during each selection function.
        # Format: list of 10 * [Sequence's index in population, Sequence string, Fitness]
        self.elite = get_elites(self.population, population_fitness(self.population, land), land)

        # The highest scoring sequence of the generation based on fitness, novelty, or a combination.
        # Will be updated during each selection function.
        # Format: [Sequence's index in population, Sequence string, Fitness]
        self.elite1 = max(self.elite, key=lambda x: x[2])

        # Archive contains sequences that were novel. Used during novelty and combo search.
        self.archive = np.empty((0, pep_len), dtype=str)
        
    def mutate(self, mut_rate: float = 1) -> None:
        """
        Mutate a random number of residues in each sequence in the population.
        The number is chosen by a Poisson distribution where the expected value is the mutation rate.

        Args:
            mut_rate: float, default = 1.
                Average number of mutations for a sequence each generation.
        """
        # Increase generation counter.
        self.generation += 1

        # For each sequence in the population, mutate a random number of residues.
        for n in range(self.pop_size):
            # Get the number of sites to mutate in an amino acid string.
            num_mut = self.rng.poisson(mut_rate)
            if num_mut > 0:
                # Ensure that number of mutations is at most the peptide length.
                num_mut = min(num_mut, self.pep_len)
                # Select random indices to mutate.
                mut_indices = self.rng.choice(self.pep_len, num_mut, replace=False)
                # Mutate amino acids at the chosen indices.
                for i in mut_indices:
                    # possible_aa = amino_acids.copy()
                    possible_aa = self.aa_sets[i].copy()
                    possible_aa.remove(self.population[n, i])
                    self.population[n, i] = self.rng.choice(possible_aa)
        return

    def replace_pop(self) -> None:
        """
        Replaces the current population of amino acids with a new population containing only the elite sequences.
        Takes each sequence in the elite and broadcasts it to be 1/10 of the population size.
        The elites are determined by the selection method (e.g., objective_selection) that calls this method.
        """
        # Sequences in numpy.ndarray format.
        elite_nparray = [self.population[e[0]] for e in self.elite]

        # Clear population array.
        self.population = np.empty((0, self.pep_len), dtype=str)
        # Broadcast elites to 10% of pop size, then append.
        for seq in elite_nparray:
            self.population = np.append(self.population, np.broadcast_to(seq, (self.pop_size // 10, self.pep_len)),
                                        axis=0)
        return

    def objective_selection(self) -> None:
        """
        Selects the top 10 sequences by fitness value, then replaces the population with the elites (fittest).
        """
        pop_fitness = population_fitness(self.population, self.land)
        self.elite = get_elites(self.population, pop_fitness, self.land)
        self.elite1 = max(self.elite, key=lambda x: x[2])
        self.replace_pop()
        return

    def population_novelty(self, nearest_neighbors: int = 10, archive_method: int = 1, prob: float = 0.10) -> list:
        """
        Evaluates each sequence's novelty and updates the archive. Then returns a list containing only the novelty
        values of the sequences in the population.
        Novelty values have the same index in this list as the index of the sequence in the population.

        Novelty:
        A sequence (or individual) is considered to be novel based solely on how different it is to all the other
        sequences evaluated so far. Novelty is calculated as the sequence's average distance from the k-nearest
        neighbors in both the current population and an archive of past sequences.

        Distance between sequences can be evaluated based on the distance between genotypes or the behavior of the
        sequences. Behaviors are represented using a behavior characterization, a set of M features, which can be based
        on the genotype, the phenotype, or characteristics of the problem domain. Then the distance between two
        behaviors is calculated using a distance metric (e.g., Euclidean distance, Hamming distance). The final novelty
        score is the average distance to the k-nearest neighbors in behavior space.

        For this simulation, the measured behavior is based on the physico-chemical properties of the constituent amino
        acids, specifically polarity (Grantham, 1974) and hydropathy index (Kyte-Doolittle, 1982).
        Values are taken from http://diverge.hunter.cuny.edu/~weigang/code-wheel/

        Args:
            nearest_neighbors: int, default = 10.
                The number of nearest neighbors to evaluate when calculating the novelty score.

            archive_method: int {1, 2}
                Method 1: Using fixed size random archive.
                    (Lehman and Stanley 2010, Efficiently Evolving Programs through the Search for Novelty)
                    Each individual has a 10% chance of being added to the archive, regardless of their novelty value.
                    The archive size stays fixed at 2 * pop_size. When adding a new sequence, remove the oldest.

                    Keeping the archive to a fixed size limits the maximum number of calculations needed to determine a
                    sequence's novelty.
                    Creating the archive by probabilistic sampling

                    The effect of such probabilistic sampling is that the permanent archive approximately characterizes
                    the distribution of prior solutions in behavior space without pushing search away from newly
                    discovered areas.

                Method 2: Unrestricted archive size using adaptive threshold.
                    (Kistemaker and Whiteson 2011, Critical Factors in the Performance of Novelty Search)
                    The archive only contains individuals that were considered novel when discovered. If the novelty of
                    a new sequence is higher than the threshold (ρ(x) > ρ_min), it is added to the archive.
                    The threshold's value changes in order to keep the growth of the archive approximately constant.
                    ρ_min is increased by a fixed fraction if the number of added behaviors exceeds the add_max
                    threshold in a certain number of evaluations. If the number of added behaviors is lower than add_min
                    in a certain number of evaluations, ρ_min is decreased by a fixed fraction.

            prob: float between 0 and 1. default = 0.1.
                Probability that a sequence is added to the archive in method 1. If using method 2, this value has no
                use.

        Returns:
            list
        """
        pop_novelty = []

        # Parameters for novelty search method 2
        threshold = 10
        add_max = 15
        add_min = 5
        fraction = 0.25

        num_eval = 0
        counter = 0

        # Calculate behavior of population + archive.
        global sequence_behavior
        for index in range(self.pop_size):
            # Check to see if the sequence's behavior was previously calculated. If not, then calculate the
            # sequence's polarity and hydropathy and record values.
            seq = arr_to_str(self.population[index])
            if seq in sequence_behavior:
                continue
            else:
                behavior = {}
                '''
                polarity_value, hydropathy_value, iso_value = 0, 0, 0
                for residue in self.population[index]:
                    polarity_value += polarity[residue]
                    hydropathy_value += hydropathy[residue]
                    iso_value += iso[residue]
                behavior['polarity'] = polarity_value
                behavior['hydropathy'] = hydropathy_value
                behavior['iso'] = iso_value
                '''
                # behavior includes euclidean distances to the global peak (based on polarity, hydropathy, and isoelectricity)
                sequence_behavior[seq] = euclidean_distance(seq, self.global_peak)

        for index in range(self.pop_size):
            # Find each sequence's k-nearest neighbors in the population and archive combined.
            seq = arr_to_str(self.population[index])
            compare_pop = np.concatenate((self.archive, np.delete(self.population, index, 0)), axis=0)
            distances_to_compare_pop = []
            for n in range(compare_pop.shape[0]):
                # Distance is evaluated by normalized chem and hamming distances (see top)
                compare_seq = arr_to_str(compare_pop[n])
                distances_to_compare_pop.append(
                    #abs(sequence_behavior[seq]['polarity'] - sequence_behavior[compare_seq]['polarity'])
                    #get_distance(seq, compare_seq)
                    abs(sequence_behavior[seq]['all'] - sequence_behavior[compare_seq]['all'])
                )
            distances_to_compare_pop.sort() #  small -> large distances
            #print(distances_to_compare_pop[0:10])

            # Calculate the sequence's average distance to the k-nearest neighbors (i.e., the sparsity).
            sparsity = 0
            for dist in distances_to_compare_pop[:nearest_neighbors]:
                sparsity += dist
            pop_novelty.append(sparsity/nearest_neighbors)

            # Update the archive.
            if archive_method == 1:
                # Method 1: Fixed size random archive.
                if self.rng.random() < prob:
                    # If the archive is full, remove the oldest individual.
                    if self.archive.shape[0] > self.pop_size * 2:
                        self.archive = np.delete(self.archive, 0, axis=0)
                    self.archive = np.concatenate((self.archive, np.array([self.population[index]])), axis=0)
            else:
                # Method 2: Unrestricted archive size using adaptive threshold.
                num_eval += 1
                if sparsity > threshold:
                    self.archive = np.concatenate((self.archive, np.array([self.population[index]])), axis=0)
                    counter += 1
                if num_eval == 20:
                    if counter > add_max:
                        threshold += fraction
                    elif counter < add_min:
                        threshold -= fraction
                    # Reset counters
                    num_eval = 0
                    counter = 0

        # Return list of novelty values.
        # print(pop_novelty)
        return pop_novelty

    def novelty_selection(self, nearest_neighbors: int = 10, archive_method: int = 1, prob: float = 0.10) -> None:
        """
        Selects the top 10 sequences by novelty value, then replaces the population with the elites (most novel).
        """
        pop_novelty = self.population_novelty(nearest_neighbors=nearest_neighbors, archive_method=archive_method,
                                              prob=prob)
        self.elite = get_elites(self.population, pop_novelty, self.land)
        self.elite1 = max(self.elite, key=lambda x: x[2])
        self.replace_pop()
        return

    def population_combo(self, weight: float = 0.5, nearest_neighbors: int = 10, archive_method: int = 1,
                         prob: float = 0.10) -> list:
        """
        Evaluate the sequences in a population based on both fitness and novelty.
        This combination score is calculated by score(i) = (1 − ρ) · fit(i) + ρ · nov(i), where ρ in [0, 1] controls
        the relative importance of fitness and novelty. High p biases towards novelty.
        Fitness and novelty are normalized according to:
            norm_metric(i) = (metric(i) - metric_min) / (metric_max - metric_min)
        where metric_min and metric_max are the lowest and highest fitness/novelty in the population.

        Returns a list containing only the scores of the sequences in the population.
        Combo scores have the same index in this list as the index of the sequence in the population.

        Returns:
            list
        """
        if weight < 0 or weight > 1:
            raise ValueError('The weight must be a number between 0 and 1.')
        pop_combo = []
        objective = population_fitness(self.population, self.land)
        novelty = self.population_novelty(nearest_neighbors=nearest_neighbors, archive_method=archive_method, prob=prob)

        # Normalize metrics, then calculate the score.
        max_fitness, min_fitness = max(objective), min(objective)
        max_novelty, min_novelty = max(novelty), min(novelty)
        for n in range(self.pop_size):
            norm_fitness = (objective[n] - min_fitness) / (max_fitness - min_fitness)
            norm_novelty = (novelty[n] - min_novelty) / (max_novelty - min_novelty)
            score = (1 - weight) * norm_fitness + weight * norm_novelty
            pop_combo.append(score)
        return pop_combo

    def combo_selection(self, weight: float = 0.5, nearest_neighbors: int = 10, archive_method: int = 1,
                        prob: float = 0.10) -> None:
        """
        Selects the top 10 sequences by combo value, then replaces the population with the elites.
        """
        pop_combo = self.population_combo(weight=weight, nearest_neighbors=nearest_neighbors,
                                          archive_method=archive_method, prob=prob)
        self.elite = get_elites(self.population, pop_combo, self.land)
        self.elite1 = max(self.elite, key=lambda x: x[2])
        self.replace_pop()
        return
