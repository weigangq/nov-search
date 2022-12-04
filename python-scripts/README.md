# Scripts
This folder contains python scripts for generating fitness landscapes and simulating evolution on fitness landscapes.

## Required Python packages
All Python packages required for the scripts in this folder are listed in the `requirements.txt` file and are listed individually in the following sections. To install all packages, use one of the commands based on your Python package manager:
- pip: `pip install -r requirements.txt`.
- conda: `conda install --file requirements.txt -c conda-forge -c bioconda`

## Landscapes
nk_landscape.py produces a fitness landscape of binary sequences, along with the normalized fitness and other information about the landscape's peaks.
- Information about NK landscapes can be found at https://doi.org/10.1016/s0022-5193(89)80019-0 .

land-stats.py gives information about a landscape.

### Required Python packages for land-stats.py
- matplotlib
- networkx
- numpy
- pandas
- scipy
- scikit-learn

### Running nk_landscape.py
Run `python nk_landscape.py [args]` to produce the landscape.
The script will print the landscape in the console, with columns separated by tabs.
- Save the landscape using `python nk_landscape.py [args] > nk_landscape.tsv`

#### Args
- `-N` - Number of variable sites. The fitness landscape will contain every binary sequence that can be made with N bits.
- `-K` - Number of interacting sites, from 0 to N-1, inclusive.
- `-t --tag` - Adds a tag to the output. The tag is a string that will appear in a column in the output for every haplotype.
- `-e --exp_fit` - Converts normalized fitness score to a z-score (# of standard deviations from the mean). Without this, the normalized score will be in the range [0, 1] with the global peak having a normalized fitness value of 1.

e.g.: `python nk_landscape.py -N 10 -K 2`

### Running land-stats.py
Run `python land-stats.py [args]` to get the landscape information.

#### Args
- `-land --landscape_file` - (Required) Landscape file in tsv format. The name of the column containing the sequences must be named "haplotype".
- `-l` - Prints the number of variable sites in the haplotypes.
- `-n` - Prints the total number of haplotypes in the landscape.
- `-p` - Prints the number of peaks in the landscape.
- `-rs` - Prints the roughness to slope (r/s) ratio.
- `b` - Prints the basin of attraction.
- `e` - Prints the escape edges.

e.g.: `python land-stats.py -land ../data/nk-landscape.tsv -l -b`

## Simulations
aa_search.py and binary_search.py simulate the evolution of a population of peptides or binary sequences based on a fitness landscape.
- aa_sim_population.py, aa_sim_population_cython.pyx, and binary_sim_population.py contain the code for the functions used in the simulation.

aa_search can optionally be run using the Cython package, which will make novelty search run significantly faster. This requires a C compiler to be installed separately (GCC, MSVC, etc.)

### Required Python packages for aa_search.py and binary_search.py
- numpy
- pandas
- Cython (optional)

### Fitness landscape
A file containing the fitness landscape is required. 
- The landscape file should be a file that can be read into the Python pandas package.
- One column should contain the sequences, and another column should contain the corresponding fitness values.
	- The sequences should all be of the same length.
- An example fitness landscape file is found in ../data/gb1_fitness.tsv, taken from https://doi.org/10.7554/eLife.16965 .

### Running the simulation
To run the simulation, run the command `python aa_search.py [args]`
- `--land` (Required) Landscape file in tsv format.
- `-c` (Optional) Used to make the simulation run faster using Cython.
- The other args can be changed to modify the parameters of the simulation.
	- `-alg [1, 2, or 3]` Pick a number to select the search algorithm.
		- `1` - Objective search
		- `2` - Novelty search
		- `3` - Combo search
	- `-pop --pop_size` - Number of sequences in the population.
	- `-gen --generation` - Number of generations to mutate the population. 
	- `-mut --mutation_rate` - Average number of mutations per generation.
	- `-s --rng_seed` - Set an RNG seed to get reproducible RNG results. If not provided, the results will be unpredictable and not reproducible.
	- `t --tag` Set a tag for the simulation. The tag is a string that will appear in a column in the output for every generation.
	
- Novelty/combo specific parameters
	- `-beh --behavior` - Choose a behavior measure of haplotypes to use for novelty search.
  		- `pol` Polarity
    	- `hydro` Hydropathy index
        - `iso` Isoelectric index
        - `all` Euclidean distance between polarity, hydropathy, and isoelectric.
        - `pol_norm, hydro_norm, iso_norm, all_norm`Normalized (0 to 1) versions of polarity, hydropathy, isoelectric, and all.
        - `blosum` BLOSUM62 matrix
	- `-n --nearest_neighbors` - Number of nearest neighbors to use when calculating novelty search (Novelty = average distance from k-nearest_neighbors).
	- `-arch [1 or 2]` - Type of archive to use during novelty search. 
      - `1` - Fixed size, random archive. 
      - `2` - Adaptive threshold.
	- `-prob --prob_arch` - Probability for an individual to get randomly added to the archive during novelty search method 1.
	- `-w --weight` - Weight for the combo search algorithm. The weight value is between 0 and 1. Closer to 1 means more bias towards novelty.

Example run: `python aa_search.py -land ../data/gb1_fitness.tsv -c -alg 3 -pop 50 -gen 50 -mut 1.5 -w 0.3 -beh all_norm`

The script will print information about the fitness landscape's global peak, followed by the results of the simulation, with columns separated by tabs.
- The results information contains each generation's highest fitness sequence, along with the fitness value.

