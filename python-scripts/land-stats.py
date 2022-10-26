#!/usr/bin/env python

"""NK model and RMF model
Referenced paper:
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9018209/

Original file is located at
    https://github.com/song88180/fitness-landscape-error
    https://github.com/Mac13kW/NK_model/blob/master/NK%20landscapes%20-%20a%20hands-on%20exercise%202019.pdf
"""
import numpy as np
import numpy.random as nrand
import pandas as pd
import multiprocessing
import argparse
from sklearn.linear_model import LinearRegression
import sys
import logging
import copy
import networkx as nx

# Parameters -------------------------------------------------------------------
parser = argparse.ArgumentParser(description='Landscape complexity scores')
parser.add_argument('-land', '--landscape_file', required=True, help='landscape file. Required')
parser.add_argument('-l', '--hap_len', action = 'store_true', help='Number of variable sites')
parser.add_argument('-n', '--hap_num', action = 'store_true', help='Number of haplotypes')
parser.add_argument('-p', '--peaks', action = 'store_true', help='Maximum number of peaks')
parser.add_argument('-rs', action = 'store_true', help='Roughness')
parser.add_argument('-b', action = 'store_true', help='Basin of attraction')

args = parser.parse_args()
logging.basicConfig(level=logging.DEBUG)

df = pd.read_csv(args.landscape_file, sep="\t", dtype={'haplotype': str})
seq_len = len(df.at[0, 'haplotype'])
recs = df.to_dict('records')

if args.hap_len is True:
    print(seq_len)

if args.hap_num is True:
    print(df.shape)

haps = []
for rec in recs:
    ar = list(rec['haplotype'])
    ar.append(rec['fitness'])
    haps.append(ar)

landscape = np.array(haps, dtype=np.float64)
#print(landscape.shape)

# Ruggedness -------------------------------------------------------
# Number of maxima (N max)
def get_N_max(landscape):
    N = landscape.shape[1] - 1
    N_max = 0
    for gt in landscape:
        seq = gt[0:N]
        fit = gt[N]
        flag = True
        for i,_ in enumerate(seq):
            seq_ = copy.deepcopy(seq)
            seq_[i] = 1 - seq_[i]
            tmp = ''.join(seq_.astype(int).astype(str))
            idx = int(tmp, 2)
            fit_ = landscape[idx,N]
            if fit < fit_:
                flag = False
                break
        if flag == True:
            N_max += 1
    return N_max

# Roughness to slope ratio (r/s)
def cal_r_s(landscape):
    N = landscape.shape[1] - 1
    X = landscape[:,:N]
    y = landscape[:,-1]
    reg = LinearRegression().fit(X, y) # fit_intercept default=True
    y_predict = reg.predict(landscape[:,:N])
    roughness = np.sqrt(np.mean(np.square(y - y_predict)))
    slope = np.mean(np.abs(reg.coef_))
    return round(roughness/slope,6)

if args.peaks is True:
    print(get_N_max(landscape))

if args.rs is True:
    print(cal_r_s(landscape))

# Calculate basin of attraction
def fit_diff(G, nd, suc):
    return G.nodes[suc]['fit'] - G.nodes[nd]['fit'] # must be positive

# recursive function
def greedy_climb(G, nd, pt):
    if nd in peaks: # a peak, no successors; but inclusive as a basin hap
        pt.append(nd)
        return
    best = sorted( [suc for suc in G.successors(nd)], key = lambda x: fit_diff(G, nd, x), reverse = True) # for each m -> n
    pt.append(best[0])
    greedy_climb(G, best[0], pt)

if args.b is True:
    # Construct graph
    recs = df.to_dict("records")
    dHap = {}
    for rec in recs:
        dHap[rec['haplotype']] = {'id': rec['id'], 'fit': rec['fitness']}

    DG = nx.DiGraph() 
    for hap in dHap: # hap as key
        DG.add_node(hap)
        for k in dHap[hap]: # k includes 'id' and 'fit'
            DG.nodes[hap][k] = dHap[hap][k]
        DG.nodes[hap]['alpha'] = dHap[hap]['fit'] # assign fitness as alpha (transparency)
    
    def add_fit_edge(G, hap, seen, level):
        num_zeros = len([x for x in hap if x == '0'])
        if num_zeros == 0: # reached the end hap (all 1's)
            seen[hap] = 1
            return
        if hap in seen_nodes: # node already reached
            return
        else: # contains at least one zero & not seen
            seen[hap] = 1
            G.nodes[hap]['subset'] = level # networkx uses node feature 'subset' to store levels for multipartite layout
            nabes = [] # start a list of 1-hamming neighbors
            for i in range(len(hap)): # mutate each 0 to 1
                if hap[i] == '0':
                    mut = hap[0:i] + '1' + hap[(i+1):]
                    if G.nodes[hap]['fit'] > G.nodes[mut]['fit']:
                        G.add_edge(mut, hap)
                    else:
                        G.add_edge(hap, mut)
                    nabes.append(mut)
            level += 1
            for x in nabes: # recurse on each child node
                if x not in seen: # skip if seen
                    G.nodes[x]['subset'] = level
                    add_fit_edge(G, x, seen, level)
    seen_nodes = {}
    add_fit_edge(DG, '0000000000', seen_nodes, 0) # start point

    # node coloring by fitness 
    node_alphas = [ DG.nodes[x]['fit'] for x in DG.nodes]
    peaks = [x[0] for x in DG.out_degree if x[1] == 0]
    vals = [x[0] for x in DG.in_degree if x[1] == 0]

    # Basin of attraction (should partition the haplotype space)
    basins = {} # nuumber of paths to a peak (through greedy hill climbing)
    for p in peaks: # initalize dict
        basins[p] = []

    for nd_i in DG.nodes: # 1024 starting nodes
        #start = np.random.default_rng().choice(DG.nodes)
        path = []
        greedy_climb(DG, nd_i, path) # recursive until a peak
        peak = path[-1] # last item is the peak
        basins[peak].append(path)
        #print(path)

    print("id", "\t", "haplotype", "\t", "basin_size", "\t", "fitness")
    for p in basins:
        print(dHap[p]['id'], "\t", p, "\t", len(basins[p]), "\t", round(dHap[p]['fit'],6))

sys.exit()