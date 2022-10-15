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

# Parameters -------------------------------------------------------------------
parser = argparse.ArgumentParser(description='Landscape complexity scores')
parser.add_argument('-land', '--landscape_file', required=True, help='landscape file. Required')
parser.add_argument('-l', '--hap_len', action = 'store_true', help='Number of variable sites')
parser.add_argument('-n', '--hap_num', action = 'store_true', help='Number of haplotypes')
parser.add_argument('-p', '--peaks', action = 'store_true', help='Maximum number of peaks')
parser.add_argument('-rs', action = 'store_true', help='Roughness')

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
    
sys.exit()