import pandas as pd
import numpy as np
import umap 
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.mixture import GaussianMixture
import time 
import seaborn as sns
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from sklearn.metrics import pairwise_distances
from joblib import Parallel, delayed
import os
import psutil
import sys

k=2
subsamples = 100


n_runs = int(subsamples/10)
n = int(sys.argv[1])
start = time.time()
process = psutil.Process(os.getpid())
np.random.seed(42)
n_jobs = int(os.environ.get("PBS_NCPUS", 8))

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

# Uploading & scaling the data

path = "../../Data/main_hypertension_dataset.csv"

df = pd.read_csv(path)

prot = df.loc[:,'AARSD1':] # For full dataset

# path = "../../../Data/toy_data_all.csv"
# df = pd.read_csv(path)
# prot = df.loc[:,'SRC':] # For full_toy dataset

prot_scaled = StandardScaler().fit_transform(prot)
prot_scaled = pd.DataFrame(prot_scaled, columns=prot.columns)

# Clustering on subsamples

n_samples = prot_scaled.shape[0]
subsample_prop = 0.5
consensus_matrices = {}

# Parallelisaiton:

def run_gmm_iteration(k, run_idx):
    subsample_idx = np.random.choice(n_samples, size=int(subsample_prop * n_samples), replace=False)
    prot_subsample = prot_scaled.iloc[subsample_idx, :]
    gmm = GaussianMixture(n_components=k, n_init=5)
    clusters = gmm.fit_predict(prot_subsample)
    
    subsample_matrix = np.zeros((n_samples, n_samples)) # n of times observations are selected together
    cluster_matrix = np.zeros((n_samples, n_samples)) # n of times observations are clustered together

    for idx1, sample1 in enumerate(subsample_idx):
        for idx2, sample2 in enumerate(subsample_idx):
            subsample_matrix[sample1,sample2] += 1
            if clusters[idx1] == clusters[idx2]:
                cluster_matrix[sample1, sample2] += 1
    if run_idx==0:
        bic = gmm.bic(prot_scaled)
        aic = gmm.aic(prot_scaled)
        print(f"BIC: {bic}")
        print(f"AIC: {aic}")
    else:
        pass

    print(f"Iteration {run_idx} complete")
    return cluster_matrix, subsample_matrix

#  For one k=3:

results = Parallel(n_jobs=8, backend="loky")(delayed(run_gmm_iteration)(k, i) for i in range(n_runs))

total_cluster_matrix = np.sum([res[0] for res in results], axis=0)
total_subsample_matrix = np.sum([res[1] for res in results], axis=0)

consensus_matrix = np.divide(total_cluster_matrix, total_subsample_matrix)
consensus_matrix[np.isnan(consensus_matrix)] = 0

np.save(f"k{k}_ConsMatrix_{n}.npy", consensus_matrix)
print(f"Runtime: {time.time()-start}s")
