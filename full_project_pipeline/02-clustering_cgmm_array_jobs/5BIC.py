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

ks = [2,3,4,5,6]

path = "../Data/full_dataset.csv"
df = pd.read_csv(path)
df = df.loc[df["Systolic.blood.pressure..automated.reading.0.0"]>=float(140)] # Filter for hypertension
prot = df.loc[:,'AARSD1':] # For full dataset
prot_scaled = StandardScaler().fit_transform(prot)
prot_scaled = pd.DataFrame(prot_scaled, columns=prot.columns)

results = pd.DataFrame({
    "k": np.zeros(4, dtype=float),
    "BIC": np.zeros(4, dtype=float), 
    "AIC": np.zeros(4, dtype=float)
})

for i, k in enumerate(ks):
    gmm = GaussianMixture(n_components=k, n_init=5)
    gmm.fit(prot_scaled)
    bic = gmm.bic(prot_scaled)
    
    aic = gmm.aic(prot_scaled)

    results.iloc[i,0] = k
    results.iloc[i,1] = bic
    results.iloc[i,2] = aic

print(results)
results.to_csv("GMMs_AIC_BIC.csv", index=False)
