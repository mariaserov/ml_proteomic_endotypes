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
n=5 # number of consensus submatrices

start = time.time()
np.random.seed(42)

path = "../../Data/main_hypertension_dataset.csv"
df = pd.read_csv(path)

matrices = []

for i in range(1,n+1):
    print(f"Subset of matrix {i}:")
    m = np.load(f"k{k}_ConsMatrix_{i}.npy")
    print(m[0:5,0:5])
    matrices.append(m)

matrix = matrices[0]

for i in range(1,n):
    matrix = matrix + matrices[i]

print("Added matrix:")
print(matrix[0:5,0:5])

matrix = matrix/n

print("Finalised matrix:")
print(matrix[0:5,0:5])

distance_matrix = 1 - matrix
linkage_matrix = linkage(distance_matrix, method='average')
dendro = dendrogram(linkage_matrix, no_plot=True)
ordered_indices = dendro['leaves']
sorted_consensus_matrix = matrix[ordered_indices, :][:, ordered_indices]

final_clusters = fcluster(linkage_matrix, t=2, criterion='maxclust')
clust = pd.DataFrame({
    "index": df.index,
    "Cluster": final_clusters
})

clust.to_csv("cluster_assignments_k3.csv", index=False)
np.save("CombinedSorted_ConsMatrix.npy", sorted_consensus_matrix)