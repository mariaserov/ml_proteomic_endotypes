import numpy as np
import matplotlib.pyplot as plt
import time 
import seaborn as sns

k=2
n=5 
start = time.time()
np.random.seed(42)

m = np.load("CombinedSorted_ConsMatrix.npy")
sns.heatmap(m, cmap="Blues", square=True, cbar_kws={'shrink': 0.6})
plt.xticks([])
plt.yticks([])
plt.title(f'Consensus matrix, k={k}', fontsize = 11)
plt.savefig(f"ConsGMM_3clust.png", dpi=300, bbox_inches='tight')

print(f"Runtime: {time.time()-start}s")