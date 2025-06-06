import numpy as np
import matplotlib.pyplot as plt
import time 
import pandas as pd
from scipy.stats import gaussian_kde

k=3
criterion = 0.9

start = time.time()
np.random.seed(42)

m = np.load("CombinedSorted_ConsMatrix.npy")

cons = m[np.triu_indices_from(m, k=1)]

cons_sorted = np.sort(cons)
kde = gaussian_kde(cons_sorted)
x_grid = np.linspace(cons_sorted.min(), cons_sorted.max(), 1000)
pdf = kde(x_grid)
cdf = np.cumsum(pdf)
cdf = cdf / cdf[-1] 

plt.plot(x_grid, cdf, label="KDE-based CDF")
plt.xlabel("Consensus Index")
plt.ylabel("Cumulative Proportion")
plt.title(f"Consensus CDF ({k} Clusters)")
plt.savefig(f"CDF_k{k}.png", dpi=300, bbox_inches='tight')

df = pd.DataFrame({})
df['X_grid'] = x_grid
df["CDF"] = cdf

df.to_csv(f"CDF_k{k}.csv", index=False)
auc_trapz = np.trapezoid(cdf, x_grid)
print(f"Area under CDF: {auc_trapz}")