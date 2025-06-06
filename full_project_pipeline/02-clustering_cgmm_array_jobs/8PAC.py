import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

criterion = 0.9
lower_bound = 1 - criterion  # Lower threshold
upper_bound = criterion  # Upper threshold
ks = [2,3,4,5,6]

results = pd.DataFrame({
    "k": ks,
    "PAC": [0.0, 0.0, 0.0, 0.0, 0.0]
})

for i, k in enumerate(ks):
    df = pd.read_csv(f"CGGM_array_job_k{k}/CDF_k{k}.csv")
    data = df.to_numpy()


    x_grid = data[:,0]
    cdf = data[:,1]

    # Find closest indices in x_grid
    lower_idx = np.searchsorted(x_grid, lower_bound)
    upper_idx = np.searchsorted(x_grid, upper_bound)
    

    # Get the corresponding CDF values
    PAC_score = cdf[upper_idx] - cdf[lower_idx]

    results.iloc[i,1] = PAC_score

    print(f"PAC Score: {PAC_score}")

bar_width = 0.3
plt.figure(figsize=(8, 5))
plt.bar(ks, results["PAC"], width=bar_width, color='darkblue', alpha=0.6)
plt.xlabel("Number of clusters (k)")
plt.xticks(ks)
plt.ylabel("PAC score")
plt.title("PAC Scores across numbers of clusters (k)", fontsize=14)

plt.savefig(f"PAC_plot.png", dpi=300, bbox_inches='tight')
plt.show()

print(results)
results.to_csv("PAC_scores.csv")