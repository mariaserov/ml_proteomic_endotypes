import numpy as np
import pandas as pd 
import os
import matplotlib.pyplot as plt

# print(os.getcwd())

ks = [3,4,5,6]
x_grid = "X_grid"
df = pd.DataFrame({})
for i, k in enumerate(ks):
    y = pd.read_csv(f"CGGM_array_job_k{k}/CDF_k{k}.csv")
    if i ==1:
        df[x_grid] = y[x_grid]
    else:
        pass
    df[f"{k}"]= y['CDF']

print(df.head())
df.to_csv("CDFs_combined.csv")

# df = pd.read_csv("../CDFs_combined.csv")


for col in ["3", "4", "5", "6"]:
    if col != 'X_grid':  # Exclude the x-axis column
        plt.plot(df["X_grid"], df[col], label=f"k={col}")

plt.xlabel("Consensus index value")
plt.ylabel("CDF")
plt.title("Consensus matrix CDFs for k=3 to k=6", fontsize=14)
plt.legend()

plt.savefig(f"CDF_k2-k6.png", dpi=300, bbox_inches='tight')

plt.show()

print(os.getcwd())