import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt

bic = pd.read_csv("GMMs_AIC_BIC.csv")
pac = pd.read_csv("PAC_scores.csv")


# x = np.array(bic["k"])
# bic = np.array(bic["BIC"])
# pac = np.array(pac["PAC"])
# aic = np.array(bic["AIC"])

# bar_width = 0.2
# fig, ax1 = plt.subplots(figsize=(8, 5))

# ax1.bar(x - bar_width/2, bic, width=bar_width, color='darkblue', alpha=0.6, label="BIC")
# ax1.set_xlabel("Number of Clusters (k)")
# ax1.set_xticks(x)
# ax1.set_ylabel("BIC", color='darkblue')
# ax1.tick_params(axis='y', labelcolor='darkblue')

# ax2 = ax1.twinx()  
# ax2.bar(x + bar_width/2, pac[1:], width=bar_width, color='red', alpha=0.6, label="PAC")
# ax2.set_ylabel("PAC", color='red')
# ax2.tick_params(axis='y', labelcolor='red')

# plt.title("BIC and PAC vs Number of Clusters")

# plt.savefig(f"BIC_PAC_plot.png", dpi=300, bbox_inches='tight')
# plt.show()

x = np.array(bic["k"])
bic = np.array(bic["BIC"])
bic = bic/bic.max()
pac = np.array(pac["PAC"])
pac = pac/pac.max()
aic = np.array(bic["AIC"])
aic = aic/aic.max()

bar_width = 0.2
fig, ax1 = plt.subplots()

plt.figure(figsize=(8, 5))



plt.bar(x - bar_width/3, bic, width=bar_width, color='darkblue', alpha=0.6, label="BIC")
plt.bar(x + bar_width/3, aic, width=bar_width, color='yellow', alpha=0.6, label="AIC")
plt.bar(x + bar_width, pac[1:], width=bar_width, color='red', alpha=0.6, label="PAC")

plt.set_xlabel("Number of Clusters (k)")
plt.set_xticks(x)
plt.set_ylabel("Value relative to max", color='darkblue')
plt.legend()

plt.title("BIC, AIC, PAC vs Number of Clusters")

plt.savefig(f"BIC_PAC_plot.png", dpi=300, bbox_inches='tight')
plt.show()