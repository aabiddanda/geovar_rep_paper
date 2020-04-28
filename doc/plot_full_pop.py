# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# +
import numpy as np 
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import os
from scipy.cluster.hierarchy import dendrogram, linkage

plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['figure.facecolor'] = "w"
plt.rcParams['figure.autolayout'] = True

from mpl_toolkits.axes_grid1 import make_axes_locatable

# Deboxing a particular axis
def debox(ax):
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
figdir = '../plots/figure1/'
os.makedirs(figdir, exist_ok=True)

# %matplotlib inline
# -

# # Frequencies in the 26 Population Setting

# +
# TODO : sort by global MAF or DAF and plot heat as the log10 value
data_dir = '../data/freq/'
freq_file = 'new_1kg_nyc_hg38_filt_chr22.biallelic_snps.pops.freq.txt.gz'

# Just arbitrarily get 500000 rows of the file
freq_chr22_df = pd.read_csv(data_dir + freq_file, sep='\t', nrows=500000)
freq_chr22_df.head()
pops = list(freq_chr22_df.columns[6:])

# using the original 1KG populations
# pops = ['YRI', 'LWK', 'GWD', 'MSL', 'ESN', 'ASW', 'ACB', 'CEU', 'TSI', 'FIN', 'GBR',
#        'IBS', 'GIH', 'PJL', 'BEB', 'STU', 'ITU', 'CHB', 'JPT', 'CHS',
#        'CDX', 'KHV', 'MXL', 'PUR', 'CLM', 'PEL']

maf = freq_chr22_df['MAF'].values
af_mat = freq_chr22_df[pops].values

# Adding in a little psuedo-frequency
af_mat = af_mat + 1e-4

# +
# choose N random variants
np.random.seed(420)
n = 100
idx = np.random.choice(maf.shape[0], size=n)
maf_subsets = maf[idx]

# choose the sorted indices
idx_sorted = idx[np.argsort(maf_subsets)[::-1]]

fig, ax = plt.subplots(1,1,figsize=(6,8))

cmap = mpl.colors.ListedColormap(['white', '#94C4DF','#64AAD3', '#1F6EB3', '#0A509A', '#09306B', 'black'])
bounds = np.array([0, 1e-3, 1e-2, 0.05, 0.25, 0.5, 0.75, 1.0])
norm = mpl.colors.BoundaryNorm(bounds, cmap.N, clip=True)

im = ax.imshow(af_mat[idx_sorted], aspect='auto', cmap=cmap, norm=norm)

# Adding in the hlines
for i in range(n):
    ax.axhline(y=i+0.5,color='gray', lw=0.25)

# Putting in the nice vertical lines
ax.axvline(x=-0.5,color='gray', lw=0.25)
for i in range(len(pops)):
    ax.axvline(x=i+0.5,color='gray', lw=0.25)
ax.axvline(x=i+0.5,color='gray', lw=0.25)

# adding in the superpopulation lines
superpop_lns = [4, 8, 13, 18]
for i in superpop_lns:
    ax.axvline(x = i + 0.5, color='gray', lw=1.0)

# Setting x-locations for the superpopulation labels
superpop_lbls = ['AFR', 'EUR', 'SAS', 'EAS', 'AMR']
xpts_superpops = [0.11, 0.23, 0.36, 0.51, 0.69]
for i in range(len(superpop_lbls)):
    fig.text(x=xpts_superpops[i] , y=-0.01, s=superpop_lbls[i], fontsize=16,  va='center',  ha='center')

# setting the ticks + the underlying population labels
ax.set_xticks(np.arange(len(pops)))
ax.set_xticklabels(pops, fontsize=10, rotation=90)
ax.set_yticks([]);
debox(ax);

# Setting the colorbar 
cmap = mpl.colors.ListedColormap(['white', '#94C4DF','#64AAD3', '#1F6EB3', '#0A509A', '#09306B', 'black'])
ticklocs = [0, 0.05, 0.1, 0.175, 0.35, 0.55, 0.75, 1.0]
bounds = np.array(ticklocs)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm._A = []
cbar = fig.colorbar(sm,
                    fraction=0.035, 
                    pad=0.05,
                    shrink=0.6,
                    spacing='proportional',
                    orientation='vertical')

tick_locator = mpl.ticker.FixedLocator(locs=ticklocs)
cbar.locator = tick_locator
cbar.update_ticks()
cbar.ax.set_yticklabels([0, 1e-3, 1e-2, 0.05, 0.25, 0.5, 0.75,  1.0],fontsize=10)
cbar.set_label(label="Globally minor allele frequency", fontsize=16)


plt.tight_layout()
plt.savefig(figdir + 'fig1.pdf', dpi=300, bbox_inches='tight')


