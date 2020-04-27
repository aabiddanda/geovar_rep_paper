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
import matplotlib.pyplot as plt
import matplotlib as mpl
import gzip as gz

# Importing the geodist plot 
import os
import sys
sys.path.append('../src/')
from plotting import GeoDistPlot, plot_multiple_geodist, plot_geodist_w_percentages

# %matplotlib inline

# +
# Setting plotting parameters
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['figure.facecolor'] = "w"
plt.rcParams['figure.autolayout'] = True
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Deboxing a particular axis
def debox(ax):
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
# Setting the x-size and y-size for the entire figure...
figdir = '../plots/overall_figs/'
os.makedirs(figdir, exist_ok=True)
# -

# # Plotting GeoDist (Superpops)
# ## Using MAF = 5% as the cutoff
# ### Using All Variants

# +
# %%time

#1. Loading superpopulation data
superpops_geodist_file = '../data/geodist/counts/total/new_1kg_nyc_hg38_filt_total.biallelic_snps.superpops_amended2.ncat3x.filt_0.geodist_cnt.txt.gz'
geodist_superpops = GeoDistPlot()
geodist_superpops._add_text_data(superpops_geodist_file, filt_unobserved=False)
geodist_superpops._add_poplabels_manual(np.array(['AFR','EUR','SAS','EAS','AMR']))

# 3. Filtering Data (largely for plotting speed) 
geodist_superpops._filter_data(max_freq=0.0001)

# Set colormaps for geodist
cmap_name = 'Blues'
geodist_superpops._add_cmap(cmap_name, str_labels=['u', 'R','C'])
title = 'hg38_filtered'
pops_str = 'superpops'
ncat_str = '3x'

fig,axs = plot_geodist_w_percentages(geodist_superpops, subset='Biallelic SNVs', 
                                     ylabel='Cumulative fraction of 1000G variants', 
                                     xsize=1.5, ysize=4)
plt.savefig(figdir + 'superpops.5percent.all.pdf', bbox_inches='tight')
# -

# ### Removing Singletons from the Data

# +
#1. Loading superpopulation data
superpops_geodist_file = '../data/geodist/counts/total/new_1kg_nyc_hg38_filt_total.biallelic_snps.superpops_amended2.ncat3x.filt_1.geodist_cnt.txt.gz'
geodist_superpops = GeoDistPlot()
geodist_superpops._add_text_data(superpops_geodist_file, filt_unobserved=False)
geodist_superpops._add_poplabels_manual(np.array(['AFR','EUR','SAS','EAS','AMR']))

# 3. Filtering Data (largely for plotting speed) 
geodist_superpops._filter_data(max_freq=0.0001)

# Set colormaps for geodist
cmap_name = 'Blues'
geodist_superpops._add_cmap(cmap_name, str_labels=['u', 'R','C'])
title = 'hg38_filtered'
pops_str = 'superpops'
ncat_str = '3x'

fig,axs = plot_geodist_w_percentages(geodist_superpops, subset='Singletons Removed', 
                                     ylabel='Cumulative fraction of 1000G variants', 
                                     xsize=1.5, ysize=4.5)
plt.savefig(figdir + 'superpops.5percent.nosingletons.pdf', bbox_inches='tight')
# -

# ## Plots with an MAF > 1% boundary

# +
#1. Loading superpopulation data
superpops_geodist_file = '../data/geodist/counts/total/new_1kg_nyc_hg38_filt_total.biallelic_snps.superpops_amended2.ncat3.filt_0.geodist_cnt.txt.gz'
geodist_superpops = GeoDistPlot()
geodist_superpops._add_text_data(superpops_geodist_file, filt_unobserved=False)
geodist_superpops._add_poplabels_manual(np.array(['AFR','EUR','SAS','EAS','AMR']))

# 3. Filtering Data (largely for plotting speed) 
geodist_superpops._filter_data(max_freq=0.0001)

# Set colormaps for geodist
cmap_name = 'Blues'
geodist_superpops._add_cmap(cmap_name, str_labels=['u', 'R','C'])
title = 'hg38_filtered'
pops_str = 'superpops'
ncat_str = '3x'

fig,axs = plot_geodist_w_percentages(geodist_superpops, subset='Biallelic SNVs', 
                                     ylabel='Cumulative fraction of 1000G variants', 
                                     xsize=1.5, ysize=4)
plt.savefig(figdir + 'superpops.1percent.all.pdf', bbox_inches='tight')

# +
#1. Loading superpopulation data
superpops_geodist_file = '../data/geodist/counts/total/new_1kg_nyc_hg38_filt_total.biallelic_snps.superpops_amended2.ncat3.filt_1.geodist_cnt.txt.gz'
geodist_superpops = GeoDistPlot()
geodist_superpops._add_text_data(superpops_geodist_file, filt_unobserved=False)
geodist_superpops._add_poplabels_manual(np.array(['AFR','EUR','SAS','EAS','AMR']))

# 3. Filtering Data (largely for plotting speed) 
geodist_superpops._filter_data(max_freq=0.0001)

# Set colormaps for geodist
cmap_name = 'Blues'
geodist_superpops._add_cmap(cmap_name, str_labels=['u', 'R','C'])
title = 'hg38_filtered'
pops_str = 'superpops'
ncat_str = '3x'

fig,axs = plot_geodist_w_percentages(geodist_superpops, subset='Singletons Removed', 
                                     ylabel='Cumulative fraction of 1000G variants', 
                                     xsize=1.5, ysize=4)
plt.savefig(figdir + 'superpops.1percent.nosingletons.pdf', bbox_inches='tight')
# -

# # Full Population Plots

# +
# %%time
#1. Loading superpopulation data
pops_geodist_file = '../data/geodist/counts/total/new_1kg_nyc_hg38_filt_total.biallelic_snps.pops.ncat3x.filt_0.geodist_cnt.txt.gz'
pops_panel = '../params/parfiles/poplists/pops_panel.txt'

geodist_pops = GeoDistPlot()
geodist_pops._add_text_data(pops_geodist_file, filt_unobserved=False)
geodist_pops._add_poplabels(pops_panel)
geodist_pops._filter_data(max_freq = 0.001)

# Set colormaps for geodist
cmap_name = 'Blues'
geodist_pops._add_cmap(cmap_name, str_labels=['u', 'R','C'])
geodist_pops.fontsize = 10

# +
fig,axs = plot_geodist_w_percentages(geodist_pops, subset='Biallelic SNVs', 
                                     ylabel='Cumulative fraction of 1000G variants', 
                                     xsize=4, ysize=7, superpops=[5, 9, 14, 19])

plt.savefig(figdir + 'pops.5percent.all.pdf', bbox_inches='tight')
