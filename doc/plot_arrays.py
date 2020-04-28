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
import matplotlib.patches as patches
import matplotlib as mpl
import pandas as pd

import os
import sys
sys.path.append('../src/')
from plotting import GeoDistPlot, plot_multiple_geodist

# %matplotlib inline

# +
# Plotting preferences
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['figure.facecolor'] = "w"
plt.rcParams['figure.autolayout'] = True

from mpl_toolkits.axes_grid1 import make_axes_locatable

# Deboxing a particular axis
def debox(ax):
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

figdir = '../plots/array_figs/'
os.makedirs(figdir, exist_ok=True)
# -

# # Plotting the Array Figure (MAF 5% Boundary)

# +
# %%time
# Relevant Parameters
array_file_lists = ["Affy6", 
                    "HumanOrigins", 
                    "OmniExpress", 
                    "Omni25Exome", 
                    "MEGA"]

pops_str = 'superpops_amended2'
# popfile = '../../../params/poplists/%s_panel.txt' % pops_str
ncat_str = '3x'
prefix = 'new_1kg_nyc_hg38_filt'
cmap_name = 'Blues'
n_arrays = len(array_file_lists)

# Storing GeoDist Objects
array_geodists = []

for i in range(n_arrays):
    subset = array_file_lists[i]
    input_file = '../data/geodist/counts/subsets/%s/%s_%s_%s.geodist_cnt.txt.gz' % (subset, prefix, pops_str, ncat_str)
    cur_geodist = GeoDistPlot()
    cur_geodist._add_text_data(input_file, filt_unobserved=True)
    cur_geodist._add_poplabels_manual(np.array(['AFR','EUR','SAS','EAS','AMR']))
    cur_geodist._filter_data(max_freq=0.005)
    cur_geodist._add_cmap(cmap_name, str_labels=['u','R','C'])
    cur_geodist.fontsize = 8
    array_geodists.append(cur_geodist)
# -

f, axs = plot_multiple_geodist(array_geodists, array_file_lists, ylabel = 'Cumulative fraction of 1000G variants', hwidth=0.3) 
plt.savefig(figdir + 'array_fig.5percent.pdf', bbox_inches='tight')
