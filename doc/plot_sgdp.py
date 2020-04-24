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
# import matplotlib.patches as patches
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

figdir = '../plots/sgdp_figs/'
os.makedirs(figdir, exist_ok=True)
# -

# # Reading SGDP Project Data

# +
# Developing pairs of SGDP individuals to test through
SGDP_ID_FILE = '../params/parfiles/sgdp_id_pop.txt'
SGDP_DF = pd.read_csv(SGDP_ID_FILE, sep='\s+')
SGDP_IDS = SGDP_DF['Illumina_ID']
POP_IDS = SGDP_DF['Population_ID']

Yoruba_idxs = np.where(POP_IDS.values == 'Yoruba')[0]
Han_idxs = np.where(POP_IDS.values == 'Han')[0]
French_idxs = np.where(POP_IDS.values == 'French')[0]

Yoruba_names = SGDP_IDS.values[Yoruba_idxs]
Han_names = SGDP_IDS.values[Han_idxs]
French_names = SGDP_IDS.values[French_idxs]

sgdp_paired_total = []
sgdp_paired_total_alt = []
sgdp_pop_pairs = []
for i in Yoruba_names[0:1]:
    for j in Yoruba_names[1:]:
        sgdp_paired_total.append('%s_%s' % (i,j))
        sgdp_paired_total_alt.append('%s\n%s' % (i,j))
        sgdp_pop_pairs.append('Yoruba/Yoruba')

for i in Yoruba_names[0:1]:
    for j in Han_names:
        sgdp_paired_total.append('%s_%s' % (i,j))
        sgdp_paired_total_alt.append('%s\n%s' % (i,j))
        sgdp_pop_pairs.append('Yoruba/Han')

for i in Yoruba_names[0:1]:
    for j in French_names:
        sgdp_paired_total.append('%s_%s' % (i,j))
        sgdp_paired_total_alt.append('%s\n%s' % (i,j))
        sgdp_pop_pairs.append('Yoruba/French')

for i in French_names[0:1]:
    for j in Han_names:
        sgdp_paired_total.append('%s_%s' % (i,j))
        sgdp_paired_total_alt.append('%s\n%s' % (i,j))
        sgdp_pop_pairs.append('French/Han')

for i in French_names[0:1]:
    for j in French_names[1:]:
        sgdp_paired_total.append('%s_%s' % (i,j))
        sgdp_paired_total_alt.append('%s\n%s' % (i,j))
        sgdp_pop_pairs.append('French/French')

for i in Han_names[0:1]:
    for j in Han_names[1:]:
        sgdp_paired_total.append('%s_%s' % (i,j))
        sgdp_paired_total_alt.append('%s\n%s' % (i,j))
        sgdp_pop_pairs.append('Han/Han')
        
print(sgdp_paired_total)
print(sgdp_pop_pairs)
print(len(sgdp_pop_pairs))

# +
# Defining meta-variables
sgdp_file_lists = sgdp_paired_total
pops_str = 'superpops_amended2'
ncat_str = '3x'
prefix = 'new_1kg_nyc_hg38_filt'
cmap_name = 'Blues'
n_sgdp = len(sgdp_file_lists)

# Accumulator for all the PGP GeoDist Objects
sgdp_geodists = []

for i in range(n_sgdp):
    subset = sgdp_file_lists[i]
    input_file = '../data/geodist/counts/subsets/sgdp_paired/%s/%s_%s_%s.geodist_cnt.txt.gz' % (subset, prefix, pops_str, ncat_str)

    cur_geodist = GeoDistPlot()
    cur_geodist._add_text_data(input_file, filt_unobserved=True)
    cur_geodist._add_poplabels_manual(np.array(['AFR','EUR','SAS','EAS','AMR']))
    cur_geodist._filter_data()
    cur_geodist._add_cmap(cmap_name, str_labels=['u','R','C'])
    cur_geodist.fontsize = 8
    sgdp_geodists.append(cur_geodist)


# +
# Choosing indices for relevant population comparisons
idxs = [5,3,0,9,11,13]
subset_geodist = [sgdp_geodists[x] for x in idxs]
subset_pop_pairs = [sgdp_pop_pairs[x] for x in idxs]
plot_multiple_geodist(subset_geodist, subset_pop_pairs, ylabel='Cumulative fraction of 1000G variants', hwidth=0.3);

# Testing the plots out here
plt.savefig(figdir + 'sgdp_panel_geodist.filt_lvl1.pdf', bbox_inches='tight', dpi=300)

# +
# print multiple panels
plot_multiple_geodist(sgdp_geodists[0:5], 
                      sgdp_paired_total_alt[0:5], 
                      ylabel='Cumulative fraction of 1000G variants', hwidth=0.5);
# plt.savefig(figdir + 'paired_sgdp_row1.filt_lvl1.pdf', bbox_inches='tight', dpi=300)


plot_multiple_geodist(sgdp_geodists[5:10], 
                      sgdp_paired_total_alt[5:10], 
                      ylabel='Cumulative fraction of 1000G variants', hwidth=0.5);
# plt.savefig(figdir + 'paired_sgdp_row2.filt_lvl1.pdf', bbox_inches='tight', dpi=300)

plot_multiple_geodist(sgdp_geodists[10:], 
                      sgdp_paired_total_alt[10:], 
                      ylabel='Cumulative fraction of 1000G variants', hwidth=0.5);
# plt.savefig(figdir + 'paired_sgdp_row3.filt_lvl1.pdf', bbox_inches='tight', dpi=300)

# +
# sgdp_file_lists = sgdp_paired_total
# pops_str = 'pops'
# popfile = '../../../params/poplists/%s_panel.txt' % pops_str
# ncat_str = '3x'
# prefix = 'new_1kg_nyc_hg38_filt'
# cmap_name = 'Blues'
# n_sgdp = len(sgdp_file_lists)

# # Accumulator for all the PGP GeoDist Objects
# sgdp_geodists_pops = []

# for i in range(n_sgdp):
#     subset = sgdp_file_lists[i]
#     input_file = '../../../data/geodist_counts_subsets/sgdp_paired/%s/%s_%s_%s.geodist_cnt.txt.gz' % (subset, prefix, pops_str, ncat_str)

#     cur_geodist = GeoDistPlot()
#     cur_geodist._add_text_data(input_file, filt_unobserved=True)
#     cur_geodist._add_poplabels(popfile)
#     cur_geodist._filter_data()
#     cur_geodist._add_cmap(cmap_name, str_labels=['u','R','C'])
#     cur_geodist.fontsize = 12
#     sgdp_geodists_pops.append(cur_geodist)

# +
# idxs = [5,3,0,9,11,13]
# subset_geodist_pops = [sgdp_geodists_pops[x] for x in idxs]
# subset_pop_pairs_pops = [sgdp_pop_pairs[x] for x in idxs]
# plot_multiple_geodist(subset_geodist_pops, 
#                       subset_pop_pairs_pops, 
#                       superpops=[5, 9, 14, 19],
#                       xsize=3, 
#                       ysize=4.5, 
#                       ylabel='Cumulative fraction of 1000G variants',
#                       hwidth=0.2);
# # plt.savefig(figdir + 'sgdp_paired_full26pops.filt_lvl1.pdf', bbox_inches='tight', dpi=300)
# # plt.savefig(figdir + 'sgdp_paired_full26pops.filt_lvl1.png', bbox_inches='tight', dpi=300)
# -

# Calculating the proportion of "globally widespread alleles across these filtered level 1 sets" 
def compute_any_present(cats, ns):
    """ Computes the number of globally present variants  """
    idx = [i for i in range(cats.size) if ('0' not in cats[i])]
    return(cats[idx], ns[idx])


idxs = [5,3,0,9,11,13]
for i in idxs:
    _, n = compute_any_present(sgdp_geodists[i].orig_geodist, sgdp_geodists[i].orig_ngeodist)
    ns = np.sum(n)
    tot_n = np.sum(sgdp_geodists[i].orig_ngeodist) + sgdp_geodists[i].orig_missing
    print(sgdp_pop_pairs[i], ns, tot_n, ns/tot_n)
