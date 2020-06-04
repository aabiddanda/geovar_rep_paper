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

import os
import sys
sys.path.append('../src/')
from theory import geodist_counts
from plotting import GeoDistPlot

# %matplotlib inline

# +
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['figure.facecolor'] = "w"
plt.rcParams['figure.autolayout'] = True

from mpl_toolkits.axes_grid1 import make_axes_locatable

# Deboxing a particular axis
def debox(ax):
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

figdir = '../plots/figure4/'
os.makedirs(figdir, exist_ok=True)
# -

# Setting up simulation parameters for divergence times ...
n = 200; j=8
t1=0.05
t2=0.5
theta=1e-3

# %%time 
# Compute Joint SFS with a divergence time of 0.1
jsfs1 = geodist_counts(t=t1, n=n, k_rare=j, c=0.0, theta=theta)
cur_geodist1 = GeoDistPlot()
cur_geodist1._add_data_jsfs(jsfs=jsfs1)
cur_geodist1._add_poplabels_manual(poplabels=np.array(['A','B']))
cur_geodist1.ncat = cur_geodist1.ncat -1
cur_geodist1._add_cmap(str_labels=['u','R','C'])
cur_geodist1._filter_data(max_freq=0.0)

# %%time
# Compute Joint SFS with a divergence time of 0.5
jsfs2 = geodist_counts(t=t2, n=n, k_rare=j, c=0.0, theta=theta)
cur_geodist2 = GeoDistPlot()
cur_geodist2._add_data_jsfs(jsfs=jsfs2)
cur_geodist2._add_poplabels_manual(poplabels=np.array(['A','B']))
cur_geodist2.ncat = cur_geodist2.ncat -1
cur_geodist2._add_cmap(str_labels=['u','R','C'])
cur_geodist2._filter_data(max_freq=0.0)

# %%time 
# Compute Joint SFS with a divergence time of 0.1
jsfs3 = geodist_counts(t=t1, n=n, k_rare=j, c=0.025, theta=theta)
cur_geodist3 = GeoDistPlot()
cur_geodist3._add_data_jsfs(jsfs=jsfs3)
cur_geodist3._add_poplabels_manual(poplabels=np.array(['A','B']))
cur_geodist3.ncat = cur_geodist3.ncat -1
cur_geodist3._add_cmap(str_labels=['u','R','C'])
cur_geodist3._filter_data(max_freq=0.0)


# ## Developing Plotting Functions for Box 2

# Plotting Divergence Model
def plot_div_model(ax, tdiv=0.3, admix=0.0, vbar=False, vbarlabel=None,  **kwargs):
    # Plotting the underlying functions
    ax.plot([0.25,0.35], [0.0,tdiv], **kwargs)
    ax.plot([0.3,0.4], [0.0,tdiv], **kwargs)
    ax.plot([0.5,0.4], [0.0,tdiv],  **kwargs)
    ax.plot([0.55,0.45], [0.0,tdiv], **kwargs)
    ax.plot([0.35,0.35], [tdiv,0.6], **kwargs)
    ax.plot([0.45,0.45], [tdiv,0.6], **kwargs)
    
    # Plotting the limits
    if vbar:
        x = 0.2
        capsize=0.01
        ax.plot([x,x], [0.0, tdiv], color='black')
        ax.plot([x-capsize, x+capsize], [0.0,0.0], color='black')
        ax.plot([x-capsize, x+capsize], [tdiv,tdiv], color='black')
        ax.text(0.2, tdiv/2, vbarlabel, 
            horizontalalignment='right', verticalalignment='center', fontsize=10)
    if admix > 0:
        # Drawing an admixture bar
        ax.arrow(0.325, tdiv/4., 0.475-0.325, 0, 
                 head_width=0.01, head_length=0.01, linewidth=1, color='b', length_includes_head=True)
        ax.arrow(0.475, tdiv/4., -(0.475-0.325), 0, 
                 head_width=0.01, head_length=0.01, linewidth=1, color='b', length_includes_head=True)
        ax.text(0.4, -0.005, r'%0.0f %%' % (100*admix), ha='center', va='top', fontsize=10)

    # Defining a vertical bar for the divergence times
    ax.set_xlim(0.1,0.7)
    ax.axis('off')
    


# +
# Using a GridSpec for this plot
fig = plt.figure(figsize=(3.3,6))
grid = plt.GridSpec(3, 3, hspace=0.1, wspace=0.2)

top_left_ax = fig.add_subplot(grid[0, 0])
top_center_ax = fig.add_subplot(grid[0, 1])
top_right_ax = fig.add_subplot(grid[0, 2])

plot_div_model(top_right_ax, tdiv=t1, vbar=False, lw=2, admix=0.025, color='black', solid_capstyle='round')
plot_div_model(top_center_ax, tdiv=t1, vbar=True, vbarlabel = '0.05 ', lw=2, color='black', solid_capstyle='round')
plot_div_model(top_left_ax, tdiv=t2, vbar=True, vbarlabel = r'$T/2N= 0.5$ ', lw=2, color='black', solid_capstyle='round')

top_center_ax.set_title('Recent\n Divergence', fontsize=12)
top_right_ax.set_title('Recent\n Admixture', fontsize=12)
top_left_ax.set_title('Deep\n Divergence', fontsize=12)

# # Plotting GeoDist Results 
bottom_left_ax = fig.add_subplot(grid[1:, 0])
bottom_center_ax = fig.add_subplot(grid[1:, 1])
bottom_right_ax = fig.add_subplot(grid[1:, 2])

cur_geodist1.plot_geodist(bottom_center_ax);
cur_geodist2.plot_geodist(bottom_left_ax);
cur_geodist3.plot_geodist(bottom_right_ax);

bottom_left_ax.set_ylabel(r'Cumulative fraction of variants', fontsize=14);
bottom_right_ax.set_xticklabels([]); bottom_right_ax.set_xticks([]);
bottom_center_ax.set_xticklabels([]); bottom_center_ax.set_xticks([]);
bottom_left_ax.set_xticklabels([]); bottom_left_ax.set_xticks([]);
bottom_center_ax.set_yticklabels([]);
bottom_right_ax.set_yticklabels([]);

plt.savefig(figdir + 'fig4A.pdf', bbox_inches='tight')
