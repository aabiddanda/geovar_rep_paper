# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
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

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os

# setting Font to Arial
import matplotlib as mpl
mpl.rcParams['font.family'] = "Arial"
mpl.rcParams['font.sans-serif'] = "Arial"
mpl.rcParams['font.size'] = "12"
mpl.font_manager._rebuild()


def no_spines(ax):
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    return ax


def wf_sim(N, t_max, n_traj, f0=None):
    '''
    Simulate n_traj neutral Wright-Fisher frequency trajectories, starting at frequency f0.
    If f0==None, simulate new mutations with a reflective lower bound at zero.
    '''
    f = np.zeros((t_max, n_traj))
    if f0:
        f[0] = f0
    for i in range(1, t_max):
        if f0:
            f[i] = np.random.binomial(N, f[i-1]) / N
        else:
            f[i] = np.clip(np.random.binomial(N, f[i-1]), 1, N) / N
    return f


def bounds(f0, t):
    '''Bounds on the typical frequency evolution starting at frequency f0.'''
    if f0 == 0:
        return (0, t)
    else:
        envelope = np.sqrt(f0*(1-f0)*t)
        return (np.clip(f0 - envelope,0,1), np.clip(f0 + envelope,0,1))


def p_ext(f0, t):
    '''Approximate extinction probability starting at frequency f0'''
    return np.piecewise(t, [t <= 0, t > 0], [0, lambda x: np.exp(-f0/x)])


# ## Run simulations

f_common = 0.25
N = 1000
t_max = 0.5
gens = int(N*t_max)
t = np.arange(0, gens)/N
n_traj = 5

np.random.seed(101)
traj_common = wf_sim(N, gens, n_traj, f_common)
plt.plot(t, traj_common)

np.random.seed(102)
traj_new = wf_sim(N, gens, n_traj)
plt.plot(t, traj_new)


# ## Test subplots

def plot_envelopes(ax, alpha=0.15, downsample=4, lw=1.5):
    c_new = 'C1'
    c_common = 'C0'
    
    ax.fill_between(t, *bounds(0,t), alpha=alpha, color=c_new)
    ax.plot(t[::downsample], traj_new[::downsample], color=c_new, lw=lw)

    ax.fill_between(t, *bounds(f_common,t), alpha=alpha, color=c_common)
    ax.plot(t[::downsample], traj_common[::downsample], color=c_common, lw=lw)
    
    ax.set_ylabel('Allele frequency')
    ax.set_xlabel('Scaled time, $t/2N$')
    
    return ax


# +
fig = plt.figure(figsize=(3,3))
ax = fig.add_subplot(111)
plot_envelopes(ax, lw=1.5, alpha=0.1)

ax.set_xlim([0, 0.51]) 
ax.set_xticks([0, 0.25, 0.5])


# -

def plot_extinction(ax, lw=1.5):
    ax.plot(t, p_ext(f_common, t))
    ax.set_ylabel('Extinction prob.')
    ax.set_xlabel('Scaled time, $T/2N$')
    return ax



fig = plt.figure(figsize=(3,1.5))
ax = fig.add_subplot(111)
plot_extinction(ax)
ax.set_xlim([0, 0.51]) 
ax.set_xticks([0, 0.25, 0.5])

# ## Full plot

# +
fig = plt.figure(figsize=(2.75, 4.5))

lw = 1.5

spec = gridspec.GridSpec(ncols=1, nrows=2, figure=fig, height_ratios=[2,1])
ax1 = fig.add_subplot(spec[0])
ax2 = fig.add_subplot(spec[1])

plot_envelopes(ax1, lw=lw)
ax1.set_xticklabels([])
ax1.set_xlabel(None)
ax1.set_ylim([0, 0.5])
ax1.set_yticks(np.linspace(0, 0.5, 3))


plot_extinction(ax2, lw=lw)
ax2.set_ylim([0, 0.7])
ax2.set_yticks(np.linspace(0,0.7,2))
ax2.set_yticklabels(map(lambda x: f"{x:.2f}", ax2.get_yticks()))

for ax in fig.axes:
    ax.vlines([0.05,0.5], 0, 1, linestyle='dashed', lw=lw, zorder=3)
    ax.set_xlim([0, 0.51]) 
    ax.set_xticks([0.05,0.5])
    ax.xaxis.set_tick_params(width=lw)
    no_spines(ax)

figdir = '../plots/figure4/'
os.makedirs(figdir, exist_ok=True)
plt.savefig(figdir + 'fig4BC.pdf', bbox_inches='tight')
# -


