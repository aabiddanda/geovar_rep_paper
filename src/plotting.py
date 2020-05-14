import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib as mpl

plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['figure.facecolor'] = "w"
plt.rcParams['figure.autolayout'] = True

import click

# Hex-Colors for GeoDist Mappings
# Blues_HSL = [188,45,81]
# Blues_Hex [#B9DFE4, #061784, #FFFFFF]

# cmap_test = mpl.colors.LinearSegmentedColormap.from_list('Custom Blues', ['#B9DFE4', '#061784', '#FFFFFF'], 3)

# # define the bins and normalize
# bounds = np.linspace(0, 20, 21)
# norm = mpl.colors.BoundaryNorm(bounds, cmap.N)


class GeoDistPlot:

  def __init__(self):
    """ Initializing plotting object """
    # Defining general parameters
    self.fontsize = 14
    self.title_fontsize = 12
    self.h_orient ='center'
    self.v_orient = 'center'
    self.bar_color = 'gray'
    self.border_color = 'gray'
    self.line_weight = 0.5
    self.alpha = 0.5
    self.x_lbl_fontsize = 16
    self.y_lbl_fontsize = 12
    # Original Data for GeoDist
    self.orig_missing = None
    self.orig_geodist = None
    self.orig_ngeodist = None
    self.orig_fgeodist = None
    self.orig_npops = None
    self.orig_ncat = None
    # Actual plotting data for GeoDist
    self.geodist = None
    self.ngeodist = None
    self.fgeodist = None
    self.npops = None
    self.poplist = None
    self.ncat = None
    # Colormap parameters
    self.colors = None
    self.str_labels = None
    self.lbl_colors = None

  def __str__(self):
    """ print all active parameters for plotting """
    test_str = 'Fontsize : %d\n' % self.fontsize
    test_str += 'Bar Color : %s\n' % self.bar_color
    test_str += 'Border Color : %s\n' % self.border_color
    test_str += 'Line Weight : %0.2f\n' % self.line_weight
    test_str += 'Alpha : %0.2f\n' % self.alpha
    # Printing data level parameters
    if self.ngeodist is not None:
      test_str += 'Number of SNPS : %d\n' % np.sum(self.ngeodist)
      test_str += 'Number of Populations : %d\n' % self.npops
      test_str += 'Number of Categories : %d\n' % self.ncat
    else:
      test_str += 'Number of SNPS : 0\n'
      test_str += 'Number of Populations : NA\n'
      test_str += 'Number of Categories : NA\n'

    return(test_str)

  def _add_text_data(self, inputfile, filt_unobserved = True):
    """ Add/replace data for a GeoDistPlot object """
    df = np.loadtxt(inputfile, dtype=str)
    geodist = df[:,0]
    ngeodist = df[:,1].astype(int)
    assert(geodist.size == ngeodist.size)
    npops = np.array([len(x) for x in geodist])
    ncat = np.array([max(list(x)) for x in geodist], dtype=int)
    # Testing out a filterign 
    if filt_unobserved:
      unobserved_geodist = '0' * npops[0]
      idx = np.where(geodist == unobserved_geodist)[0]
      if idx.size >= 1:
        self.orig_missing = ngeodist[idx][0]
        geodist = np.delete(geodist, idx)
        ngeodist = np.delete(ngeodist, idx)
      else:
        raise ValueError('No excluded variants!')
    assert(len(np.unique(npops)) <= 1)
    self.orig_npops = npops[0]
    self.orig_geodist = geodist
    self.orig_ngeodist = ngeodist
    self.orig_ncat = max(ncat) + 1
    self.orig_fgeodist = (self.orig_ngeodist / np.sum(self.orig_ngeodist))
    # Setting all the plotting variables
    self.npops = npops[0]
    self.geodist = geodist
    self.ngeodist = ngeodist
    self.ncat = max(ncat) + 1
    self.fgeodist = self.orig_fgeodist

  def _add_data_jsfs(self, jsfs):
    """ Add data from a joint SFS via numpy   """
    npops = len(jsfs.shape)
    ncats = np.array(jsfs.shape)
    # Assert that all of the populations have the same categories
    assert(np.all(ncats == ncats[0]))
    # iterate through them all...
    geo_codes = []
    ngeo_codes = []
    for i,x in np.ndenumerate(jsfs):
      # generate category
      cat = ''.join([str(v) for v in list(i)])
      geo_codes.append(cat)
      ngeo_codes.append(x)
    geodist = np.array(geo_codes)
    ngeodist = np.array(ngeo_codes)
    assert(geodist.size == ngeodist.size)
    self.orig_npops = npops
    self.orig_geodist = geodist
    self.orig_ngeodist = ngeodist
    self.orig_ncat = ncats[0] + 1
    self.orig_fgeodist = (self.orig_ngeodist / np.sum(self.orig_ngeodist))
    # Setting all the plotting variables
    self.npops = npops
    self.geodist = geodist
    self.ngeodist = ngeodist
    self.ncat = ncats[0] + 1
    self.fgeodist = self.orig_fgeodist

  def _filter_data(self, max_freq=0.005, rare=False):
    """ Filter geodist data for easier plotting """
    assert(self.orig_geodist is not None)
    assert(self.orig_ngeodist is not None)
    assert(self.orig_fgeodist is not None)
    if rare:
      idx = np.where(self.orig_fgeodist < max_freq)[0]
    else:
      idx = np.where(self.orig_fgeodist >= max_freq)[0]
    assert(idx[0].size != 0)
    self.orig_fgeodist_alt = self.orig_fgeodist[idx]
    self.ngeodist = self.orig_ngeodist[idx]
    self.geodist = self.orig_geodist[idx]
    self.fgeodist = self.ngeodist / np.sum(self.ngeodist)
    # sort to try to order based on frequency (high to low)
    sorted_idx = np.argsort(self.fgeodist)[::-1]
    self.ngeodist = self.ngeodist[sorted_idx]
    self.geodist = self.geodist[sorted_idx]
    self.fgeodist = self.fgeodist[sorted_idx]
    self.orig_fgeodist_alt = self.orig_fgeodist_alt[sorted_idx]

  def _add_cmap(self, base_cmap='Blues',
                str_labels=['U', 'R', 'C'],
                lbl_colors=['black', 'black', 'white']):
    """ Create a colormap object to use """
    assert(self.ncat is not None)
    assert(len(str_labels) == len(lbl_colors))
    assert(len(str_labels) == self.ncat)
    # Generating a discrete mapping
    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, self.ncat))
    cmap_name = base.name + str(self.ncat)
    test_cmap = base.from_list(cmap_name, color_list, self.ncat)
    # normalizing to 1
    norm = mpl.colors.Normalize(vmin=0, vmax=self.ncat)
    colors = [mpl.colors.rgb2hex(test_cmap(norm(i))) for i in range(self.ncat)]
    self.colors = colors
    self.str_labels = str_labels
    self.lbl_colors = lbl_colors

  def _set_colors(self, colors):
    """ Method to add custom hex colors during plotting """
    assert(self.ncat is not None)
    assert(self.colors is not None)
    assert(len(colors) == len(self.colors))
    self.colors = colors

  def _add_poplabels(self, popfile):
    """ Adding population labels here """
    assert(self.geodist is not None)
    assert(self.ngeodist is not None)
    assert(self.npops is not None)
    assert(self.ncat is not None)
    # Reading the appropriate population file
    pops = np.loadtxt(popfile, dtype=str)
    assert(pops.size == self.npops)
    self.poplist = pops
  
  def _reorder_pops(self, new_poplist):
    """ Reordering populations within a geodist instance """
    assert(new_poplist.size == self.poplist.size)
    new_pop_idx = np.hstack([np.where(self.poplist == i)[0] for i in new_poplist])
    acc = []
    for j in range(self.orig_ngeodist.size):
        acc.append(''.join([self.orig_geodist[j][i] for i in new_pop_idx]))
    new_geodist = np.array(acc)
    self.orig_geodist = new_geodist
    self.poplist = new_poplist
    
  def _add_poplabels_manual(self, poplabels):
    """ Adding population labels here """
    assert(self.geodist is not None)
    assert(self.ngeodist is not None)
    assert(self.npops is not None)
    assert(self.ncat is not None)
    assert(poplabels.size == self.npops)
    self.poplist = poplabels

  def plot_geodist(self, ax, superpops=None, superpop_lbls=None):
    """ Final function to call to generate a geodist plot on a particular axis """
    # Starting assertions to make sure we can call this
    assert(self.geodist is not None)
    assert(self.ngeodist is not None)
    assert(self.npops is not None)
    assert(self.poplist is not None)
    assert(self.ncat is not None)
    assert(self.colors is not None)
    assert(self.str_labels is not None)
    assert(self.lbl_colors is not None)
    assert(self.poplist is not None)
    # Setting up the codes here
    x_limits = ax.get_xlim()
    xbar_pts = np.linspace(x_limits[0], x_limits[1], num=self.npops + 1)
    xpts_shifted = (xbar_pts[1:] + xbar_pts[:-1])/2.0
    ax.set_xticks(xpts_shifted)
    # setting up the vertical lines
    for x in xbar_pts:
      ax.axvline(x=x, color=self.bar_color, lw=self.line_weight, alpha=self.alpha);
    
    if superpops is not None:
      for i in superpops:
        ax.axvline(x = xbar_pts[i], color='gray', lw=1.0)
    
    # changing the border color
    for spine in ax.spines.values():
      spine.set_edgecolor(self.border_color);

    # Plotting the horizontal bars in this case
    ylims = ax.get_ylim()
    y_pts = np.cumsum(self.fgeodist) / (ylims[1] - ylims[0])
    n_y = y_pts.size
    cur_y = ylims[0]
    nsnps = np.sum(self.orig_ngeodist)
    for i in range(n_y):
      y = y_pts[i]
      y_dist = y - cur_y
      ax.axhline(y = y, color=self.bar_color, lw=self.line_weight, alpha=self.alpha);
      cur_code = list(self.geodist[i])
      for j in range(self.npops):
        # Defining the current category
        cur_cat = int(cur_code[j])
        # Drawing in the rectangle
        cur_xy = xbar_pts[j], cur_y
        rect = patches.Rectangle(xy = cur_xy,
                                 width=(xbar_pts[j+1] - cur_xy[0]),
                                 height=y_dist,
                                 facecolor=self.colors[cur_cat]);
        ax.add_patch(rect);
        # Drawing in the text
        y_frac = y_dist/(ylims[1] - ylims[0])
        fontscale = 1.0
        freq_thresh = 0.06
        if y_frac < freq_thresh:
          fontscale = y_frac*fontscale/freq_thresh
        ax.text(x=xpts_shifted[j],
                y=(cur_y + y_dist/2.0),
                ha=self.h_orient, va=self.v_orient,
                s=self.str_labels[cur_cat],
                color=self.lbl_colors[cur_cat],
                fontsize=self.fontsize*fontscale);
      cur_y = y
    ax.set_ylim(ylims)
    return(ax, nsnps, y_pts)
  
  def plot_percentages(self, ax):
    """ Generates a plot with the percentages  """
    ns = self.ngeodist
    fracs = ns / np.sum(ns)
    cum_frac = np.cumsum(self.fgeodist)
    # Setting the border here
    for spine in ax.spines.values():
      spine.set_edgecolor(self.border_color);
    prev = 0.0
    for i in range(cum_frac.size):
      ax.axhline(cum_frac[i], lw=0.5, color=self.bar_color)
      # Get the midpoint 
      ydist = (cum_frac[i] - prev)/2.
      fontscale = min(1., self.fontsize*fracs[i])
      nstr = '{:,}'.format(ns[i])
      ax.text(x=0.01, y=prev+ydist, s = '%s (%d%%)' % (nstr,round(self.orig_fgeodist_alt[i]*100)), va='center', fontsize=self.fontsize*fontscale)
      prev = cum_frac[i]
    ax.set_ylim(0,1)
    return(ax)
      

    
def plot_multiple_geodist(geodist_obj_list, subsets, xsize=1, ysize=4, hwidth=0.1, top_buff=0.5, bot_buff = 0.5, left_buff = 0.75, ylabel='Cumulative fraction of variants', superpops=None, superpop_lbls=None):
    """
      Function to plot a list of geodist objects 
      NOTE : this is a preliminary function currently 
    """
    assert(len(geodist_obj_list) == len(subsets))
    n = len(geodist_obj_list)
    
    fig_width = left_buff + xsize*n + hwidth*n 
    fig_height = bot_buff + ysize + top_buff
    
    # start the figure
    fig = plt.figure(constrained_layout=True, figsize=(fig_width,fig_height))
    
    # Setting the proportional arguments
    left_start = left_buff / fig_width
    left_increment = (fig_width - left_buff)/fig_width / n
    bot_prop = bot_buff/fig_height
    w_prop = left_increment-(hwidth/fig_width)
    top_prop = ysize/fig_height
    
    # Setting all the 
    axs = [fig.add_axes([left_start + i*left_increment, bot_prop, w_prop, top_prop]) for i in range(n)]
    [a.set_yticks([]) for a in axs[1:]];
    
    for i in range(n):
        subset = subsets[i]
        cur_geodist = geodist_obj_list[i]
        cur_geodist.fontsize = 14
        _, nsnps, _ = cur_geodist.plot_geodist(ax=axs[i], superpops=superpops, superpop_lbls=superpop_lbls)
        axs[i].set_xticklabels(cur_geodist.poplist, fontsize=10, rotation=90, ha='center')
        if cur_geodist.orig_missing is not None:
          
          n = np.sum(cur_geodist.orig_ngeodist) + cur_geodist.orig_missing
          nstr = '{:,}'.format(n)
          nmiss_str = '{:,}'.format(cur_geodist.orig_missing)
          axs[i].set_title('%s\n $S$ = %s\n$S_u$ = %s (%d%%)' % (subset, nstr, nmiss_str, int(cur_geodist.orig_missing/n*100)), fontsize=10)
        else:
          n = np.sum(cur_geodist.orig_ngeodist)
          nstr = '{:,}'.format(n)
          axs[i].set_title('%s\n S = %s' % (subset,nstr), fontsize=10)
    axs[0].set_ylabel(ylabel, fontsize=14)
    return(fig, axs)    

def plot_geodist_w_percentages(geodist_obj, subset="", ylabel='Cumulative fraction of variants', xsize=1, ysize=4, hwidth=0.1, top_buff=0.5, bot_buff = 0.5, left_buff = 0.75, superpops=None, superpop_lbls=None):
    """
      Function to plot a list of geodist objects 
      NOTE : this is a preliminary function currently 
    """
    n = 2
    
    fig_width = left_buff + xsize*n + hwidth*n 
    fig_height = bot_buff + ysize + top_buff
    
    # start the figure
    fig = plt.figure(constrained_layout=True, figsize=(fig_width,fig_height))
    
    # Setting the proportional arguments
    left_start = left_buff / fig_width
    left_increment = (fig_width - left_buff)/fig_width / n
    bot_prop = bot_buff/fig_height
    w_prop = left_increment-(hwidth/fig_width)
    top_prop = ysize/fig_height
    
    # Setting all the 
    axs = [fig.add_axes([left_start + i*left_increment, bot_prop, w_prop, top_prop]) for i in range(n)]
    [a.set_yticks([]) for a in axs[1:]];
    

    geodist_obj.fontsize = 14
    _, nsnps, _ = geodist_obj.plot_geodist(ax=axs[0], superpops=superpops, superpop_lbls=superpop_lbls)
    axs[0].set_xticklabels(geodist_obj.poplist, fontsize=10, rotation=90, ha='center')
    n = np.sum(geodist_obj.orig_ngeodist)
    nstr = '{:,}'.format(n)
    #  Setting the title figures
    axs[0].text(1.0 + (hwidth/fig_width), 1.01, '%s\n (S = %s)' % (subset,nstr), 
                va='bottom', ha='center', fontsize=14)
    geodist_obj.fontsize=12
    geodist_obj.plot_percentages(axs[1])
    axs[1].set_xticks([]);
    axs[0].set_ylabel(ylabel, fontsize=14)
    return(fig, axs)  
  
  
    
@click.command()
@click.option('--geodist_cnt', help='Geodist Categories File', required=True)
@click.option('--pop_panel', help='Cluster Labels', required=True)
@click.option('--max_freq', help='Minimal global minor allele frequency',type=float, default=0.001)
@click.option('--rare', help='Upper or lower bound on variant filtering', type=bool, default=False)
@click.option('--cmap_name', default='Blues')
@click.option('--figsz_x', help='Figure width', type=float, default = 4)
@click.option('--figsz_y', help='Figure width', type=float, default = 8)
@click.option('--title', help='title', default = None)
@click.option('--dpi', help='Output resolution', default=300)
@click.option('--out', required=True, help='Output File (PNG, JPG, PDF)')
def main(geodist_cnt, pop_panel, max_freq, rare, cmap_name, figsz_x, figsz_y, title, dpi, out):
    #1.Creating the geodist object
    cur_geodist = GeoDistPlot()
    cur_geodist._add_data(geodist_cnt)
    cur_geodist._add_poplabels(pop_panel)
    cur_geodist._filter_data(max_freq=max_freq, rare=rare)
    # TODO : this is what
    if cur_geodist.ncat == 3:
      cur_geodist._add_cmap(base_cmap=cmap_name)
    else:
      cur_geodist._add_cmap(base_cmap=cmap_name,
                            str_labels=['U','R','L','C'],
                            lbl_colors=['black','black','white','white'])
    # Setting up the plot
    f, ax = plt.subplots(1,2,figsize=(figsz_x*2, figsz_y))
    _, nsnps, _ = cur_geodist.plot_geodist(ax[0])
    ax[0].set_yticks([0, 0.25,0.5,0.75,1.0])
    ax[0].set_yticklabels(['0.0', '0.25', '0.5', '0.75', '1.0'],fontsize=cur_geodist.y_lbl_fontsize)
    ax[0].set_xticklabels(cur_geodist.poplist, fontsize=cur_geodist.x_lbl_fontsize, rotation=90)
    ax[0].set_xlabel('Pop', fontsize=cur_geodist.x_lbl_fontsize)
    # NOTE : this scaling is very adhoc...
    cur_geodist.fontsize = 0.65*cur_geodist.fontsize
    _ = cur_geodist.plot_percentages(ax[1])
    ax[1].set_yticks([]);
    ax[1].set_xticks([]);
    ax[1].set_ylim(0,1);
    f.suptitle('%s\n N = %d' % (title,nsnps), fontsize=cur_geodist.title_fontsize, y=1.05)
    plt.tight_layout()
    plt.savefig(out, dpi=dpi, bbox_inches='tight')
    plt.close('all')


if __name__ == '__main__':
  main()
