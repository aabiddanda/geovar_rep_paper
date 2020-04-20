#!/usr/local/bin/python3

'''
  Naive categorization of geodist categories
'''

import sys
import click
import numpy as np
import gzip as gz
from tqdm import tqdm


COMP_BASES = {'A': 'G', 'G': 'A', 'C': 'T', 'T': 'C'}

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def generate_bins(endpts):
    assert(np.all(np.array(endpts) < 1.0))
    b = 0.0
    y = len(endpts)
    bins = []
    for x in endpts:
      bins.append((b, x))
      b = x
    bins.append((b, 1.0))
    return(bins)

def geodist_binning(X, bins=[(0, 0), (0, 0.01), (0.01, 0.05), (0.05, 1.0)]):
    """
      Calculate the **simple** version of Geodist labelings
    """
    output = np.zeros(shape=X.shape, dtype=int)
    i = 1
    for b in bins[1:]:
      idx = np.where((X > b[0]) & (X <= b[1]))
      output[idx] = i
      i += 1
    return(output)


def flip_af(x, maj, anc):
    """
    Flip allele frequency vector if REF and ANC don't line up
    """
    maj = maj.rstrip()
    anc = anc.rstrip()
    # defining a dictionary for the complementary bases
    x_new = x
    if (maj != anc) and (maj != COMP_BASES[anc]):
        x_new = 1. - x
    return(x_new)

@click.command()
@click.option('--freqs', required=True, help='Minor allele count file')
@click.option('--bins', required=False,
              type=str, default="[0.0,  0.05]",
              help='Bin endpoints for allele frequency')
@click.option('--anc', required=False, type=bool, default=False, help='Using the ancestral allele!')
@click.option('--ancq', required=False, type=bool, default=True, help='Required quality for the ancestral calls')
@click.option('--strand', required=False, type=bool, default=False, help='Remove strand-ambiguous alleles')
def main(freqs, bins, anc, ancq, strand):
    # 1. Generating the bins
    bin_endpts = list(map(float, bins.strip('[]').split(',')))
    bins = generate_bins(bin_endpts)
    # 2. Reading through the frequency file line by line and printing it out
    # 2a. First setting up the parameters
    idx = 6
    if anc:
        idx = 7
    # streaming through the file now
    with gz.open(freqs, 'rt') as f:
      header_ln = f.readline().split()
      print('\t'.join(header_ln[0:idx] + ['ID']))
      
    with gz.open(freqs, 'rt') as f:
        next(f)
        for line in tqdm(f):
            spltln = line.split()
            freq_values = np.array(spltln[idx:], dtype=np.float32)
            test_values = geodist_binning(freq_values, bins=bins)
            code = ''.join([str(i) for i in test_values])
            print('\t'.join(spltln[0:idx] + [code]))
#             i += 1    
#             if i % 100000 == 0:
#                 eprint("Finished %d00 k Variants" % int(i/100000))

if __name__ =='__main__':
    main()
    
    
#                 if anc:
#                     major_allele = spltln[3]
#                     minor_allele = spltln[2]
#                     anc_allele = spltln[4]
#                     if strand and (major_allele == COMP_BASES[minor_allele]):
#                         eprint("Filtering out variant", major_allele, minor_allele)
#                         continue
#                     elif anc_allele in ['.', '-', 'N']:
#                         continue
#                     elif ancq and (not anc_allele.isupper()):
#                         # We don't have the necessary quality
#                         continue
#                     freq_values = flip_af(freq_values, major_allele,
#                                 anc_allele.upper())
#                     test_values = simple_geodist_binning(freq_values, bins=bins)
#                     code = ''.join([str(i) for i in test_values])
#                     print('\t'.join(spltln[0:idx] + [code]))