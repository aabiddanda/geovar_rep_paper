#!/usr/local/bin/python3

'''
 Generate a frequency table for use in downstream GeoDist applications
'''

from __future__ import print_function
import sys

import click
import numpy as np 
import pandas as pd
import gzip as gz
import json

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

@click.command()
@click.option('--mac', required=True, help='Minor allele count file')
@click.option('--total', required=True, help='Total samples count file')
def main(mac, total):
    #1. iterating through the files simultaneously
    with gz.open(mac, 'rt') as macfile, gz.open(total, 'rt') as total_file:
      i = 0
      for mac_line, total_line in zip(macfile,total_file):
        if i == 0:
          spltln = mac_line.split()
          header_list = spltln[0:4] + ['MAC','MAF'] + spltln[4:]
          print('\t'.join(header_list))
        else:
          # Splitting each of the lines to compute the allele frequency
          spltlnmac = mac_line.split()
          spltlntotal = total_line.split()
          macs = np.array(spltlnmac[4:], dtype=np.float32)
          totals = np.array(spltlntotal[4:], dtype=np.float32)
          allele_freqs = np.nan_to_num(macs/totals)
          global_maf = np.sum(macs)/np.sum(totals)
          total_mac = np.int64(np.sum(macs))
          total_list = spltlnmac[0:4] + [str(total_mac), str(global_maf)] + [str(f) for f in allele_freqs]
          print('\t'.join(total_list))
        if i % 100000 == 0:
          eprint('%d00 K snps finished' % int(i/100000))
        i += 1

if __name__ =='__main__':
  main()
  
