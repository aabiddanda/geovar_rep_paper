#!python3

'''
    Scalable way to reformat geodist tables   
'''

import click
import gzip as gz
import numpy as np


@click.command()
@click.option('--freq_file', help='Frequency File from Plink')
@click.option('--pop_labels', help='Population Labels')
@click.option('--header', help='Header or not', type=bool, default=False)
@click.option('--sep', help='Field separator string', type=str, default=None)
@click.option('--mac', help='Minor Allele Count', type=bool, default=True)
def main(freq_file, pop_labels, header, sep, mac):
  """Reformat a plink file into a reasonable frequency table """
  # 1. Getting a set of population labels 
  pops = []
  with open(pop_labels,'r') as pf:
    for line in pf:
      pop = line.split()[0]
      if pop not in pops:
        pops.append(pop)
	# 2. Reading through the variants 
  with gz.open(freq_file, 'rt') as f:
    if header:
      next(f)
    print('\t'.join(['CHR','SNP','A1','A2'] + pops))
    cur_k = None
    cur_pop_vec = ['0' for i in range(len(pops))]
    for line in f:
      fields = line.split(sep)
      k = (fields[0], fields[1], fields[3], fields[4])
      pop = fields[2]
      if k != cur_k:
        if cur_k is None:
          cur_k = k
        else:
          s = '\t'.join(list(cur_k) + cur_pop_vec)
          print(s)
          cur_k = k
          cur_pop_vec = ['0' for i in range(len(pops))]
          i = pops.index(pop)
          if mac:
            cur_pop_vec[i] = '%s' % fields[6]
          else:
            cur_pop_vec[i] = '%s' % fields[7]
      else:
        i = pops.index(pop)
        if mac:
          cur_pop_vec[i] = '%s' % fields[6]
        else:
          cur_pop_vec[i] = '%s' % fields[7]

if __name__ =='__main__':
	main()
