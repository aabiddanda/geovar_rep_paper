# geodist_rep_paper
Repository to replicate results from the GeoDist paper

## Installation Requirements

We have setup an [Anaconda](https://www.anaconda.com/distribution/) environment to ensure accurate replication of results and management of dependencies. We suggest using this with miniconda. You can create the relevant environment by running:

```
conda env create -f config/env_geodist.yml
conda activate geodist
```

## Working from intermediate data

The step of generating "Geographic distribution Codes" for the entire NYGC 1000 Genomes hg38 dataset takes ~40 minutes due to iterating over all of the variants. If you are interesting in using the same allele frequency binning that we have, we highly suggest downloading an pre-computed dataset below:

```
snakemake download_data_w_geodist --cores <number of cores> 
```

If you are interested in generating the geodist codes from scratch - remove the `data/geodist` subdirectory and then run the following command to regenerate all plots. Be warned that this can take a considerable amount of time and is best done on a HPC cluster.

## Generating visual results 

If you have the `geodist` conda environment activated, to recreate the main plots you will have to run:

```
snakemake gen_all_plots --cores <number of cores> --dryrun 
```

You can remove the `--dryrun` flag to actually run the pipeline. After running the pipeline, you should be able to see the major figures in the `plots` directory as PDFs. Note that these are somewhat different from the versions in the manuscript as they have not been annotated.  

## Questions

For any questions on this pipeline please either raise an issue or email Arjun Biddanda <abiddanda[at]uchicago.edu>.
