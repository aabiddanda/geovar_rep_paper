# geodist_rep_paper
Repository to replicate results from the GeoDist paper


## Installation Requirements

We have setup an [Anaconda](https://www.anaconda.com/distribution/) environment to ensure accurate replication of results and management of dependencies. You can create the relevant environment by running:

```
conda env create -f config/environment.yml
source activate geodist
```

## Working from intermediate data

While it is entirely possible to replicate the analyses from the original 30x WGS variant-calls, some of these steps are less idealized for replication on a laptop. To start from a reasonable intermediate datapoint we have the following pre-packaged datasets (allele frequency tables) that have been computed and can be downloaded using the following commands:

```
wget <dropbox link>
tar -xvf <data.tar.gz>
```

The prepackaged datasets contain: (1) an allele frequency table [more details below] and (2) the relevant SNP lists for comparing across datasets. 

Alternatively, the step of generating "Geographic distribution Codes" for the entire NYGC 1000 Genomes hg38 dataset takes ~40 minutes due to iterating over all of the variants. If you are interesting in using the same allele frequency binning that we have, we highly suggest downloading this dataset. If you would like to avoid this step and use the pre-specified codes we have generated in our analysis:

```
wget <dropbox link here>
tar -xvf <data_w_geodist.tar.gz>
```

You can execute either of these steps using `snakemake` as `snakemake download_data` or `snakemake download_data_w_geodist`.

## Generating visual results 

If you have the `geodist` conda environment activated, all you should need to do in order to recreate the various plots :

```
snakemake gen_all_plots --cores <number of cores> 
```

Before doing so include the `--dryrun` flag to visually see all of the various rules run throughout the pipeline. After running the pipeline, you should be able to see the major figures in the `plots` directory as PDFs. Note that these are somewhat different from the versions in the manuscript as they have not been annotated.  

## Questions

For any questions on this pipeline please either raise an issue or email Arjun Biddanda <abiddanda[at]uchicago.edu>.
