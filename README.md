# geodist_rep_paper
Repository to replicate results from the GeoDist paper


## Installation Requirements

We have setup an Anaconda environment to ensure accurate replication of results. You can create the relevant environment by using 



## Working from intermediate data

While it is entirely possible to replicate the analyses from the original 30x WGS variant-calls, some of these steps are less idealized for replication on a laptop. To start from a reasonable intermediate datapoint we have the following pre-packaged datasets (allele frequency tables) that have been computed and can be downloaded with the following commands:

```
wget <insert public dropbox link here>
tar -xvf <data.tar.gz>
```

The prepackaged datasets contain: (1) an allele frequency table [more details below] and (2) the relevant snp lists for performing comparisons across the datasets. 


## Generating Results 

If you have the geodist conda environment activated, all you should need to do in order to recreate the various plots should be to run:

```
snakemake gen_all_plots 
```

Before doing so include the `--dryrun` flag to visually see all of the various rules run within the pipeline.


