# geodist_rep_paper
Repository to replicate results from the GeoDist paper


## Installation Requirements

We have setup an Anaconda environment to ensure accurate replication of results. You can create the relevant environment by using:

```
conda source activate <test> 
```



## Working from intermediate data

While it is entirely possible to replicate the analyses from the original 30x WGS variant-calls, some of these steps are less idealized for replication on a laptop. To start from a reasonable intermediate datapoint we have the following pre-packaged datasets (allele frequency tables) that have been computed and can be downloaded using the following commands:

```
wget <insert public dropbox link here>
tar -xvf <data.tar.gz>
```

The prepackaged datasets contain: (1) an allele frequency table [more details below] and (2) the relevant snp lists for performing comparisons across the datasets. 

Alternatively, the step of generating Geographic distribution codes for the entire NYGC 1000 Genomes hg38 dataset takes ~40 minutes due to iterating over all of the variants. If you would like to avoid this step and use the pre-specified codes we have generated:

```
wget <insert public dropbox link here>
tar -xvf <data_w_geodist.tar.gz>
```

You can execute either of these steps using `snakemake` as `snakemake download_data` or `snakemake download_data_w_geodist`.


## Generating visual results 

If you have the `geodist` conda environment activated, all you should need to do in order to recreate the various plots :

```
snakemake gen_all_plots 
```

Before doing so include the `--dryrun` flag to visually see all of the various rules run throughout the pipeline.




