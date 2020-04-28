# geodist_rep_paper
Repository to replicate results from the GeoDist paper. 

## Cloning from github

To bring the repository to your local computer, please use `git clone` as follows:

```
git clone https://github.com/aabiddanda/geodist_rep_paper.git
cd geodist_rep_paper
```

## Installation Requirements

We have setup an [Anaconda](https://www.anaconda.com/distribution/) environment to ensure accurate replication of results and management of dependencies. We suggest using this with miniconda. You can create the relevant environment by running:

```
conda env create -f config/env_geodist.yml
conda activate geodist
```

## Working from intermediate data

The pipeline we have written uses the popular workflow managment system, [snakemake](https://snakemake.readthedocs.io/en/stable/). We refer users to the documentation there in order to understand the various rules and dependencies. The step of generating "Geographic distribution Codes" for the entire NYGC 1000 Genomes hg38 dataset takes ~40 minutes due to iterating over all ~92 million variants. If you are interesting in using the same allele frequency binning that we have, we highly suggest downloading an pre-computed dataset below:

```
snakemake download_minimal_data --cores 1 
```

If you are interested in generating the geodist codes from scratch - remove the `data/geodist` subdirectory and then run the command in the following section to regenerate all plots. Be warned that this can take a considerable amount of time and is best done on a HPC cluster (and has only been tested in Linux).

## Generating main plots

If you have the `geodist` conda environment activated, to recreate the main plots you will have to run:

```
snakemake gen_all_plots --cores 1 --dryrun
```

You can remove the `--dryrun` flag to actually run the pipeline. After running the pipeline, you should be able to see the major figures in the `plots` directory as PDFs. Note that these are somewhat different from the versions in the manuscript as they have not been annotated.  


## File Descriptions

### Frequency Files 

The gzipped frequency files are  tab separated files with the following columns:

  * `CHR` : chromosome
  * `SNP` : position
  * `A1` : major allele
  * `A2` : minor allele (globally)
  * `MAC` : global minor allele count
  * `MAF` : global minor allele frequency 

Then the subsequent columns represent the frequency of the globally minor allele (A2) across the defined populations. You can find these in `data/freq` for our minimal dataset.

### GeoDist Files 

The gzipped "GeoDist" files contain relevant frequency information as well as their geographic distribution "Codes" that we use in the manuscript. They have the following fields: 

  * `CHR` : chromosome
  * `SNP` : position
  * `A1` : major allele
  * `A2` : minor allele (globally)
  * `MAC` : global minor allele count
  * `MAF` : global minor allele frequency 
  * `ID` : geographic distribution code (length refers to the number of populations) 

We note that the integers correspond to the "frequency bin" that the variant falls into within that population. For further detail on the particular scheme used to bin variants please find the details in our paper:

TBD

## Questions

For any questions on this pipeline please either raise an issue or email Arjun Biddanda <abiddanda[at]uchicago.edu>.
