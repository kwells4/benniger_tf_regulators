# Analysis of TF regulators

A snakemake pipeline and R scripts to replicate the analysis in Dynamic changes in β-cell [Ca 2+ ] regulates NFAT activation, gene transcription and islet gap junction communication published in Molecular Metabolism DOI: 10.1016/j.molmet.2021.101430

This pipeline first uses homer to look for enriched motifs around the provided gene list. It next runs Meme to identify the location of specific TF motifs. To run Meme, first a bed file is generated that includes regions that are certain distances away from the promoter (specified in "DISTANCES" in the config file). This bed file is then used to generate a fasta file that is specifically the regions around the promoters of interest. Finally, fimo is run looking for the motifs specified in "MOTIFS" in the config file

There are many config files found in src/configs that were used to generate the plots found in the manuscript. You can also run this pipeline with your own gene lists by following the below instructions.

Writen by Kristen Wells 2021

To use:

1. Download and install miniconda3: For Linux
```{bash}
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh bash Miniconda3-latest-Linux-x86_64.sh
```
2. Install Snakemake:
```{bash}
conda install -c conda-forge mamba
mamba install -c bioconda -c conda-forge snakemake
```

3. Update the config file (config.yaml) 
>* RESULTS: The path to the output directory. It must already exist, you can make it with
```{bash}
mkdir results
```
>* BED_FILE: Bath to the referece bed file
>* FASTA_FILE: The path to the genome fasta file
>* CHR_SIZE: The path to a file containing the size of all chromosomes
>* MOTIFS: A path to jaspar style motifs to use
>* SPECIES: What species you want to analyze
>* DISTANCES: A list of distances from the promoter to check with homer and meme.
>* GENES: A list of named gene lists to test. Multiple lists can be provided.

4. Update snakecharmer.sh to your specific cluster specs. 
>* change the -q argument to the queue you want to use 

5. submit the job using `bsub < snakecharmer.sh`
