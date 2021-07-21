#! /usr/bin/env bash

#BSUB -J tf_binding
#BSUB -o logs/tf_binding_%J.out
#BSUB -e logs/tf_binding_%J.err
#BSUB -R "select[mem>4] rusage[mem=4]" 
#BSUB -q rna

set -o nounset -o pipefail -o errexit -x

mkdir -p logs

module load meme
module load bedtools

config_file=src/configs/ca_targets_all_tfs_human.yaml
#config_file=src/configs/ca_targets_all_tfs_mouse.yaml
#config_file=src/configs/ca_targets_ca_tfs_human.yaml
#config_file=src/configs/ca_targets_ca_tfs_mouse.yaml
#config_file=src/configs/nfat_targets_all_tfs.yaml
#config_file=src/configs/nfat_targets_nfat_tfs.yaml

# Function to run snakemake
run_snakemake() {
    local num_jobs=$1
    local config_file=$2

    args='
        -o {log}.out 
        -e {log}.err 
        -J {params.job_name} 
        -R "{params.memory} span[hosts=1] " 
        -n {threads} 
        -q rna '

    snakemake \
        --snakefile Snakefile \
        --drmaa "$args" \
        --jobs $num_jobs \
        --latency-wait 60 \
        --rerun-incomplete \
        --configfile $config_file
}

run_snakemake 12 $config_file