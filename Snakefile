""" Snake pipeline for running tf analysis """


# Configure shell for all rules 
shell.executable("/bin/bash")
shell.prefix("set -o nounset -o pipefail -o errexit -x; ")

BED_FILE = config["BED_FILE"]
FASTA_FILE = config["FASTA_FILE"]
GENES = config["GENES"]
RESULTS = config["RESULTS"]
DISTANCES = config["DISTANCES"]
CHR_SIZE = config["CHR_SIZE"]
MOTIFS = config["MOTIFS"]
SPECIES = config["SPECIES"]


rule all:
    input: 
        expand(
            "{results}/meme_{gene_list}/{distance}_test/subset_promoter.fa",
            results = RESULTS, gene_list = GENES, distance = DISTANCES
            ),

        expand(
            "{results}/meme_{gene_list}/{distance}_test/fimo_out/",
            results = RESULTS, gene_list = GENES, distance = DISTANCES
            ),
        expand(
            "{results}/homer_{gene_list}/{distance}_test/",
            results = RESULTS, gene_list = GENES, distance = DISTANCES
            )

include: "src/rules/run_meme.snake"
include: "src/rules/run_homer.snake"
