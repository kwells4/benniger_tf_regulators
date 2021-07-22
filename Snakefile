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
# OTHER_SPECIES = config["OTHER_SPECIES"]

# OTHER_SPECIES[SPECIES] = {"bed":BED_FILE, "fasta":FASTA_FILE, "chr_size":CHR_SIZE}

# bed_file = OTHER_SPECIES["hg38"]["bed"]

# # Create gene list
# GENE_LIST = []
# for gene_set in GENES:
#     with open(GENES[gene_set][0], "r") as gene_file:
#         for line in gene_file:
#             GENE_LIST.append(line.strip())

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
            ),
        # expand(
        #     "{results}/conservation_{gene_list}/{distance}_test/{gene}_species_promoter_fasta.fa",
        #     results = RESULTS, gene_list = GENES, distance = DISTANCES, gene = GENE_LIST
        #     ),
        # expand(
        #     "{results}/conservation_{gene_list}/{distance}_test/{gene}_alignment_fasta.fa",
        #     results = RESULTS, gene_list = GENES, distance = DISTANCES, gene = GENE_LIST
        #     )

include: "src/rules/run_meme.snake"
include: "src/rules/run_homer.snake"
# include: "src/rules/generate_fastas.snake"
