library(tidyverse)
library(biomaRt)

working_dir <- 
  "/Users/wellskr/Documents/Analysis/Richard_Benniger/benniger_tf_regulators/results/"

gene_list <- "nfat_targets_ca_tfs"

meme_dir <- "meme_NFAT_LIST"

conservation_dir <- "UCSC_phast"

conservation_path <- paste0(working_dir, gene_list, "/",
                            conservation_dir)

distances <- c("300_50_test", "2000_2000_test", "2000_0_test", "5000_5000_test")

# Don't have conservation information here
skip_genes <- c("ENSMUSG00000005232", "ENSMUSG00000050931")

# Functions --------------------------------------------------------------------
get_conservation <- function(fimo_row, conservation_path,
                             skip_rows = 9){
  gene <- fimo_row$gene_id
  conservation_file <- read.table(paste0(conservation_path, "/",
                                         gene), skip = skip_rows,
                                  header = FALSE)
  colnames(conservation_file) <- c("position", "score")
  
  motif_positions <- fimo_row$motif_start:fimo_row$motif_end
  
  conservation_short <- conservation_file %>%
    dplyr::filter(position %in% dplyr::all_of(motif_positions))
  
  average_conservation <- mean(conservation_short$score)
  
  all_conservation <- paste(conservation_short$score,
                            collapse = ",")
  
  return(data.frame(average_conservation = average_conservation,
                    all_conservation = all_conservation))
}

# Read in data -----------------------------------------------------------------

invisible(lapply(distances, function(x){
  fimo_path <- paste0(x, "/fimo_out/")
  
  fimo_file <- "fimo.tsv"
  
  fimo_res <- read.table(paste0(working_dir, gene_list, "/", meme_dir,
                                "/", fimo_path, fimo_file), sep = "\t",
                         header = TRUE)
  
  
  
  # Seperate sequence_name column to individual columns of information and fnd
  # region information of the motifs
  fimo_res <- fimo_res %>%
    dplyr::mutate(gene_id = sub("::.*", "", sequence_name),
                  chromosome = sub(".*::(.*?):.*", "\\1", sequence_name),
                  region_start = sub(".*:([0-9]*?)-.*", "\\1",
                                     sequence_name),
                  region_end = sub(".*-", "", sequence_name)) %>%
    dplyr::mutate(motif_start = as.numeric(region_start) + start,
                  motif_end = as.numeric(region_start) + stop)
  
  # Add gene ids
  mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
  genes <- unique(fimo_res$gene_id)
  G_list <- getBM(filters = "ensembl_gene_id",
                  attributes = c("ensembl_gene_id", "external_gene_name"),
                  values = genes, mart = mart)
  
  gene_names <- G_list$external_gene_name
  names(gene_names) <- G_list$ensembl_gene_id
  
  fimo_res$gene_name <- gene_names[fimo_res$gene_id]
  
  # Add no conservation score for these
  fimo_res_remove <- fimo_res %>%
    dplyr::filter(gene_id %in% skip_genes) %>%
    mutate(average_conservation = "NA",
           all_conservation = "NA")
  
  fimo_res <- fimo_res %>%
    dplyr::filter(!(gene_id %in% skip_genes))
  
  # Determine conservation scores and add them to the fimo output
  conservation_scores <- lapply(1:nrow(fimo_res), function(x){
    new_res <- fimo_res %>%
      dplyr::slice(x)
    return(get_conservation(new_res,
                            conservation_path))
  })
  
  conservation_scores <- do.call(rbind, conservation_scores)
  
  fimo_res <- cbind(fimo_res, conservation_scores)
  
  fimo_res <- rbind(fimo_res, fimo_res_remove)
  
  # Save data frame
  write.table(fimo_res, file = paste0(working_dir, gene_list, "/", meme_dir,
                                      "/", fimo_path, 
                                      "fimo_conservation.tsv"),
              sep = "\t", row.names = FALSE, col.names = TRUE)
}))
