library(tidyverse)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

base_path <- "/Users/wellskr/Documents/Analysis/Richard_Benniger/benniger_tf_regulators"

read_dir <- file.path(base_path, "results/nfat_g_gd_targets_all_tfs/homer_NFAT_LIST/300_50_test")

save_dir <- file.path(base_path, "results/R_analysis/files")

read_file <- file.path(read_dir, "knownResults.txt")

results <- read.table(read_file, sep = "\t", skip = 1)
colnames(results) <- c("Motif_name", "consensus", "p_value", "log_p", "q_value",
                       "target_count", "target_percent", "background_count",
                       "background_percent")

tf_list <- c("NFAT", "CREB", "Fos", "Mef2",
             "Npas4", "NFkB-p65-Rel", "NeuroD")

tf_results <- lapply(tf_list, function(x){
  result_return <- results[grepl(x, results$Motif_name), ]
})

tf_results <- do.call(rbind, tf_results)

tf_results <- tf_results[order(tf_results$p_value), ]

write.csv(tf_results, file = file.path(save_dir, "homer_pvals.csv"))
