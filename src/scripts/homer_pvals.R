library(tidyverse)

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 10))

base_path <- "/Users/wellskr/Documents/Analysis/Richard_Benniger/benniger_tf_regulators/"

save_dir <- paste0(base_path, "results/")

image_dir <- paste0(save_dir, "R_analysis/images/")

homer_dirs <- c("ca_targets_all_tfs_human",
                "ca_targets_all_tfs_mouse",
                "nfat_g_gd_targets_all_tfs",
                "nfat_gd_gfk_targets_all_tfs",
                "nfat_targets_all_tfs")


read_in_files <- lapply(homer_dirs, function(x){
  full_dir <- dir(paste0(save_dir, x), pattern = "homer_.*",
                  full.names = TRUE)
  return(read.csv(file.path(full_dir, "300_50_test", "knownResults.txt"),
                  sep = "\t"))
})

names(read_in_files) <- homer_dirs

barplots <- lapply(names(read_in_files), function(x){
  plot_data <- read_in_files[[x]]
  plot_data$neg_logP <- -(plot_data$Log.P.value)
  plot_data <- plot_data %>%
    dplyr::filter(P.value < 0.05) %>%
    dplyr::arrange(neg_logP)
  plot_data$Motif.Name <- factor(plot_data$Motif.Name,
                                 levels = unique(plot_data$Motif.Name))
  plot_1 <- ggplot2::ggplot(plot_data, ggplot2::aes(y = Motif.Name,
                                                    x = neg_logP)) +
    ggplot2::geom_bar(stat = "identity")
  plot_2 <- plot_1 +
    ggplot2::theme(axis.title.y=element_blank(),
                    axis.text.y=element_blank(),
                    axis.ticks.y=element_blank())
  ggplot2::ggsave(filename = paste0(image_dir, x, "_names.pdf"),
                  plot = plot_1)
  ggplot2::ggsave(filename = paste0(image_dir, x, "_no_lab.pdf"),
                  plot = plot_2)
  return(list(plot_1, plot_2))
})

