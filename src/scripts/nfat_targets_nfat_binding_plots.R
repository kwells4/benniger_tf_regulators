library(tidyverse)

working_dir <- 
  "/Users/wellskr/Documents/Analysis/Richard_Benniger/benniger_tf_regulators/results"

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 12))


gene_list <- "nfat_targets_nfat_tfs"

meme_dir <- "meme_NFAT_LIST"

conservation_dir <- "UCSC_phast"

conservation_path <- paste0(working_dir, gene_list, "/",
                            conservation_dir)

# Distances to plot
distances <- c("300_50_test", "2000_2000_test", "2000_0_test", "5000_5000_test")

# genes to plot
plot_list <- c("Npas4", "Fos", "Nr4a1", "Nr4a2") # NFAT binding

# Bed file so we can plot the gene as well
bed_file <- read.table("files/Mus_musculus.GRCm38.96.bed",
                       header = FALSE, sep = "\t",
                       stringsAsFactors = FALSE, quote = "")

colnames(bed_file) <- c("chr", "start", "end", "gene_id",
                        "other", "strand", "source", "type",
                        "score", "gene_info")

# Color by score
all_plots <- lapply(distances, function(x){
  print(x)
  
  # Load in fimo file
  fimo_path <- file.path(x, "fimo_out")
  
  fimo_file <- "fimo_conservation.tsv"
  
  fimo_res <- read.table(file.path(working_dir, gene_list, meme_dir,
                                   fimo_path, fimo_file), sep = "\t",
                         header = TRUE)
  
  gene_plot <- lapply(plot_list, function(y){
    # Select only the gene and keep the highest score if there are duplicates
    if(y %in% fimo_res$gene_name){
      gene_df <- fimo_res %>%
        dplyr::filter(gene_name == y) %>%
        dplyr::mutate(log10p_val = -log10(p.value)) %>%
        dplyr::group_by(start) %>%
        dplyr::top_n(1, log10p_val)
    
      # Identify region and chromosome information for limits
      region_start <- gene_df$region_start[1]
      region_end <- gene_df$region_end[1]
      chromosome <- gene_df$chromosome[1]
      ens_id <- gene_df$gene_id[1]
    
      bed_df <- bed_file %>%
        dplyr::filter(gene_id == ens_id)
    
      gene_start <- bed_df[ , "start"]
      gene_end <- bed_df[, "end"]
    
      color_limits = c(min(gene_df$score) - 1, max(gene_df$score) + 2)
      #y_limits = c(min(gene_df$log10p_val) - 0.1, max(gene_df$log10p_val) + 0.01)
      #y_val = min(gene_df$log10p_val) - 0.05
      y_limits = c(3.9, 5.9)
      y_val = 3.85
    
      make_plot <- ggplot2::ggplot(data = gene_df,
                                   mapping = ggplot2::aes(x = motif_start,
                                                          y = log10p_val,
                                                          color = score)) +
        ggplot2::geom_point(size = 3) +
        viridis::scale_color_viridis(option = "magma", limits = color_limits) +
        ggplot2::geom_segment(x = gene_start, xend = gene_end,
                              y = y_val, yend = y_val, color = "black") +
        ggplot2::xlim(c(region_start, region_end)) +
        ggplot2::ylim(y_limits) +
        ggplot2::labs(x = paste0("position on mm10 chromosome ", chromosome),
                      y = "-log10 FIMO p-value",
                      color = "FIMO score",
                      title = paste0("NFAT binding around ", y)) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                           vjust = 0.5,
                                                           hjust=0.5))
      return(make_plot)
    } else {
      return(NA)
    }
  })
  # Only make a plot if any values were found
  if(FALSE %in% names(table(is.na(gene_plot)))){
    nplots <- table(is.na(gene_plot))["FALSE"]
    save_plot <- cowplot::plot_grid(plotlist = gene_plot,
                                    nrow = ceiling(nplots/2),
                                    ncol = 2)
    pdf(file.path(working_dir, "R_analysis", "images",
                  paste0(x, "_nfat_binding_plots.pdf")),
        height = 8, width = 10)
    print(save_plot)
    dev.off()
    return(save_plot)
  }
})
