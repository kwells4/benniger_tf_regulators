library(tidyverse)

working_dir <- 
  "/Users/wellskr/Documents/Analysis/Richard_Benniger/benniger_tf_regulators/results"

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 12))


gene_list <- "ca_targets_ca_tfs_human"

meme_dir <- "meme_CA_LIST"

conservation_dir <- "UCSC_phast"

conservation_path <- paste0(working_dir, gene_list, "/",
                            conservation_dir)

distances <- c("300_50_test", "2000_2000_test", "2000_0_test", "5000_5000_test")

plot_gene <- "GJD2"

plot_list <- c("Npas4", "FOS", "FOS::JUN", "FOS::JUND",
               "NR4A1", "NR4A2") #For GJD2 indicate Npas4, Fos and Nr4a1/2 binding? In ca motifs

bed_file <- read.table("files/Homo_sapiens.GRCh38.95.bed",
                       header = FALSE, sep = "\t",
                       stringsAsFactors = FALSE, quote = "")

colnames(bed_file) <- c("chr", "start", "end", "gene_id",
                        "other", "strand", "source", "type",
                        "score", "gene_info")

colors <- RColorBrewer::brewer.pal(7, "Set1")
colors <- c(colors[1:5], colors[7])
names(colors) <- plot_list

# Color by score
all_plots <- lapply(distances, function(x){
  print(x)
  fimo_path <- file.path(x, "fimo_out")
  
  fimo_file <- "fimo_conservation.tsv"
  
  fimo_res <- read.table(file.path(working_dir, gene_list, meme_dir,
                                   fimo_path, fimo_file), sep = "\t",
                         header = TRUE)
  
  # First subset to only "Gjd2"
  gene_df <- fimo_res %>%
    dplyr::filter(gene_name == plot_gene &
                    motif_alt_id %in% plot_list) %>%
    dplyr::mutate(log10p_val = -log10(p.value)) %>%
    dplyr::group_by(start) %>%
    dplyr::top_n(1, log10p_val)
  
  if(nrow(gene_df > 0)){
    
    region_start <- gene_df$region_start[1]
    region_end <- gene_df$region_end[1]
    chromosome <- gene_df$chromosome[1]
    ens_id <- gene_df$gene_id[1]
    
    bed_df <- bed_file %>%
      dplyr::filter(gene_id == ens_id)
    
    gene_start <- bed_df[ , "start"]
    gene_end <- bed_df[, "end"]
    
    #y_limits = c(min(gene_df$log10p_val) - 0.1, max(gene_df$log10p_val) + 0.01)
    #y_val = min(gene_df$log10p_val) - 0.05
    y_limits = c(3.9, 5.9)
    y_val = 3.85
    
    make_plot <- ggplot2::ggplot(data = gene_df,
                                 mapping = ggplot2::aes(x = motif_start,
                                                        y = log10p_val,
                                                        color = motif_alt_id)) +
      ggplot2::geom_point(size = 3) +
      ggplot2::scale_color_manual(values = colors) +
      ggplot2::geom_segment(x = gene_start, xend = gene_end,
                            y = y_val, yend = y_val, color = "black") +
      ggplot2::xlim(c(region_start, region_end)) +
      ggplot2::ylim(y_limits) +
      ggplot2::labs(x = paste0("position on mm10 chromosome ", chromosome),
                    y = "-log10 FIMO p-value",
                    color = "TF motif",
                    title = paste0("TF binding around ", plot_gene)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                         vjust = 0.5,
                                                         hjust=0.5))
    
    pdf(file.path(working_dir, "R_analysis", "images",
                  paste0(x, "_", plot_gene, "_binding_plots_human.pdf")),
        height = 4, width = 5)
    print(make_plot)
    dev.off()
    return(make_plot)
  }
  
})
