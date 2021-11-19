library(tidyverse)
library(cowplot)

working_dir <- 
  "/Users/wellskr/Documents/Analysis/Richard_Benniger/benniger_tf_regulators/results"

# Set theme
ggplot2::theme_set(ggplot2::theme_classic(base_size = 12))

# NFAT targets -----------------------------------------------------------------
gene_list <- "nfat_targets_nfat_tfs"

meme_dir <- "meme_NFAT_LIST"

conservation_dir <- "UCSC_phast"

conservation_path <- paste0(working_dir, gene_list, "/",
                            conservation_dir)

# Distances to plot
distance <- "5000_5000_test"

# genes to plot
plot_list <- c("Npas4", "Fos", "Nr4a1", "Nr4a2") # NFAT binding

# Bed file so we can plot the gene as well
bed_file <- read.table("files/Mus_musculus.GRCm38.96.bed",
                       header = FALSE, sep = "\t",
                       stringsAsFactors = FALSE, quote = "")

colnames(bed_file) <- c("chr", "start", "end", "gene_id",
                        "other", "strand", "source", "type",
                        "score", "gene_info")


# Load in fimo file
fimo_path <- file.path(distance, "fimo_out")
  
fimo_file <- "fimo_conservation.tsv"
  
fimo_res <- read.table(file.path(working_dir, gene_list, meme_dir,
                                 fimo_path, fimo_file), sep = "\t",
                       header = TRUE)
  
gene_plot <- lapply(plot_list, function(y){
  # Select only the gene and keep the highest score if there are duplicates
  if(y %in% fimo_res$gene_name){
    gene_df <- fimo_res %>%
      dplyr::filter(gene_name == y) %>%
      dplyr::mutate(log10p_val = -log10(p.value),
                    motif_pos = start - 5000) %>%
      dplyr::group_by(start) %>%
      dplyr::top_n(1, log10p_val)
    
    # Identify ens id to get strand information, negative = left arrow, positive = right
    ens_id <- gene_df$gene_id[1]
    
    bed_df <- bed_file %>%
      dplyr::filter(gene_id == ens_id)
    
    strand <- bed_df[, "strand"]
    
    if(strand == "+"){
      segment_start <- 0
      segment_end <- 1000
    } else if(strand == "-"){
      segment_start <- 0
      segment_end <- -1000
    }
    
    color_limits = c(10, 15)
    y_limits = c(3.9, 5.9)
    y_val = 3.85
    
    make_plot <- ggplot2::ggplot(data = gene_df,
                                 mapping = ggplot2::aes(x = motif_pos,
                                                        y = log10p_val,
                                                        color = score)) +
      ggplot2::geom_point(size = 3) +
      viridis::scale_color_viridis(option = "magma", limits = color_limits) +
      ggplot2::geom_segment(x = segment_start, xend = segment_end,
                            y = y_val, yend = y_val, color = "black",
                            arrow = arrow(length = unit(0.03, "npc"))) +
      ggplot2::xlim(c(-5000, 5000)) +
      ggplot2::ylim(y_limits) +
      ggplot2::labs(x = paste0("position relative to ", y, " TSS"),
                    y = "-log10 FIMO p-value",
                    color = "FIMO score",
                    title = paste0("NFAT motifs near ", y, " TSS")) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                         vjust = 0.5,
                                                         hjust=0.5))
    return(make_plot)
  }
})

legend <- get_legend(
  gene_plot[[1]] + theme(legend.box.margin = margin(0, 0, 0, 12))
)

# Remove the legend
gene_plot <- lapply(gene_plot, function(x){
  return(x + ggplot2::theme(legend.position = "None"))
})

save_plot <- cowplot::plot_grid(plotlist = gene_plot,
                                nrow = 2,
                                ncol = 2)

save_plot_legend <- plot_grid(save_plot, legend, rel_widths = c(3, .4))

pdf(file.path(working_dir, "R_analysis", "images",
              "nfat_binding_plots_final.pdf"),
    height = 5, width = 10)
print(save_plot_legend)
dev.off()

# Gjd2 targets -----------------------------------------------------------------
gene_list <- "ca_targets_ca_tfs_mouse"

meme_dir <- "meme_CA_LIST"

conservation_dir <- "UCSC_phast"

conservation_path <- paste0(working_dir, gene_list, "/",
                            conservation_dir)

distance <- "5000_5000_test"

Npas4_site <- 8400

plot_gene <- "Gjd2"

plot_list <- c("FOS", "FOS::JUND",
               "NR4A1", "NR4A2") #For GJD2 indicate Npas4, Fos and Nr4a1/2 binding? In ca motifs

bed_file <- read.table("files/Mus_musculus.GRCm38.96.bed",
                       header = FALSE, sep = "\t",
                       stringsAsFactors = FALSE, quote = "")

colnames(bed_file) <- c("chr", "start", "end", "gene_id",
                        "other", "strand", "source", "type",
                        "score", "gene_info")

colors <- RColorBrewer::brewer.pal(4, "Set1")
#colors <- c(colors[1:4], colors[7])
names(colors) <- plot_list

fimo_path <- file.path(distance, "fimo_out")

fimo_file <- "fimo_conservation.tsv"
  
fimo_res <- read.table(file.path(working_dir, gene_list, meme_dir,
                                 fimo_path, fimo_file), sep = "\t",
                       header = TRUE)
  
  # First subset to only "Gjd2"
gene_df <- fimo_res %>%
  dplyr::filter(gene_name == plot_gene &
                  motif_alt_id %in% plot_list) %>%
  dplyr::mutate(log10p_val = -log10(p.value),
                motif_pos = start - 5000) %>%
  dplyr::group_by(start) %>%
  dplyr::top_n(1, log10p_val)

ens_id <- gene_df$gene_id[1]
    
bed_df <- bed_file %>%
  dplyr::filter(gene_id == ens_id)
    
strand <- bed_df[, "strand"]
  
if(strand == "+"){
  segment_start <- 0
  segment_end <- 1000
} else if(strand == "-"){
  segment_start <- 0
  segment_end <- -1000
}
    

y_limits = c(3.9, 5.9)
y_val = 3.85

make_plot <- ggplot2::ggplot(data = gene_df,
                             mapping = ggplot2::aes(x = motif_pos,
                                                    y = log10p_val,
                                                    color = motif_alt_id)) +
  ggplot2::geom_point(size = 3) +
  ggplot2::scale_color_manual(values = colors) +
  ggplot2::geom_segment(x = segment_start, xend = segment_end,
                        y = y_val, yend = y_val, color = "black",
                        arrow = arrow(length = unit(0.03, "npc"))) +
  ggplot2::geom_segment(x = Npas4_site, xend = Npas4_site,
                        y = 4.3, yend = 3.9, color = "red",
                        arrow = arrow(length = unit(0.03, "npc"))) +
  ggplot2::geom_vline(xintercept = 5000, linetype = "dashed",
                      color = "black", size = 1) +
  ggplot2::xlim(c(-5000, 9000)) +
  ggplot2::ylim(y_limits) +
  ggplot2::labs(x = paste0("position relative to ", plot_gene, " TSS"),
                y = "-log10 FIMO p-value",
                color = "TF motif",
                title = paste0("TF motifs near ", plot_gene, " TSS")) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                     vjust = 0.5,
                                                     hjust=0.5))

pdf(file.path(working_dir, "R_analysis", "images",
              "gjd2_mouse_binding_plots_final.pdf"),
    height = 3, width = 7)
print(make_plot)
dev.off()

# Gjd2 targets human -----------------------------------------------------------
gene_list <- "ca_targets_ca_tfs_human"

meme_dir <- "meme_CA_LIST"

conservation_dir <- "UCSC_phast"

conservation_path <- paste0(working_dir, gene_list, "/",
                            conservation_dir)

distance <- "5000_5000_test"

plot_gene <- "GJD2"

plot_list <- c("FOS", "FOS::JUND",
               "NR4A1", "NR4A2") #For GJD2 indicate Npas4, Fos and Nr4a1/2 binding? In ca motifs

bed_file <- read.table("files/Homo_sapiens.GRCh38.95.bed",
                       header = FALSE, sep = "\t",
                       stringsAsFactors = FALSE, quote = "")

colnames(bed_file) <- c("chr", "start", "end", "gene_id",
                        "other", "strand", "source", "type",
                        "score", "gene_info")

colors <- RColorBrewer::brewer.pal(4, "Set1")
#colors <- c(colors[1:4], colors[7])
names(colors) <- plot_list

fimo_path <- file.path(distance, "fimo_out")

fimo_file <- "fimo_conservation.tsv"

fimo_res <- read.table(file.path(working_dir, gene_list, meme_dir,
                                 fimo_path, fimo_file), sep = "\t",
                       header = TRUE)

# First subset to only "Gjd2"
gene_df <- fimo_res %>%
  dplyr::filter(gene_name == plot_gene &
                  motif_alt_id %in% plot_list) %>%
  dplyr::mutate(log10p_val = -log10(p.value),
                motif_pos = start - 5000) %>%
  dplyr::group_by(start) %>%
  dplyr::top_n(1, log10p_val)

ens_id <- gene_df$gene_id[1]

bed_df <- bed_file %>%
  dplyr::filter(gene_id == ens_id)

strand <- bed_df[, "strand"]

if(strand == "+"){
  segment_start <- 0
  segment_end <- 1000
} else if(strand == "-"){
  segment_start <- 0
  segment_end <- -1000
}


y_limits = c(3.9, 5.9)
y_val = 3.85

make_plot <- ggplot2::ggplot(data = gene_df,
                             mapping = ggplot2::aes(x = motif_pos,
                                                    y = log10p_val,
                                                    color = motif_alt_id)) +
  ggplot2::geom_point(size = 3) +
  ggplot2::scale_color_manual(values = colors) +
  ggplot2::geom_segment(x = segment_start, xend = segment_end,
                        y = y_val, yend = y_val, color = "black",
                        arrow = arrow(length = unit(0.03, "npc"))) +
  ggplot2::xlim(c(-5000, 5000)) +
  ggplot2::ylim(y_limits) +
  ggplot2::labs(x = paste0("position relative to ", plot_gene, " TSS"),
                y = "-log10 FIMO p-value",
                color = "TF motif",
                title = paste0("TF motifs near ", plot_gene, " TSS")) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                     vjust = 0.5,
                                                     hjust=0.5))

pdf(file.path(working_dir, "R_analysis", "images",
              "gjd2_human_binding_plots_final.pdf"),
    height = 3.5, width = 7)
print(make_plot)
dev.off()

