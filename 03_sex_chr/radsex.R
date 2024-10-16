# R codes to plot radsex output

#library(devtools)
#install_github("biodray/sgtr")

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(sgtr)
library(ggplot2)
library(dplyr)

# tile plot 
#png("60indv/distribution.png", units = "in", width = 8, height = 5, res=400)
radsex_distrib("60indv/distribution.tsv", 
               group_labels = c("Males", "Females"),
               width=13,height=13)
#dev.off()


# Marker depths plot with specified groups, labels, and group info to color
# individual labels
#png("60indv/signif_marker.png", units = "in", width = 8, height = 5, res=400)
radsex_marker_depths("60indv/signif_markers.tsv",
                     group_info_file = "60indv/popmap.txt",
                     group_labels = c("Females", "Males"),
                     label_colors = c("firebrick2", "dodgerblue3"))
#dev.off()


# Manhattan plot
radsex_map_manhattan("60indv/alignment_results.tsv",chromosomes_file="60indv/chr.txt")
# need to modify the y label 
radsex_map_manhattan_ylabel <- function (input_file, chromosomes_file = NA, detect_chromosomes = TRUE, 
          unplaced_label = "U.", output_file = NA, width = 12, height = 6, 
          res = 300, colors = c("dodgerblue3", "darkgoldenrod2"), bg_colors = c("grey85", 
                                                                                "white"), point_size = 0.5, x_title = NA, show_chromosomes_names = TRUE, 
          chromosomes_as_numbers = FALSE, show_signif_line = TRUE, 
          signif_line_color = "red") 
{
  chromosomes <- load_chromosome_names(chromosomes_file)
  data <- load_genome_metrics(input_file, chromosomes = chromosomes, 
                              detect_chromosomes = detect_chromosomes, unplaced_label = unplaced_label, 
                              comment_char = "#", comment_sep = ";", comment_internal_sep = ":")
  data$data$P <- -log(data$data$P, 10)
  y_label <- expression(bold(paste("-log"[10], "(p"["Association with sex"],")")))
  tracks <- list(single_metric_track("P", colors = colors, 
                                     point_size = point_size, alpha = 1, type = "points", 
                                     label = y_label, bg_colors = bg_colors, ylim = NA, major_lines_x = FALSE, 
                                     major_lines_y = TRUE, legend_position = "none"))
  if (show_signif_line) {
    s <- as.numeric(data$properties$signif_threshold)
    n_markers <- as.numeric(data$properties$n_markers)
    signif_threshold <- -log(s/n_markers, 10)
    tracks[[1]]$h_lines <- list(h_line(signif_threshold, 
                                       label = paste0("p<", s), color = signif_line_color, 
                                       type = 2, size = 0.75, label_font_size = 5))
  }
  m <- draw_manhattan_plot(data$data, data$lengths, tracks, 
                           output_file = output_file, width = width, track_height = height, 
                           res = res, x_title = x_title, show_chromosomes_names = show_chromosomes_names, 
                           chromosomes_as_numbers = chromosomes_as_numbers)
  return(m)
}
#png("60indv/manhattan.png", units = "in", width = 8, height = 5, res=400)
radsex_map_manhattan_ylabel("60indv/alignment_results.tsv",chromosomes_file="60indv/chr.txt")
#dev.off()



