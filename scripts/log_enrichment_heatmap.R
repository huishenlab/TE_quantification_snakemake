suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(optparse))
suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library(ggh4x))

plot_log_enrich <- function(log_enrich_file, output_log_enrich_pdf){
    log_enrichment <- fread(log_enrich_file, header = TRUE) %>% as.data.frame()
    
    log_enrichment$samples <- rep(unique(log_enrichment$samples), nrow(log_enrichment))

    log_enrichment$repClass.x <- factor(log_enrichment$repClass.x, 
                                                 levels = str_sort(unique(log_enrichment$repClass.x)))
    log_enrichment$samples <- factor(log_enrichment$samples, 
                                      levels = unique(log_enrichment$samples))
                                      
    log_enrichment_plt <- ggplot(log_enrichment, aes(repFamily, cells)) + 
                            geom_tile(aes(fill = log_enrichment)) +
                            facet_grid(samples ~ repClass.x, switch = "y",
                            scales = "free", space = "free_x") +
                            force_panelsizes(rows = c(0.3, 1,1,1,1)) +
                            scale_fill_gradientn(name = "log2(Enrichment)", 
                                                colours = c("dodgerblue4","whitesmoke","firebrick3"), 
                                                values = rescale(c(-3,0,3)),
                                                guide = "colorbar", limits=c(-3,3)) +
                            ggtitle(paste0(unique(log_enrichment$samples), " TE log2(enrichment)")) +
                            theme(plot.title = element_text(hjust = 0.5, color="black", size=14, face="bold"),
                                  strip.text.y.left = element_text(angle = 0),
                                  strip.background = element_blank(),
                                  strip.text.x = element_text(size=12, face="bold"),
                                  panel.border = element_rect(color = "black", fill = NA, linewidth = 1), 
                                  panel.grid.major = element_blank(), 
                                  panel.grid.minor = element_blank(), 
                                  panel.background = element_blank(), 
                                  axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
                                  axis.text.y=element_blank(),
                                  axis.ticks.y = element_blank()) +
                            xlab("Repeat Family") +
                            ylab("Condition")

     pdf(file=output_log_enrich_pdf,  width=12, height=3)
     print(log_enrichment_plt)
     dev.off()
}

log_enrich_file <- snakemake@input[['log_enrich_file']]
output_log_enrich_pdf <- snakemake@output[['output_log_enrich_pdf']]

## Generate log enrichment heatmap

message("\n\n** Generating log enrichment heatmap **")

plot_log_enrich(log_enrich_file, output_log_enrich_pdf)

message("UPDATE: Saved log enrichment heatmap as a pdf\n\n")