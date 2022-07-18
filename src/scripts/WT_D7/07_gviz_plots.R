library(Gviz)
library(biomaRt)
library(GenomicFeatures)
library(tidyverse)
library(here)

server_dir <- "/Users/wellskr/Documents/sshfs"

sample <- "WT_D7"

analysis_dir <- file.path(server_dir,
                          "Analysis/Lori_Sussel/Ferret")

start_dir <- file.path(analysis_dir, "sc_220429")


bam_file <- file.path(start_dir, "results", sample,
                      "outs/possorted_genome_bam.bam")

# gtf_file <- file.path(server_dir, 
#                       "Analysis/ref/annotation/ferret",
#                       "ASM1176430v1.1/GCF_011764305.1_ASM1176430v1.1_genomic.gtf")

gtf_file <- here("files/GCF_011764305.1_ASM1176430v1.1_genomic.gtf")

save_dir <- file.path(here(), "results", sample, "R_analysis")

gtf <- rtracklayer::import(gtf_file)


#gtf <- renameSeqlevels(gtf, value = paste0("chr", seqlevels(gtf)))

gtf_transcript <- gtf[gtf$type == "exon",]

test <- AlignmentsTrack(bam_file)

options(ucscChromosomeNames=FALSE)


plot_gviz <- function(gtf_all, gtf_transcript, gene, bam_path,
                      upstream = 0, downstream = 0, ...){
  
  # pull chromosome, start and end from gtf file
  gene_position <- gtf_all[gtf_all$type == "gene" & gtf_all$gene_id == gene,]
  chromosome <- as.character(seqnames(gene_position))
  range <- data.frame(gene_position@ranges)
  afrom <- min(range$start, range$end) - upstream
  ato <- max(range$start, range$end) + downstream
  
  # Make a track of all transcripts
  gtf_gene <- gtf_transcript[gtf_transcript$gene_id == gene, ]
  
  grtrack <- GeneRegionTrack(gtf_gene, name = paste0(gene, "_transcripts"))
  
  
  # Make plot of all data
  plotTracks(c(grtrack, bam_path), 
             chromosome=chromosome, 
             from=afrom,to=ato, ...)
  
  
}

# Make plots
pdf(file.path(save_dir, "images/pax6_coverage.pdf"))

plot_gviz(gtf_all = gtf, gtf_transcript = gtf_transcript, gene = "PAX6",
          bam_path = test, upstream = 1000,
          downstream = 1000)

dev.off()

pdf(file.path(save_dir, "images/pdx1_coverage.pdf"))

plot_gviz(gtf_all = gtf, gtf_transcript = gtf_transcript, gene = "PDX1",
          bam_path = test, upstream = 20000,
          downstream = 20000)

dev.off()


pdf(file.path(save_dir, "images/ins_coverage.pdf"))

plot_gviz(gtf_all = gtf, gtf_transcript = gtf_transcript, gene = "INS",
          bam_path = test, upstream = 1000, downstream = 1000)

dev.off()

# Make plots
png(file.path(save_dir, "images/pax6_coverage.png"))

plot_gviz(gtf_all = gtf, gtf_transcript = gtf_transcript, gene = "PAX6",
          bam_path = test, upstream = 1000,
          downstream = 1000)

dev.off()

png(file.path(save_dir, "images/pdx1_coverage.png"))

plot_gviz(gtf_all = gtf, gtf_transcript = gtf_transcript, gene = "PDX1",
          bam_path = test, upstream = 20000,
          downstream = 20000)

dev.off()


png(file.path(save_dir, "images/ins_coverage.png"))

plot_gviz(gtf_all = gtf, gtf_transcript = gtf_transcript, gene = "INS",
          bam_path = test, upstream = 1000, downstream = 1000)

dev.off()


png(file.path(save_dir, "images/nkx6.1_coverage.png"))

plot_gviz(gtf_all = gtf, gtf_transcript = gtf_transcript, gene = "NKX6-1",
          bam_path = test, upstream = 5000, downstream = 5000)

dev.off()


png(file.path(save_dir, "images/sox9_coverage.png"))

plot_gviz(gtf_all = gtf, gtf_transcript = gtf_transcript, gene = "SOX9",
          bam_path = test, upstream = 5000, downstream = 5000)

dev.off()

png(file.path(save_dir, "images/neurog3_coverage.png"))

plot_gviz(gtf_all = gtf, gtf_transcript = gtf_transcript, gene = "NEUROG3",
          bam_path = test, upstream = 5000, downstream = 5000)

dev.off()

