library(Gviz)
library(biomaRt)
library(GenomicFeatures)
library(tidyverse)


bam_file_wt <- file.path("path", "to", "wt_bam.bam")

bam_file_cfko <- file.path("path", "to", "cfko_bam.bam")

gtf_file <- file.path("path", "to", "gtf")

gtf <- rtracklayer::import(gtf_file)

# I used this because I also used ensembl
gtf <- renameSeqlevels(gtf, value = paste0("chr", seqlevels(gtf)))

gtf_transcript <- gtf[gtf$type == "exon",]

wt_bam_track <- AlignmentsTrack(bam_file_wt)

cfko_bam_track <- AlignmentsTrack(bam_file_cfko)

# Get ensembl data too
# mart <- useEnsembl("genes")
# all_marts <- listDatasets(mart)
# all_marts[grepl("Ferret", all_marts$description),]
# 112 mpfuro_gene_ensembl Ferret genes (MusPutFur1.0) MusPutFur1.0
ferret <- useEnsembl("ensembl", dataset = "mpfuro_gene_ensembl") 

options(ucscChromosomeNames=FALSE)


plot_gviz <- function(gtf_all, gtf_transcript, gene, bam_path_wt, bam_path_cfko,
                      ensembl_mart, ...){
  
  # pull chromosome, start and end from gtf file
  gene_position <- gtf_all[gtf_all$type == "gene" & gtf_all$gene_name == gene,]
  chromosome <- as.character(seqnames(gene_position))
  range <- data.frame(gene_position@ranges)
  afrom <- min(range$start, range$end)
  ato <- max(range$start, range$end)
  
  # Make a track of all transcripts
  gtf_gene <- gtf_transcript[gtf_transcript$gene_name == gene, ]
  
  grtrack <- GeneRegionTrack(gtf_gene, name = paste0(gene, "_transcripts"))
  
  # Make a track for ensembl info ('m not sure that the "genome" is correct)
  bmt <- BiomartGeneRegionTrack(biomart = ensembl_mart, symbol = gene,
                                filter = list(with_refseq_mrna = TRUE),
                                stacking = "dense",
                                genome = "MusPutFur1.0", name = "ensembl")
  
  # rename chromosomes to be consistent
  bmt@range <- renameSeqlevels(bmt@range,
                               value = paste0("chr", seqlevels(bmt@range)))
  
  chromosome(bmt) <- chromosome
  
  
  # Make plot of all data
  plotTracks(c(bmt, grtrack, bam_path_wt, bam_path_cfko), 
             chromosome=chromosome, 
             from=afrom,to=ato, ...)
  
  
}

# Make plots
plot_gviz(gtf_all = gtf, gtf_transcript = gtf_transcript, gene = "Pax6",
          bam_path_wt = wt_bam_track, bam_path_cfko = cfko_bam_track,
          ensembl_mart = mouse, type = "coverage")

# type = "coverage" goes directly to `plotTracks` and will make it so the genes
# aren't plotted
