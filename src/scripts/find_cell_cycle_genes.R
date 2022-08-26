library(Seurat)
library(here)
library(tidyverse)

s_genes <- cc.genes$s.genes
g2m_genes <- cc.genes$g2m.genes

mapping_file <- read.csv(here("files/species_mapping_file.csv"))

s_genes_mapping <- mapping_file %>%
  dplyr::filter(Human.gene.name %in% s_genes) %>%
  dplyr::select(Human.gene.name, gene_id) %>%
  dplyr::distinct()


g2m_genes_mapping <- mapping_file %>%
  dplyr::filter(Human.gene.name %in% g2m_genes) %>%
  dplyr::select(Human.gene.name, gene_id) %>%
  dplyr::distinct()

saveRDS(list("s_genes" = s_genes_mapping,
             "g2m_genes" = g2m_genes_mapping),
        file = here("files/cell_cycle_genes.rds"))
