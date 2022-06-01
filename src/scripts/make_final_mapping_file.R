library(here)
library(openxlsx)
library(tidyverse)

orthologs_file <- read.csv(here("files/orthologs/orthologs_all.csv"))

species_file <- read.csv("files/Complete_Multiple_species_comparison_feret_ref_genome.csv")


orthologs_short <- orthologs_file %>%
  dplyr::filter(ferret_orth_id %in% all_of(species_file$ENS_prot_id))


species_short <- species_file %>%
  dplyr::filter(ENS_prot_id %in% all_of(orthologs_file$ferret_orth_id))

# Add gene id to species file
species_short <- merge(species_short, orthologs_short, by.x = "ENS_prot_id",
                      by.y = "ferret_orth_id", all.x = TRUE, all.y = TRUE)

write.csv(species_short, "files/species_mapping_file.csv")
