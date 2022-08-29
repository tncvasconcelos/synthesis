# rm(list=ls())

source("00_utility_functions_synthesis.R")

# Reading trees from Smith and Brown (2018)
# These trees have been previously taxized for the GBIF taxonomic backbone
whole_taxized_tree <- read.tree("data/ALLMB.taxized.noauthor.tre")
gb_taxized_tree <- read.tree("data/GBMB.taxized.noauthor.tre")

# tip name cleaning:
# remove indet, unmatched, etc
whole_taxized_tree$tip.label <- gsub("_"," ", whole_taxized_tree$tip.label) 
gb_taxized_tree$tip.label <- gsub("_"," ", gb_taxized_tree$tip.label) 

cleaned_tree <- drop.tip(whole_taxized_tree, grep("UNMATCHED",whole_taxized_tree$tip.label))
cleaned_tree <- drop.tip(cleaned_tree, which(duplicated(cleaned_tree$tip.label)))
cleaned_tree <- drop.tip(cleaned_tree, grep("tip_to_drop", simplify.names.taxize(cleaned_tree$tip.label)))

cleaned_tree_gb <- drop.tip(gb_taxized_tree, grep("UNMATCHED",gb_taxized_tree$tip.label))
cleaned_tree_gb <- drop.tip(cleaned_tree_gb, which(duplicated(cleaned_tree_gb$tip.label)))
cleaned_tree_gb <- drop.tip(cleaned_tree_gb, grep("tip_to_drop", simplify.names.taxize(cleaned_tree_gb$tip.label)))

####
# making a "trait" dataset for species with no genetic data
gb_gap_table <- data.frame(species=cleaned_tree$tip.label, gb_data=NA)
gb_gap_table$gb_data[which(cleaned_tree$tip.label %in% cleaned_tree_gb$tip.label)] <- "has_data"
gb_gap_table$gb_data[which(!cleaned_tree$tip.label %in% cleaned_tree_gb$tip.label)] <- "no_data"

write.csv(gb_gap_table, file="gb_gap_table.csv", row.names=F)

####
# making a "trait" dataset for species with no seed mass data (data from BIEN)
list_of_one_trait <- BIEN_trait_trait(trait="seed mass")
list_of_one_trait0 <- subset(list_of_one_trait, !is.na(list_of_one_trait$scrubbed_species_binomial))
list_of_one_trait0 <- subset(list_of_one_trait0, list_of_one_trait0$access=="public")
seed_mass <- unique(list_of_one_trait0$scrubbed_species_binomial)
seed_mass_taxized <- gbif.taxize(seed_mass)

seed_mass_taxized <- unname(seed_mass_taxized)
seed_mass_taxized <- fix.names.taxize(seed_mass_taxized)
seed_mass_taxized <- subset(seed_mass_taxized, !grepl("UNMATCHED",seed_mass_taxized))
seed_mass_taxized <- simplify.names.taxize(seed_mass_taxized)

seed_gap_table <- data.frame(species=cleaned_tree$tip.label, seed_mass=NA)
seed_gap_table$seed_mass[which(seed_gap_table$species %in% seed_mass_taxized)] <- "has_data"
seed_gap_table$seed_mass[which(!seed_gap_table$species %in% seed_mass_taxized)] <- "no_data"

write.csv(seed_gap_table, file="seed_mass_gap_table.csv", row.names=F)

######
# making a "trait" dataset for species with no ploidy data (data from Rice et al., 2019)
ploidy <- read.csv("all_ploidy_single_df_taxized.csv")
ploidy_taxized <- fix.names.taxize(ploidy$Taxon)
ploidy_taxized <- subset(ploidy_taxized, !grepl("UNMATCHED",ploidy_taxized))
ploidy_taxized <- simplify.names.taxize(ploidy_taxized)

ploidy_gap_table <- data.frame(species=cleaned_tree$tip.label, ploidy=NA)
ploidy_gap_table$ploidy[which(ploidy_gap_table$species %in% ploidy_taxized)] <- "has_data"
ploidy_gap_table$ploidy[which(!ploidy_gap_table$species %in% ploidy_taxized)] <- "no_data"

write.csv(ploidy_gap_table, file="ploidy_gap_table.csv", row.names=F)
