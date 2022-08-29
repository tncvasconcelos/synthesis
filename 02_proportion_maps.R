# rm(list=ls())

source("00_utility_functions_synthesis.R")

# Making maps for Figure 2

#######################################################################################
#######################################################################################
##### IMPORTANT ##### 
# The following data is under an embargo and will probably be publicly available only in 2023.
# The use of the POWO data for this synthesis was kindly made available by Rafael Govaerts. 
# Please contact R.Govaerts@kew.org to have access to this dataset before 2023.
# Or access https://powo.science.kew.org for more information.
source("/Users/thaisvasconcelos/Desktop/WCVP_special_issue/WCVPtools/WCVPtools_functions.R")
dist_sample <- read.table("/Users/thaisvasconcelos/Desktop/WCVP_special_issue/wcvp_names_and_distribution_special_edition_2022/wcvp_distribution.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
names_sample <- read.table("/Users/thaisvasconcelos/Desktop/WCVP_special_issue/wcvp_names_and_distribution_special_edition_2022/wcvp_names.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
#######################################################################################
#######################################################################################

all_vars <- merge(dist_sample, names_sample, by="plant_name_id")

#########################
# reference table for taxized names
#-----------------------------
reference_table <- list.files("data/taxized_reference_tables", full.names = T)
reference_table <- do.call(rbind, lapply(reference_table, read.csv))

#-----------------------------

###########
# Loading TDWG shape files (see README in data/wgsrpd-master/ for more information)
#-----------------------------
path="data/wgsrpd-master/level3/level3.shp"
#-----------------------------

twgd_data <- suppressWarnings(maptools::readShapeSpatial(path))
twgd_data01 <- sf::st_as_sf(twgd_data)

# Figure 2A: molecular data for phylogenetic reconstruction
trait_table_A <- read.csv("gb_gap_table.csv")
proportion_table_A <- map.for.synthesis(trait_table=trait_table_A, reference_table, all_vars, twgd_data)
twgd_data_A <- merge(twgd_data01, proportion_table_A, by.x="LEVEL3_COD", by.y="one_area")
twgd_data_A <- subset(twgd_data_A , twgd_data_A$LEVEL3_NAM!="Antarctica")
figure2_map_A <- ggplot(data = twgd_data_A) +
  geom_sf(aes(fill = sp_rich_prop)) +
  scale_fill_viridis_c(option = "plasma") +
  theme_classic()

# Figure 2B: ploidy data
trait_table_B <- read.csv("ploidy_gap_table.csv")
proportion_table_B <- map.for.synthesis(trait_table=trait_table_B, reference_table, all_vars, twgd_data)
twgd_data_B <- merge(twgd_data01, proportion_table_B, by.x="LEVEL3_COD", by.y="one_area")
twgd_data_B <- subset(twgd_data_B , twgd_data_B$LEVEL3_NAM!="Antarctica")
figure2_map_B <- ggplot(data = twgd_data_B) +
  geom_sf(aes(fill = sp_rich_prop)) +
  scale_fill_viridis_c(option = "plasma") +
  theme_classic()

# Figure 2C: seed mass
trait_table_C <- read.csv("seedmass_gap_table.csv")
proportion_table_C <- map.for.synthesis(trait_table=trait_table_C, reference_table, all_vars, twgd_data)
twgd_data_C <- merge(twgd_data01, proportion_table_C, by.x="LEVEL3_COD", by.y="one_area")
twgd_data_C <- subset(twgd_data_C , twgd_data_C$LEVEL3_NAM!="Antarctica")
figure2_map_C <- ggplot(data = twgd_data_C) +
  geom_sf(aes(fill = sp_rich_prop)) +
  scale_fill_viridis_c(option = "plasma") +
  theme_classic()

#pdf("figure2_raw.pdf",height=15,width=10)
#grid.arrange(figure2_map_A, figure2_map_B, figure2_map_C, ncol=1, nrow = 3)
#dev.off()

