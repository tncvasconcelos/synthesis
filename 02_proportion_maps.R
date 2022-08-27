# rm(list=ls())

source("00_utility_functions_synthesis.R")

setwd("/Users/thaisvasconcelos/Desktop/manuscripts/evol_potential/synthesis/synthesis")

# Making maps for Figure 2
#########################
# If local
#setwd("~/Desktop/WCVP_special_issue/James_perennial_annual/life_history_houwie")
source("/Users/thaisvasconcelos/Desktop/WCVP_special_issue/WCVPtools/WCVPtools_functions.R")
dist_sample <- read.table("/Users/thaisvasconcelos/Desktop/WCVP_special_issue/wcvp_names_and_distribution_special_edition_2022/wcvp_distribution.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
names_sample <- read.table("/Users/thaisvasconcelos/Desktop/WCVP_special_issue/wcvp_names_and_distribution_special_edition_2022/wcvp_names.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")

#########################
all_vars <- merge(dist_sample, names_sample, by="plant_name_id")

#########################
# reference table for taxized names
#-----------------------------
# If local
reference_table <- list.files("/Users/thaisvasconcelos/Desktop/WCVP_special_issue/WCVPtools/taxized_reference_tables", full.names = T)
reference_table <- do.call(rbind, lapply(reference_table, read.csv))

#-----------------------------

###########
# Looking at the WCVP table and TDWG to clean GBIF points
#-----------------------------
# If local
path="/Users/thaisvasconcelos/Desktop/WCVP_special_issue/WCVPtools/wgsrpd-master/level3/level3.shp"
#-----------------------------

twgd_data <- suppressWarnings(maptools::readShapeSpatial(path))

gbif_data <- read.csv("woodness_gap_table.csv")
colnames(gbif_data) <- c("species","gbif")

organized_table_for_plot_total <- organize.bubble.plot(gbif_data, reference_table, all_vars, twgd_data)
organized_table_for_plot_no_data <- organize.bubble.plot(subset(gbif_data, gbif_data$gbif=="no_data"), reference_table, all_vars, twgd_data)

nas <- rep(NA, nrow(organized_table_for_plot_no_data))
proportion_table <- data.frame(sp_rich_prop=nas, sp_rich_total=nas, one_area=nas, lon=nas, lat=nas)
for(i in 1:nrow(organized_table_for_plot_no_data)) {
  tmp_sp_rich <- organized_table_for_plot_no_data$sp_rich[i]
  one_area <- organized_table_for_plot_no_data$one_area[i]
  total_sp_rich <- organized_table_for_plot_total$sp_rich[organized_table_for_plot_total$one_area == one_area]
  one_proportion <- round(tmp_sp_rich / total_sp_rich, 3)
  proportion_table$sp_rich_prop[i] <- one_proportion
  proportion_table$sp_rich_total[i] <- total_sp_rich
  proportion_table$one_area[i] <- one_area
  proportion_table$lon[i] <- organized_table_for_plot_no_data$lon[i]
  proportion_table$lat[i] <- organized_table_for_plot_no_data$lat[i]
}

hist(proportion_table$sp_rich_prop)
# cutting the data

#proportion_table$sp_rich_prop[which(proportion_table$sp_rich_prop > 0.5)] <- 0.4
#proportion_table$sp_rich_prop[which(between(proportion_table$sp_rich_prop, 0.2, 0.3))] <- 2
#proportion_table$sp_rich_prop[which(between(proportion_table$sp_rich_prop, 0.1, 0.2))] <- 1
#proportion_table$sp_rich_prop[which(between(proportion_table$sp_rich_prop, 0, 0.1))] <- 0

twgd_data01 <- sf::st_as_sf(twgd_data)
twgd_data01 <- merge(twgd_data01, proportion_table, by.x="LEVEL3_COD", by.y="one_area")

twgd_data01 <- subset(twgd_data01 , twgd_data01$LEVEL3_NAM!="Antarctica")

tmp_map1 <- ggplot(data = twgd_data01) +
  geom_sf(aes(fill = sp_rich_prop)) +
  scale_fill_viridis_c(option = "plasma", trans="sqrt") +
  theme_classic()
