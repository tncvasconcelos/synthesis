# Load required packages:
library(ape)
library(taxize)
library(data.table)
library(BIEN)
library(maptools)
library(raster)
library(sf)
library(ggplot2)
library(maps)
library(ggthemes)
library(viridis)
library(ggridges)
library(gridExtra)

# Function to taxize species list according to GBIF taxonomy
gbif.taxize <- function (species) {
  sources <- taxize::gnr_datasources()
  gbif <- sources$id[sources$title == 'GBIF Backbone Taxonomy']
  gnr_resolve_x <- function(x, data_source_ids=gbif, best_match_only = TRUE) {
    new.name <- suppressWarnings(taxize::gnr_resolve(names=x, data_source_ids=data_source_ids, best_match_only=best_match_only)$matched_name)
    if(is.null(new.name)) {
      new.name <- paste0("UNMATCHED_",x)
    }
    return(new.name)
  }
  new.names <- pbapply::pbsapply(species, gnr_resolve_x, cl=4)
  return(new.names)
}

# Fixing taxize names 
fix.names.taxize <- function(focal_species_trees) {
  for(name_index in 1:length(focal_species_trees)){
    one_tmp_string <- focal_species_trees[name_index]
    if(any(grepl("[()]", one_tmp_string))){
      splitted_names <- strsplit(one_tmp_string," ")[[1]]
      begin_author <- which(grepl("[()]", splitted_names))[1]
      species_name <- paste0(splitted_names[1:(begin_author-1)], collapse=" ")
      author <- splitted_names[begin_author:length(splitted_names)]
      old_authors <- author[grep("[()]", author)]
      end_first_half <- floor(length(old_authors)/2)
      before <- old_authors[1:end_first_half]
      after <- old_authors[(end_first_half+1):(length(old_authors))]
      if(paste(before,collapse = " ") == paste(after, collapse=" ")) {
        author <- paste(author[1:(length(author)/2)], collapse=" ")
        focal_species_trees[name_index] <- paste0(species_name, " ", author, collapse=" ")
      } else {
        author <- paste(author, collapse=" ")
        focal_species_trees[name_index] <- paste0(species_name, " ", author, collapse=" ")
      }
    }
  }
  return(focal_species_trees)
}

# A way to simplify names in table and trees so that species names match again
simplify.names.taxize <- function(names) {
  results <- c()
  for(name_index in 1:length(names)){
    one_tmp_string <- names[name_index]
    splitted_names <- strsplit(one_tmp_string," ")[[1]]
    genus <- splitted_names[1]
    epiphet <- splitted_names[2]
    if(is.na(epiphet)) {
      full_name <- "tip_to_drop" # indet species
    } else {
    if(any(grepl("indet_sp",splitted_names))) {
      full_name <- "tip_to_drop" # indet species
    } else {
      if(stringr::str_detect(epiphet,"[[:upper:]]")) {
        full_name <- "tip_to_drop" # indet species
      } else {
        if(length(splitted_names) == 2) { 
          full_name <- paste(c(genus, epiphet), collapse = " ")
        } else {
          if(length(splitted_names) > 2) {
            complement <- splitted_names[3:length(splitted_names)]
            if(grepl("[()]", complement[1])) {
              full_name <- paste(c(genus, epiphet), collapse = " ")
            } else {
              if(stringr::str_detect(complement[1],"[[:upper:]]")) {
                full_name <- paste(c(genus, epiphet), collapse = " ")
              } else {
                complement <- subset(complement, !stringr::str_detect(complement,"[[:upper:]]"))
                complement <- subset(complement, !grepl(paste(c("[()]","&","([0-9]+).*$","^ex$"), collapse="|"), complement))
                if(length(complement)==0){
                  full_name <- paste(c(genus, epiphet), collapse = " ")
                } else {
                  full_name <- paste(c(genus, epiphet, complement), collapse = " ")
                }
              }
            } 
          }
        }  
      }
    }
    }
    results[name_index] <- full_name
  }
  return(results)
}

#########################
organize.bubble.plot2 <- function(points, twgd_data) {
  focal_areas <- as.character(twgd_data$LEVEL3_COD)
  points <- points[,c(2,3)]
  sp::coordinates(points) <- ~ lon + lat
  results <- matrix(nrow=0, ncol=4)
  for(i in 1:length(focal_areas)) {
    #i<-which(focal_areas=="BZL")
    one_area <- focal_areas[i]
    area_plus_buffer <- twgd_data[which(as.character(twgd_data$LEVEL3_COD) %in% one_area),]
    if(nrow(area_plus_buffer)>0) {
      res <- over(points, area_plus_buffer)
      if(any(!is.na(res$LEVEL1_COD))) {
        n_points <- length(which(!is.na(res$LEVEL1_COD)))
        centroids <- rgeos::gCentroid(area_plus_buffer, byid=TRUE)
        lon <- extent(centroids)[1]
        lat <- extent(centroids)[3]
        results <- rbind(results, cbind(n_points, one_area, lon, lat))
      }
      cat(i, "\r")
    }
  }
  results <- as.data.frame(results)
  results$n_points <- as.numeric(results$n_points)
  results$lon <- as.numeric(results$lon)
  results$lat <- as.numeric(results$lat)
  return(results)
}


#########################
organize.bubble.plot <- function(trait_table, reference_table, all_vars, twgd_data) {
  tmp_reference_table <- subset(reference_table, reference_table$wcvp_name %in% unique(trait_table$species))
  wcvp_subset <- subset(all_vars, all_vars$taxon_name %in% tmp_reference_table$wcvp_name)
  wcvp_subset <- subset(wcvp_subset, wcvp_subset$introduced==0)
  wcvp_subset <- subset(wcvp_subset, wcvp_subset$extinct==0)
  wcvp_subset <- subset(wcvp_subset, wcvp_subset$location_doubtful==0)
  
  focal_areas <- unique(wcvp_subset$area_code_l3)
  results <- matrix(nrow=0, ncol=5)
  for(i in 1:length(focal_areas)) {
    one_area <- focal_areas[i]
    one_subset <- subset(wcvp_subset, wcvp_subset$area_code_l3==one_area)
    sp_rich <- length(unique(one_subset$taxon_name))
    family_rich <- length(unique(one_subset$family))
    area_plus_buffer <- twgd_data[which(as.character(twgd_data$LEVEL3_COD) %in% one_area),]
    if(nrow(area_plus_buffer)>0) {
      centroids <- rgeos::gCentroid(area_plus_buffer, byid=TRUE)
      lon <- extent(centroids)[1]
      lat <- extent(centroids)[3]
      results <- rbind(results, cbind(sp_rich, family_rich, one_area, lon, lat))
    }
    cat(i, "\r")
  }
  results <- as.data.frame(results)
  results$sp_rich <- as.numeric(results$sp_rich)
  results$family_rich <- as.numeric(results$family_rich)
  results$lon <- as.numeric(results$lon)
  results$lat <- as.numeric(results$lat)
  return(results)
}


#########################
sum.twgd.trait <- function(one_dataset,twgd_data,all_vars) {
  tmp_rasters <- list()
  for(i in 1:nrow(one_dataset)) {
    wcvp_subset <- subset(all_vars, all_vars$taxon_name == one_dataset$species[i])
    occ_areas <- wcvp_subset$area_code_l3
    area_plus_buffer <- twgd_data[which(as.character(twgd_data$LEVEL3_COD) %in% occ_areas),]
    tmp0 <- raster(area_plus_buffer, res=0.1)
    tmp0[]  <- 1  
    tmp1 <- crop(tmp0, extent(area_plus_buffer))
    tmp2 <- mask(tmp1, area_plus_buffer)
    tmp_rasters[[i]] <- tmp2
    cat(paste0(i, " out of ", nrow(one_dataset)), "\r")
  }
  return(tmp_rasters)
}

#######
map.for.synthesis <- function(trait_table, reference_table, all_vars, twgd_data) {
  trait_table=trait_table
  colnames(trait_table) <- c("species","gbif")
  organized_table_for_plot_total <- organize.bubble.plot(trait_table, reference_table, all_vars, twgd_data)
  organized_table_for_plot_no_data <- organize.bubble.plot(subset(trait_table, trait_table$gbif=="no_data"), reference_table, all_vars, twgd_data)
  
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
  
  # hist(proportion_table$sp_rich_prop)
  # cutting the data

  return(proportion_table)
}
