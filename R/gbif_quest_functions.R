
GetGbifRecordOneTaxonSurviveFailure <- function(full_names){
  final.df <- c()
  return.taxon <- NA
  search.hits <- NA
  slow.time <- 60
  while(is.na(search.hits)[1]) {
    try(search.hits <- rgbif::occ_search(scientificName = taxon, return='data', limit=50000,  hasCoordinate=TRUE))
    if(is.na(search.hits)) { # warning message here
      Sys.sleep(slow.time)
      slow.time <- min(slow.time*1.2, 60*15)
    }
  }
  if((search.hits[1] == "no data found, try a different search")[1]){
    print("No records")
  }else{
    tmp.df <- data.frame(name=search.hits$acceptedScientificName, key=search.hits$key, lat=search.hits$decimalLatitude, long=search.hits$decimalLongitude, basisOfRecord=search.hits$basisOfRecord, issues=search.hits$issues)
    tmp.df.clean <- try(CleanGBIFRecords(tmp.df, base.dir))
    if(!is.null(tmp.df.clean) & class(tmp.df.clean)!="try-error") {
      #We require that each species has, at a minimum, 10 records:
      # if(dim(tmp.df.clean)[1]>9){
      names(tmp.df.clean) <- c("species", "key", "lat", "long", "basisOfRecord", "issue")
      write.table(tmp.df.clean, file=paste0(getwd(), "/gbif/", full_names, ".csv"), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE, append=FALSE)
      return.taxon <- taxon
      #}
    }
  } 
  
  return(return.taxon)
}

# preliminary cleaning
CleanGBIFRecords <- function(gbif.file, base.dir) {
  #Basic idea, is to do everything as we search by taxon:
  subset.mat <- gbif.file
  #Step 1: Parse out duplicates with respect to lat and long:
  subset.mat <- subset.mat[!duplicated(subset.mat[,3]) & !duplicated(subset.mat[,4]),]
  #Step 2: Remove human observations from the list.
  subset.mat <- subset.mat[which(!subset.mat[,5]=="HUMAN_OBSERVATION" & !subset.mat[,5]=="UNKNOWN" & !subset.mat[,5]=="LIVING_SPECIMEN"),]
  return(subset.mat)
}


DoOneTaxon <- function(full_names, base.dir) {
  #full_names <- unname(full_names)
  gbif_cleaned <- lapply(1:length(full_names), function(x) GetGbifRecordOneTaxonSurviveFailure(full_names[x], base.dir))
  gbif_cleaned <- gbif_cleaned[!is.na(gbif_cleaned)]
  #save_gbif_cleaned <- save(gbif_cleaned, file=file_out(paste0("data/", names(individual_clade), "_gbif_clean.Rsave")))
  #PRUNE TREE TO ONES WITH GBIF CLEANED 
  gbif_cleaned_prefix <- paste0(paste0(full_names, "__"), gbif_cleaned)
  return(gbif_cleaned_prefix)
}


organize <- function(data){
  species      = data[,"scientificName"]
  order        = data[,"order"]
  family       = data[,"family"]
  full_names <- paste(order$order, family$family, species$scientificName, sep="_")
  return(full_names)
}
  


resolveGBIF <- function(name) {
  gnr_resolve_x <- function(x) {
    sources <- taxize::gnr_datasources()
    tmp.name <- suppressWarnings(taxize::gnr_resolve(names=x, data_source_ids=sources$id[sources$title == "GBIF Backbone Taxonomy"], best_match_only=FALSE)$matched_name)
    ########################################################
    # workaround to find plant names in case of conflict: #
    if(length(tmp.name) > 1) { # If two or more names are matched...
      splitted_name <- strsplit(tmp.name, " ")
      year <- tail(strsplit(tmp.name, " ")[[1]], 1) 
      if(!is.na(as.numeric(year))){ # ... and the last character string of the name is a year (NCBI uses author and year of publication to identify conflicting names)
        possible.names <- sapply(1:length(splitted_name), function(x) paste0(splitted_name[[x]][-length(splitted_name[[x]])], collapse = " ")) 
        new.name.plant <- taxize::gnr_resolve(names=x, data_source_ids=sources$id[sources$title == "The International Plant Names Index"], best_match_only=TRUE)$matched_name
        new.name <- tmp.name[which(possible.names %in% new.name.plant)]
      } else { new.name <- tmp.name[1] } # else uses best match from first search
    } else { new.name <- tmp.name[1] } # else uses best match from first search
    ########################################################
    if(is.null(tmp.name)) {
      new.name <- paste0("UNMATCHED_",x)
    }
    return(new.name) 
  }
  new.names <- pbapply::pbsapply(name, gnr_resolve_x)
  return(new.names)
}


#report <- file(paste0(polygon.dir,"/","report.txt"))
#writeLines(paste0(nrow(few.points)," species had less points than the 'threshold' and ", nrow(too.wide), " species have background polygons that are possibly too wide for a natural distribution. Check 'few.points.txt' and 'too.wide.txt' for a full list"), report)
#suppressWarnings(close(report))

FullCleanGBIF <- function(gbif_file_names, invasive, crops, base.dir) {
  import.fail <- c()
  for(i in 1:length(gbif_file_names)) {
    # reading files
    subset.mat = as.data.frame(fread(paste0(base.dir, "gbif_final/",gbif_file_names[i]), h=F))
    label = sub(".csv", "", gbif_file_names[i])
      if(nrow(subset.mat) == 0) { next }
    if(grepl("\t", paste(subset.mat[,1:ncol(subset.mat)], collapse=" "))) {
      import.fail <- c(import.fail, gbif_file_names[i])
      print(gbif_file_names[i])
      next
      #splited = strsplit(subset.mat$V1, "\t")
      #subset.mat$V5 <- paste(sapply(splited, "[",6), 
      #                       subset.mat$V2,
      #                       subset.mat$V3,
      #                       subset.mat$V4,
      #                       subset.mat$V5,
      #                       subset.mat$V6, sep=" ")
      #subset.mat$V1 <- sapply(splited, "[",1)
      #subset.mat$V2 <- sapply(splited, "[",2)
      #subset.mat$V3 <- as.numeric(sapply(splited, "[",3))
      #subset.mat$V4 <- as.numeric(sapply(splited, "[",4))
    #}
    #if(grepl("\t",subset.mat$V2)[1]) {
    #  splited = strsplit(subset.mat$V2, "\t")
    #  subset.mat$V5 <- paste(sapply(splited, "[",6), subset.mat$V6, subset.mat$V7) 
    #  subset.mat$V1 <- paste(subset.mat$V1, sapply(splited, "[",1), sep=" ")
    #  subset.mat$V2 <- sapply(splited, "[",2)
    #  subset.mat$V3 <- as.numeric(sapply(splited, "[",3))
    #  subset.mat$V4 <- as.numeric(sapply(splited, "[",4))
    }
    if(!is.numeric(subset.mat$V4)) {
      import.fail <- c(import.fail, gbif_file_names[i])
      print(gbif_file_names[i])
      next
    }
    # removing centroids
    subset.mat = CoordinateCleaner::cc_cen(subset.mat, lon="V4", lat="V3", species ="V1", geod=T, buffer=25000, verbose=F)
      if(nrow(subset.mat) == 0) { next }
    # only preserved specimens
    subset.mat = subset.mat[subset.mat$V5=="PRESERVED_SPECIMEN",]
      if(nrow(subset.mat) == 0) { next }
    # removing points with no decimal cases
    no.decimal.rows = which(!is.decimal(subset.mat[,"V3"]) | !is.decimal(subset.mat[,"V4"]))
      if(length(no.decimal.rows) > 0) {  subset.mat = subset.mat[-no.decimal.rows, ]}
      if(nrow(subset.mat) == 0) { next }
    # removing outliers in distribution
    species <-  unique(subset.mat$V1) 
    cleaned_point <- data.frame()
    for(cleaning in 1:length(species)){
      sp0 <- subset.mat[subset.mat$V1==species[cleaning],]
      out_lat <- boxplot.stats(sp0$V3)$out
      out_lon <- boxplot.stats(sp0$V4)$out
      sp0 <- sp0[ ! sp0$V3 %in% out_lat, ]
      sp0 <- sp0[ ! sp0$V4 %in% out_lon, ]
      cleaned_point <- rbind(cleaned_point, sp0)
    }
    subset.mat = cleaned_point
      if(nrow(subset.mat) == 0) { next }
    # removing points that are likely the location of herbaria
    herbaria_localities <- read.csv(paste0(base.dir,"/allHerbaria_ADM1_badCoords.txt"), header=FALSE)
    # herbaria_localities <- read.csv("~/Desktop/MiSSEgradient/MiSSEGradient/Sampling/data/allHerbaria_ADM1_badCoords.txt", header=FALSE)
    for(herb.index in 1:length(herbaria_localities)){
      subset.mat <- subset.mat[which(!round(subset.mat$V3,2)==herbaria_localities[herb.index,1] | !round(subset.mat$V4,2)==herbaria_localities[herb.index,2]),]
    }
      if(nrow(subset.mat) == 0) { next }
    # crops and invasive species based on two datasets
    crop_rows = which(subset.mat$V1 %in% as.character(crops))
    if(length(crop_rows) > 0) { subset.mat = subset.mat[- crop_rows, ] }
      if(nrow(subset.mat) == 0) { next }
    invasive_rows = which(subset.mat$V1 %in% as.character(invasive))
    if(length(invasive_rows) > 0) { subset.mat = subset.mat[- invasive_rows, ] }
      if(nrow(subset.mat) == 0) { next }
    # remove points based on issues. Removed any records with the following codes: "cdout", "cdrepf", etc.
    tmp.vector <- c()
    if(dim(subset.mat)[1] > 0){
      issues.list.by.record <- strsplit(as.character(subset.mat$V6), ",")
      if(!all(is.na(sapply(issues.list.by.record, "[", 1)))) {
      for(record.index in 1:dim(subset.mat)[1]) {
        if(any(issues.list.by.record[[record.index]] == "ccm" | issues.list.by.record[[record.index]] == "conti" | issues.list.by.record[[record.index]] == "cdiv" | issues.list.by.record[[record.index]] == "cdout" | issues.list.by.record[[record.index]] == "cdrepf" | issues.list.by.record[[record.index]] == "cdreps" | issues.list.by.record[[record.index]] == "cucdmis" | issues.list.by.record[[record.index]] == "cuiv" | issues.list.by.record[[record.index]] == "cum" |  issues.list.by.record[[record.index]] == "gdativ" | issues.list.by.record[[record.index]] == "reneglat" | issues.list.by.record[[record.index]] == "reneglon" | issues.list.by.record[[record.index]] == "preswcd" | issues.list.by.record[[record.index]] == "zerocd" | issues.list.by.record[[record.index]] == "txmatfuz" | issues.list.by.record[[record.index]] == "txmathi" | issues.list.by.record[[record.index]] == "txmatnon")){
        tmp.vector <- tmp.vector
        } else {
        tmp.vector <- c(tmp.vector, record.index)
        }
      }
    }    
      if(nrow(subset.mat) == 0) { next }
    final.mat <- subset.mat[tmp.vector,]
    write.csv(final.mat, file=paste0(base.dir,"z_filtered_gbif/",label,"_filtered.csv"))
      }
    }
  return(final.mat)
}  
  


