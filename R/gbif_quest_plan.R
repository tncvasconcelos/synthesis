# setwd("~/Desktop/MiSSEgradient/full_gbif_quest")


full_gbif_quest <- drake_plan (
  #base.dir     = "/Users/thaisvasconcelos/Desktop/full_gbif_quest/",
  base.dir = "/home/tvasconcelos/full_gbif_quest/",
  #traqueophyta = data.table::fread(paste0(base.dir, "/0031497-200221144449610.csv")),
  #data         = subset(traqueophyta, !is.na(speciesKey) & !is.na(familyKey) & !is.na(orderKey) & taxonomicStatus=="ACCEPTED" & taxonRank=="SPECIES"),
  #full_names   = sort(organize(data)),
  full_names    = as.character(read.csv(paste0(base.dir, "full_names_0031497-200221144449610.csv"))[,2]),
  gbif_cleaned_all = DoOneTaxon(full_names, base.dir)
)

#base.dir     = "/Users/thaisvasconcelos/Desktop/full_gbif_quest/"
#traqueophyta = data.table::fread(paste0(base.dir, "/0031497-200221144449610.csv"))
#data         = subset(traqueophyta, !is.na(speciesKey) & !is.na(familyKey) & !is.na(orderKey) & taxonomicStatus=="ACCEPTED" & taxonRank=="SPECIES")
#full_names   = sort(organize(data))
#write.csv(full_names, "full_names_.csv")

# unlink(".drake/.gitignore")


gbif_cleaning_plan <- drake_plan (
  base.dir     = "/Users/thaisvasconcelos/Desktop/full_gbif_quest/",
  #base.dir = "/home/tvasconcelos/full_gbif_quest/",
  gbif_file_names = list.files(paste0(base.dir, "gbif_final"), ".csv"),
  crops = resolveGBIF(as.character(read.csv(paste0(base.dir, "crops.csv"),h=F)[,1])),
  invasive = resolveGBIF(as.character(read.csv2(paste0(base.dir, "gisd_database.csv"),h=T)[,1])),
  full_cleaned_gbif = FullCleanGBIF(gbif_file_names, invasive, crops, base.dir)
)

crops_final <- unique(paste0(lapply(strsplit(as.character(crops), " "), "[", 1)," ", lapply(strsplit(as.character(crops), " "), "[", 2)))
invasive_final <- unique(paste0(lapply(strsplit(as.character(invasive), " "), "[", 1)," ", lapply(strsplit(as.character(invasive), " "), "[", 2)))

write.csv(crops_final, file="/Users/thaisvasconcelos/Desktop/full_gbif_quest/crops_taxized_gbif.csv")
write.csv(invasive_final, file="/Users/thaisvasconcelos/Desktop/full_gbif_quest/invasive_taxized_gbif.csv")

saveRDS(crops_final, file="/Users/thaisvasconcelos/Desktop/full_gbif_quest/crops_taxized_gbif.Rdata")
saveRDS(invasive_final, file="/Users/thaisvasconcelos/Desktop/full_gbif_quest/invasive_taxized_gbif.Rdata")



