Codes used to generate Figure 2 of paper "Discovering the rules of plant biogeography using a trait-based approach" (submitted)

Author: Thais Vasconcelos

Data used for the maps was retrieved from:

Seed mass data:
https://bien.nceas.ucsb.edu/bien/ 

BIEN R package:
Maitner et al. (2018). The bien r package: A tool to access the Botanical Information and Ecology Network (BIEN) database. Methods in Ecology and Evolution, 9(2), 373-379.

Ploidy data:
Rice et al. (2019). The global biogeography of polyploid plants. Nature Ecology & Evolution, 3(2), 265-273.

Genetic data for phylogenetic inference:
Smith, S. A., & Brown, J. W. (2018). Constructing a broadly inclusive seed plant phylogeny. American journal of botany, 105(3), 302-314.

Distribution data was kindly made available by:
https://powo.science.kew.org/ [Accessed February 2022]

Taxize package:
Chamberlain, S. A., & Szöcs, E. (2013). taxize: Taxonomic search and retrieval in R. F1000Research, 2, 191. https://doi.org/10.12688/ f1000research.2-191.v1



Short description of methods:

I used the tips of the complete tree of seed plants from Smith and Brown (2018) (i.e. the ALLMB.tre tree) as a reference for a complete list of species names in seed plants. This tree initially included 356,305 tips, and I used functions of the R package taxize to standardize the species name taxonomy according to the GBIF backbone as a first step. After dropping synonyms and unmatched names, 337,473 species names remained. This species list was used as a base for all subsequent analyses. For the first map (Figure 2A), I intersected the names remained in the cleaned ALLMB.tre tree with those present in the cleaned tree from the same publication that only includes species with molecular data (i.e. the GBMB.tre tree). I then used the information about natural species distribution (i.e. excluding invasive species) from the WCVP dataset to set which proportion of species in each TWDG level 3 botanical country were present only in the full species list -- that is, which proportion of species in each botanical country did not have molecular data available for phylogenetic inference at the time when the tree from Smith and Brown (2018) was inferred. I followed the same pipeline to create a list of species that are present in the cleaned ALLMB.tre tree, but that are not available in the ploidy dataset of Rice et al. (2019) nor in the BIEN online dataset (BIEN 2022) of seed mass data. Note that here I also used the R package taxize to standardize the taxonomy of the species in all datasets before intersecting them to avoid creating a mismatch due to outdated taxonomy. As in Figure 2A, I used these lists and the WCVP dataset to map the proportion of species that are in in the cleaned ALLMB.tre tree but not in the cleaned ploidy and seed mass datasets -- that is, species that are native to those botanical countries but for which there was no ploidy nor seed mass data available at the time these datasets were assembled. 

