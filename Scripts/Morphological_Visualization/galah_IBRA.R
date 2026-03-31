#remotes::install_github("AtlasOfLivingAustralia/galah-R@dev")
# requires R >4.3

################################################################################

library(galah)

galah_config(atlas="Australia")
galah_config(email = "iangbrennan@gmail.com")

################################################################################

# identify an IBRA field
View(search_fields("IBRA"))
# let's choose 'cl1048'
View(search_fields("cl1048") |> show_values())
ibra.names <- search_fields("cl1048") |> show_values()

# the filter function is currently broken and can only take a single condition, so we have to do a for-loop
all.spp <- NULL
for (k in 1:length(savannah.ibra)){
  print(savannah.ibra[[k]])
  res <- galah_call() |>
    galah_identify("Squamata") |>
    galah_filter(cl1048 == savannah.ibra[[k]]) |>
    atlas_species()
  
  all.spp <- rbind(all.spp, res)
}
unique.sp <- dplyr::distinct(all.spp)
write.csv(unique.sp, file="~/Desktop/galah_Savannah_Squamata_spp.csv", row.names=F, quote=F)

all.occ <- NULL
for (k in 1:length(savannah.ibra)){
  print(savannah.ibra[[k]])
  res <- galah_call() |>
    galah_identify("Squamata") |>
    galah_filter(cl1048 == savannah.ibra[[k]]) |>
    atlas_occurrences()
  
  all.occ <- rbind(all.occ, res)
}
write.csv(all.occ, file="~/Desktop/galah_Savannah_Squamata_occ.csv", row.names=F, quote=F)

################################################################################

# if the filter function wasn't broken, we could just do it like this
ibra.res <- NULL
for(k in 1:nrow(ibra.names)){
  ibra.res[[k]] <- galah_call() |>
    galah_identify("Agamidae") |>
    galah_filter(cl1048 == ibra.names[[1]][[k]]) |>
    atlas_species()
}
names(ibra.res) <- ibra.names[[1]]

all.ibras <- NULL
for(y in 1:length(ibra.res)){all.ibras <- append(all.ibras, ibra.res[[y]]$species_name)}
ibra.table <- table(all.ibras)
ibra.table[order(ibra.table)]


ibra.res[[1]]

galah_call() |>
  galah_identify("Agamidae") |>
  galah_filter(cl1048 == ibra.names[[1]][[2]]) |>
  atlas_species()

write.csv(res, file="~/Desktop/galah_Savannah_Squamata_spp.csv", row.names=F, quote=F)

rec <- galah_call() |>
  galah_identify("Squamata") |>
  galah_filter(cl1048 %in% savannah.ibra) |>
  atlas_occurrences()
write.csv(rec, file="~/Desktop/galah_Savannah_Squamata_occ.csv", row.names=F, quote=F)

################################################################################

# Let's plot to make sure we've got the whole region covered
library("rnaturalearth")
world <- ne_countries(scale = "medium", returnclass = "sf")
# plot the map 
ggplot(data = world) +
  geom_sf() +
  geom_point(data = all.occ, aes(x = decimalLongitude, y = decimalLatitude), size = 1, 
             shape = 21, fill = "#79cc27", alpha=0.5) +
  coord_sf(xlim = c(110, 160), ylim = c(-5, -45), expand = FALSE) +
  theme_bw()

################################################################################

spp <- read.csv("/Users/ianbrennan/Desktop/Savannah/SavSquam.csv")
spp <- dplyr::filter(spp, remove=="")

occ <- read.csv("/Users/ianbrennan/Desktop/Savannah/galah_Savannah_Squamata_occ.csv")
occ2 <- dplyr::filter(occ, scientificName %in% spp$species_name)

testo <- atlas_counts(occ2, type = "species", # get species counts
             limit = NULL)



################################################################################


library(sf)

## Australia outline
ozout <- read_sf("/Users/ianbrennan/Desktop/Savannah/Australia_SHP/Australia.shp")
ozout <- st_set_crs(ozout, 4236) # set it to WGS84 projection

## Read in the IBRA subregions shapefile
oz <- read_sf("/Users/ianbrennan/Downloads/IBRA7_subregions/IBRA7_subregions.shp")

## Read in the Reptile species shapefile from GARD
gard <- read_sf("/Users/ianbrennan/Downloads/doi_10_5061_dryad_9cnp5hqmb__v20220427/Gard_1_7_ranges.shp")
gard$geometry

# reproject the IBRA regions to match the GARD data
target_epsg <- 4326 # Example: WGS 84 geographic coordinate system
ibra_wgs84 <- st_transform(oz, crs = target_epsg)

# identify the taxa of interest and subset the GARD data down to those species
taxa <- c("Moloch horridus", "Chlamydosaurus kingii")
genera <- c("Physignathus","Hypsilurus","Chelosania","Moloch","Lophosaurus",
            "Intellagama","Cryptagama","Ctenophorus","Chlamydosaurus",
            "Lophognathus","Amphibolurus","Gowidon","Tropicagama",
            "Rankinia","Pogona","Tympanocryptis","Diporiphora")
dragon.gard <- subset(gard, binomial %in% gard$binomial[grep(paste(genera,collapse="|"), gard$binomial)])

# we can plot the range of any individual taxon
ggplot() +
  geom_sf(data=ozout, fill = "#69b3a2", color = "white") +
  geom_sf(data=subset(dragon.gard,binomial=="Ctenophorus slateri"), fill="pink") +
  theme_void()

# if you need to drop a subregion for any reason
#ibra_wgs84 <- subset(ibra_wgs84, SUB_CODE_7 !="GAW08")
#ibra_wgs84 <- subset(ibra_wgs84, SUB_CODE_7 !="SEC02")

# make a dataframe to hold the results
# columns are the IBRA subregions
# rows are the species
dragon.subregions <- NULL
dragon.subregions <- data.frame(matrix(nrow=length(dragon.gard$binomial),
                                       ncol=length(ibra_wgs84$SUB_CODE_7)))
colnames(dragon.subregions) <- ibra_wgs84$SUB_CODE_7
rownames(dragon.subregions) <- dragon.gard$binomial

# run a loop that uses sf::st_intersects to compare each species range
# to the IBRA subregions and stores that information in the dataframe
for(k in 1:length(dragon.gard$binomial)){
  print(paste("taxon:",dragon.gard$binomial[[k]]))
  curr.taxon <- subset(gard, binomial == dragon.gard$binomial[[k]])
  if(!st_is_valid(curr.taxon)==TRUE){curr.taxon <- st_make_valid(curr.taxon)}
  for(j in 1:length(ibra_wgs84$SUB_CODE_7)){
#  for(j in 1:10){
#    print(paste("subregion:", ibra_wgs84$SUB_CODE_7[[j]]))
    curr.region <- subset(ibra_wgs84, SUB_CODE_7 == ibra_wgs84$SUB_CODE_7[[j]])
    inter <- tryCatch({st_intersects(curr.taxon, curr.region)[[1]]},
                      error = function(msg){return(0)})
    # NEED TO RESOLVE THIS
    if(length(inter)==0){inter <- 0}
    dragon.subregions[k,j] <- inter
    #    dragon.subregions[k,j] <- ifelse(inter==1,1,0)
  }
}

# summarize the information by determining how many regions each
# species ranges across
dragon.subregions$range <- rowSums(dragon.subregions)
range.res <- dplyr::select(dragon.subregions, range)
View(range.res)
range.res$taxon <- rownames(range.res); range.res$taxon <- gsub(" ","_",range.res$taxon)

# these are the taxa we have to score by hand
setdiff(agam.tree$tip.label, range.res$taxon)

#############################################

ibra.ranges <- read.csv("/Users/ianbrennan/Documents/GitHub/Amphibolurinae/Data/GARD_IBRA.csv")
ranges <- dplyr::select(ibra.ranges, taxon, corrected_subregions)
ranges <- ranges[complete.cases(ranges),]
