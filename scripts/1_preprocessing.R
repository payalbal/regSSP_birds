## on boab

## Processing input data for analysis
##
## Outputs:
##  1. min non-NA set masks - raster layer [3 files: til, vn, aus + 2 gsdm files from getData()]
##  2. min non-NA set covariates - raster stack [3 files: til, vn, aus]
##  3. covariates wihtout NA values - raster stack [3 files: til, vn, aus]
##  4. bioq - raster stack [18 files: vn, aus for 0.25, 0.5 and 0.75 
##      quartiles under rcp 45, 60, 85]
##  5  bioq wihtout NA values - raster stack [18 files: vn, aus for 0.25, 
##      0.5 and 0.75 quartiles under rcp 45, 60, 85]
##  6. bioregions - raster layer [1 file: aus]
##  7. bioregions wihtout NA values - raster layer [1 file: aus]
##  8. harvested2016 - data.frame - 8 sectors x 2 columns (sector, proportion harvested)
##      prop of land harvested for each commodity/sector in 2016 [2 files: vn, aus]
##  9. gtap_landendowments - data.frame - 10 sectors x 52 years
##      prop of land harvested 'endowmnet' for each commodity/sector in each 
##      year from 2019 - 2070 [6 files: vn, aus for ssp 1, 2, 3 ]
##  10. occ - data.table of filtered gbif data for Aves 
##      [2 files: vn, aus + 2 text log files]
## 
## TOTAL FILES CREATED IN OUTPUT FOLDER (RData): xxx

## Master log file
job_start <- Sys.time()
masterlog <- paste0("./preprocessing_run.txt")
writeLines(c(""), masterlog)
cat(paste0(">> Job start = ", job_start, "\n"), file = masterlog, append = TRUE)


## Set up work environment ####
setwd("./regSSP_birds/")

rm(list = ls())
gc()
# devtools::install_github('kapitzas/WorldClimTiles')
# devtools::install_github('skiptoniam/sense')
# devtools::install_github('smwindecker/gdaltools')
x <- c('data.table','rgdal','rgeos','matrixStats',"sp",'raster',
       'WorldClimTiles','sense' , 'readxl')
lapply(x, require, character.only = TRUE)
source(file.path(".", "scripts", "0_functions.R"))
gsdms_data <- "../regSSP_data" #"/Volumes/discovery_data/gsdms_data"
regSSP_data <- "../regSSP_data" # "/Volumes/discovery_data/regSSP_data"
rdata_path <- "./RData" # file.path(regSSP_data, "RData_lulcc")
if(!dir.exists(rdata_path)){dir.create(rdata_path)}

## 1. Prepare masks ####
regions <- c("vn", "aus", "til")

for (region in regions){
  
  if (region == "til") {
    reg_mask <- raster(extent(90, 120, 0, 30),
                       crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ",
                       res = 0.008333333)
    reg_mask[] <- 1
    saveRDS(reg_mask, file.path(rdata_path, paste0("mask_", region, ".rds")))
  } else {
    reg_mask <- getData("GADM", country = region, level = 0, path = rdata_path)
    reg_mask <- gSimplify(reg_mask, tol = 0.00833)
    
    if (region == "vn") {
      mask_template <- raster(
        nrow = 1779,
        ncol = 879,
        crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ",
        ext = extent(c(102.1417, 109.4667, 8.566667, 23.39167))
      )
      mask_template[] <- 1
    }
    
    else if (region == "aus") {
      mask_template <- raster(
        nrow = 4091,
        ncol = 4990,
        crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ",
        ext = extent(c(112.4667, 154.05,-44.04167,-9.95))
      )
      mask_template[] <- 1
    }
    reg_mask <- mask(mask_template, reg_mask)
    saveRDS(reg_mask, file.path(rdata_path, paste0("mask_", region, ".rds")))
    plot(reg_mask)
  }
}
rm(mask_template)


## ***** Specify region for analysis (up till Bioclim - future) ***** ####
##  Rerun for c("vn", "aus", "til")
# region <- 'aus'
for (region in regions){
  
  
  ## 2. Prepare covariate data ####
  reg_mask <- readRDS(file.path(rdata_path, paste0("mask_", region, ".rds")))
  
  ## 2a. Topography ####
  ## Source: https://webmap.ornl.gov/ogc/wcsdown.jsp?dg_id=10008_1
  ## Select polygon including vn and aus
  srtm <- raster(file.path(regSSP_data, "sdat_10008_1_20191110_184812575.asc"))
  crs(srtm) <- crs(reg_mask)
  srtm <- projectRaster(srtm, reg_mask)
  srtm <- crop(srtm, reg_mask)
  srtm <- mask(srtm, reg_mask)
  plot(srtm)
  
  names(srtm) <- "srtm"
  elevation <- mask(srtm, reg_mask)
  slope <- terrain(srtm, opt = "slope")
  roughness <- terrain(srtm, opt = "roughness")
  terrain <- stack(elevation, slope, roughness)
  rm(elevation, roughness, slope, srtm)
  
  ## 2b. Soil ####
  ## Source: https://daac.ornl.gov/SOILS/guides/igbp-surfaces.html
  ## See gsdms for data download instructions
  soil <- stack(list.files(file.path(gsdms_data, "orders"), pattern = "*.dat", full.names = T, recursive = T))
  crs(soil) <- crs(reg_mask)
  soil <- projectRaster(soil, reg_mask)
  soil <- crop(soil, reg_mask)
  soil <- mask(soil, reg_mask)
  names(soil) <- c( "bulk", "awco", "carb",  "nitro")
  
  
  ## 2c. Distance rasters ####
  ## Shortest distance rasters calculated were distance to roads, distance to rivers, distance to built-up areas and distance to lakes.
  
  ## Sources:
  ## Global road networks: http://sedac.ciesin.columbia.edu/data/set/groads-global-roads-open-access-v1/data-download#openModal
  ## Global drainage systems: http://www.soest.hawaii.edu/wessel/gshhg/
  ## GLobal built-up areas: http://ref.data.fao.org/map?entryId=c22837d0-88fd-11da-a88f-000d939bc5d8&tab=metadata
  ## Global lakes: http://www.soest.hawaii.edu/wessel/gshhg/
  ## Protected areas: http://wcmc.io/wdpa_current_release (downloaded Feb 2018)
  
  ## Processing steps:
  ## 1. combined shapefiles (different levels) in qgis (Levels 2-9 for rivers)
  ## 2. rasterize in qgis (automatically subsetted from global data set)
  ## 3. calculate proximity in qqis
  ## 4. then load here to mask NAs and combine with rest of data
  
  ## this step is so qgis can read the attributes for some of the croppsed shapefiles
  ## only needs done once
  # name <- "WDBII_rivers_global_L2-L9"
  # region <- "vn"
  # shp <- readOGR("/Users/simon/Dropbox/PhD - Large Files/PhD - Raw Data/Global", layer = name)
  # shp@data[,2] <- as.numeric(shp@data[,2])
  # shp@data[,1] <- as.numeric(shp@data[,1])
  # writeOGR(shp, dsn = "/Users/simon/Dropbox/PhD - Large Files/PhD - Raw Data/Global", layer = name, driver = "ESRI Shapefile", overwrite = T)
  ## load qgis processed files, and set NA
  # lakes <- raster(paste0(raw_rasters_rivers_lakes_cl/lakes_", region, ".tif"))
  # coast <- raster(paste0(raw_rasters_rivers_lakes_cl/coast_", region, ".tif"))
  # rivers <- raster(paste0(raw_rasters_rivers_lakes_cl/rivers_", region, ".tif"))
  # roads <- raster(paste0(raw_rasters_rivers_lakes_cl/PA_", region, ".tif"))
  # builtup <- raster(paste0(raw_rasters_rivers_lakes_cl/PA_", region, ".tif"))
  # names <- c("distrivers", "distlakes", "distcoastline", "PA", "distbuiltup", "distroads")
  
  distances <- stack(list.files(file.path(regSSP_data, "qgis_files"), full.names = T, pattern = region)[-6])
  distances <- crop(distances, reg_mask)
  distances <- mask(distances, reg_mask)
  names(distances) <- c("dibu", "dico", "dila", "diri", "diro")
  
  
  ## 2d. Protected areas #### - to be fixed
  # file_in <- list.files(file.path(gsdms_data, "protectedareas"), pattern= "WDPA_Nov2018-shapefile-polygons.shp", full.names = T, recursive = T)
  # file_out <- file.path(regSSP_data, "qgis_files", basename(file_in))
  # file_out <- gsub(".shp", "_cropped.shp", file_out)
  # crop_shp(file_in, file_out, extent(reg_mask))
  # cropped <- readOGR(file_out)
  # subs <- cropped[which(cropped@data$IUCN_CAT%in%c("II", "Ia", "Ib")),]
  # output_subs <- "WDPA_Mar2018-shapefile-polygons_cropped_subs"
  # writeOGR(subs, file.path(regSSP_data, "qgis_files"), output_subs, driver = "ESRI Shapefile", overwrite_layer = FALSE)
  # 
  # file_in <- paste0(file.path(regSSP_data, "qgis_files", output_subs), ".shp")
  # file_out <- file.path(regSSP_data, "qgis_files", paste0("pa_raster_", region, ".tif"))
  # gdaltools::rasterize_shp(file_in, file_out, res = res(reg_mask), ext = extent(reg_mask))
  # 
  # pa <- raster(file_out)
  # pa <- crop(pa, reg_mask)
  # pa <- mask(pa, reg_mask)
  # pa[] <- (pa[] * -1) + 1
  # names(pa) <- "pa"
  # # writeRaster(pa, paste0(rdata_path, "/pa_", region, ".tif"), format = "GTiff", overwrite = TRUE)
  # rm(cropped, subs, output_subs)
  pa <- list.files(file.path(regSSP_data, "qgis_files"), full.names = T, pattern = region)
  pa <- pa[grepl("PA_", pa)]
  pa <- raster(pa)
  names(pa) <- "pa"
  pa <- crop(pa, reg_mask)
  pa <- mask(pa, reg_mask)
  pa[] <- (pa[] * -1) + 1
  
  
  
  ## 2e. Population density ####
  ## Source: http://sedac.ciesin.columbia.edu/data/set/grump-v1-population-density
  popdens <- raster(file.path(gsdms_data, "popdens", "gluds00ag.bil"))
  popdens <- crop(popdens, reg_mask, snap = "near")
  popdens <- mask(popdens, reg_mask)
  names(popdens) <- "popd"
  
  
  ## 2f. Land use ####
  ## Source: USGS - MODIS Land cover (used for lulcc analysis)
  ## link: https://archive.usgs.gov/archive/sites/landcover.usgs.gov/global_climatology.html
  ## direct link: https://archive.usgs.gov/archive/sites/landcover.usgs.gov/documents/GlobalLandCover_tif.zip
  l <- list.files(file.path(regSSP_data, "copernicus_discrete", region), full.names = TRUE, recursive = FALSE)
  lutiles_list <- list()
  for(i in 1:length(l)){
    print(i)
    r <- raster(l[[i]])
    lutiles_list[[i]] <- projectRaster(r, reg_mask, method = "ngb")
    names(lutiles_list[[i]]) <- "lu"
  }
  
  tiles <- tile_merge(lutiles_list)
  lu <- mask(tiles, reg_mask)
  
  table(lutiles_list[[2]][])
  lu[lu[] == 50] <- 1 #urban
  lu[lu[] == 40] <- 2 #crop
  lu[lu[] == 30] <- 3 #herbacious vegetation
  lu[lu[] == 20] <- 4 #shurbs
  lu[lu[]%in%c(121:126)] <- 5 #open forest
  lu[lu[]%in%c(111:116)] <- 6 #closed forest
  lu[lu[]%in%c(90, 100)] <- 7 #herbaceous wetland, moss, lichen
  lu[lu[]%in%c(60, 70)] <- 8 #bare, sparese, ice
  lu[lu[]%in%c(0, 80, 200)] <- NA #permanent water bodies (covered by distance to river and lake covariates)
  names(lu) <- "landuse"
  plot(lu)
  rm(r, lutiles_list, tiles)
  
  # ## Source: https://earthexplorer.usgs.gov/
  # ## See gsdms script for dowload instructions 
  # world_lu <- raster("~/OneDrive - The University of Melbourne/PhD - Large Files/PhD - Raw Data/Global/LCType.tif")
  # lu <- crop(world_lu, reg_mask)
  # lu <- projectRaster(lu, reg_mask, method = "ngb")
  # lu <- mask(lu, reg_mask)
  # 
  # #Reclassify, see original calsses below.
  # lu[lu[]%in%c(12,14)] <- 100
  # lu[lu[]%in%c(6, 7,9,10,8)] <- 200
  # lu[lu[]%in%c(1,2,3,4,5)] <- 300
  # lu[lu[]%in%c(13)] <- 400
  # lu[lu[]%in%c(11,15,16,0)] <- 500
  # 
  # lu <- lu/100 - 1
  # table(lu[])
  # names(lu) <- "landuse"
  # 0	Water
  # 1	Evergreen Needle leaf Forest
  # 2	Evergreen Broadleaf Forest
  # 3	Deciduous Needle leaf Forest
  # 4	Deciduous Broadleaf Forest
  # 5	Mixed Forests
  # 6	Closed Shrublands
  # 7	Open Shrublands
  # 8	Woody Savannas
  # 9	Savannas
  # 10	Grasslands
  # 11	Permanent Wetland
  # 12	Croplands
  # 13	Urban and Built-Up
  # 14	Cropland/Natural Vegetation Mosaic
  # 15	Snow and Ice
  # 16	Barren or Sparsely Vegetated
  
  
  ## 2g. Bioclim - current ####
  ## Download tiles, mosaic, crop and write
  file_in <- list.files(file.path(gsdms_data, 'bio_30s'), full.names = T)
  bioclim <- list()
  
  for (f in 1:length(file_in)){
    bioclim[[f]] <- crop(raster(file_in[f]), reg_mask)
    bioclim[[f]] <- mask(bioclim[[f]], reg_mask)
  }
  bioclim <- stack(bioclim)
  names(bioclim) <- paste0("bio", 1:19)  
  
  ## Sync NAs - I ####
  ## Find min non-NA set values across mask and covariates and sync NAs 
  covariates <- stack(terrain, soil, distances, pa, popdens, bioclim, lu)
  for(i in 1:nlayers(covariates)){
    reg_mask <- mask(reg_mask, covariates[[i]]) # to find minimum non-NA set values for mask layer
    print(i)
  }
  
  for(i in 1:nlayers(covariates)){
    covariates[[i]] <- mask(covariates[[i]], reg_mask) # mask using updated mask layer
    print(i)
  }
  
  summary(reg_mask) # summary(reg_mask[])
  summary(as.data.frame(covariates))
  saveRDS(readAll(covariates), file = paste0(rdata_path, "/covariates_", region, ".rds"))
  saveRDS(reg_mask, file = paste0(rdata_path, "/mask_", region, ".rds"))
  
  rm(terrain, soil, distances, pa, popdens, bioclim, lu)
  rm(reg_mask, covariates, file_in)
  gc()
}

## 2h. Bioclim - future ####
## Global variables
# regions <- c("vn", "aus")
# rcps <- c("45", "60", "85")
# models <- c("BC", "CC", "GS", "HD", "HE", "IP", "MI", "MR", "MC", "MG", "NO")


# ## Downlaoad GCM model predictions
# urls_bioclim <- c()
# for (model_i in models){
#   for (rcp_j in rcps){
#     fn = paste0(tolower(model_i), rcp_j, 'bi70', '.zip')
#     thisurl = paste0('https://biogeo.ucdavis.edu/data/climate/cmip5/30s/', fn)
#     urls_bioclim = c(urls_bioclim, thisurl)
#   }
# }
# library(bitops)
# library(RCurl)
# test = lapply(urls_bioclim, url.exists)
# all(unlist(test))
# 
# zipdst = file.path(gsdms_data, 'gcm_30s.zip')
# rasterdst = file.path(gsdms_data, 'gcm_30s/')
# if(!dir.exists(rasterdst)) {
#   dir.create(rasterdst)
# } # end if !dir.exists
# for (i in 1:length(urls_bioclim)){
#   download_from_url(urls = urls_bioclim[i], zipdst, rasterdst)
# }
# 
# ## Mask data for study regions & create stacks by region and gcms
# files_gcm <- list.files(file.path(gsdms_data, 'gcm_30s'), full.names = T, recursive = T)
# gcm_reg_path <- file.path(regSSP_data, 'gcm_reg')
gcm_quant_path <- file.path(regSSP_data, 'gcm_quant')
# gcm_reg_path <- file.path(regSSP_data, 'gcm_reg')
# dir.create(gcm_reg_path)
# 
# for(region in regions){
#   for(model_i in models){
#     reg_stack <- list()
#     file_mod <- files_gcm[grepl(model_i, files_gcm)]
#     for(j in 1:length(rcps)){
#       file_mod_rcp <- file_mod[grepl(rcps[j], file_mod)]
#       temp_stack <- list()
#       for(f in 1:length(file_mod_rcp)){
#         reg_mask <- readRDS(file.path(rdata_path, paste0("mask_", region, ".rds")))
#         temp_stack[[f]] <- mask(crop(raster(file_mod_rcp[f]), reg_mask), reg_mask)
#       }
#       reg_stack[[j]] <- readAll(brick(temp_stack))
#     }
#     saveRDS(reg_stack, file = paste0(gcm_reg_path, "/", region, "_", model_i, ".rds"))
#   }
# }
# rm(temp_stack, reg_stack)
# 
# ## Extract cell-wise quartiles across GCM
# quartiles <- c("q1", "q2", "q3")
# 
# for(region in regions){
#   gcm <- list.files(gcm_reg_path, full.names = T, pattern = region)
#   reg_mask <- readRDS(file.path(rdata_path, paste0("mask_", region, ".rds")))
#   r <- reg_mask
#   inds <- which(!is.na(r[]))
#   
#   j<-1
#   for(j in 1:length(rcps)){
#     saveRDS(stack(), file = paste0(gcm_quant_path, "/bio", "q1_", rcps[j], "_", region,  ".rds"))
#     saveRDS(stack(), file = paste0(gcm_quant_path, "/bio", "q2_", rcps[j], "_", region,  ".rds"))
#     saveRDS(stack(), file = paste0(gcm_quant_path, "/bio", "q3_", rcps[j], "_", region,  ".rds"))
#     print(paste0("processing rcp", rcps[j]))
#     for(k in 1:19){
#       print(paste0("processing bioclim var: ", k))
#       bio <- stack()
#       for(i in 1:length(models)){
#         print(paste0("processing model: ", i))
#         dat <- readRDS(gcm[[i]])[[j]]
#         bio <- stack(bio, dat[[k]])
#       }
#       
#       print(paste0("getting quartiles..."))
#       df1 <- na.omit(as.matrix(getValues(bio)))
#       c <-rowQuantiles(df1, probs = c(0.25, 0.5, 0.75))
#       for(m in 1:3){
#         bioclim <- readRDS(file = paste0(gcm_quant_path, "/bio", quartiles[m], "_", rcps[j], "_", region,  ".rds"))
#         r[inds] <- c[,m]
#         names(r) <- paste0("bio", k)
#         saveRDS(readAll(stack(bioclim, r)), file = paste0(gcm_quant_path, "/bio", quartiles[m], "_", rcps[j], "_", region,  ".rds"))
#       }
#     }
#   }
# }
file.copy(list.files(gcm_quant_path, full.names = TRUE), file.path(rdata_path))
rm(list=setdiff(ls(), c("rdata_path", "gsdms_data", "regSSP_data")))
gc()


## Sync NAs - II ### 
## For covariartes + bioq layers
## Find min non-NA set values across mask, covariates and bioq layers and sync NAs 
regions <- c("vn", "aus") ## or c("vn", "aus", "til")
for (region in regions){
  print(paste0("processing ", region, "... "))
  reg_mask <- readRDS(file.path(rdata_path, paste0("mask_", region, ".rds")))
  covariates <- readRDS(file.path(rdata_path, paste0("covariates_", region, ".rds")))
  
  files  <- list.files(file.path(rdata_path), pattern = region, full.names = TRUE)
  bioq <- files[grepl("bioq", files)]
  covs <- files[grepl("covariates", files)]
  files <- c(bioq, covs)
  for(j in 1:length(files)){
    print(paste0("processing file = ", j, " of ", length(files) ,":", files[j]))
    r <- readRDS(files[[j]])
    for(i in 1:nlayers(r)){
      print(paste0("processing layer = ", i, " of ", nlayers(r)))
      reg_mask <- mask(reg_mask, r[[i]])
    }
  }
  saveRDS(reg_mask, file = paste0(rdata_path, "/mask_", region, ".rds"))
  
  for(j in 1:length(files)){
    print(paste0("processing file = ", j, " of ", length(files) ,":", files[j]))
    r <- readRDS(files[[j]])
    for(i in 1:nlayers(r)){
      print(paste0("processing layer = ", i, " of ", nlayers(r)))
      r[[i]] <- mask(r[[i]], reg_mask)
    }
    saveRDS(readAll(r), file = files[[j]])
  }
  print(paste0("summary for ", region, "... "))
  print(paste0("NAs in mask = ", length(which(is.na(reg_mask[])))))
  print(paste0("NAs in stack = ", length(which(is.na(r[[1]][])))))
  ## only check againt first layer of stack r[[1]] (all layers inn stack are the same)
  print(paste0("Which NAs in mask != NAs in layers: ", length(which(is.na(reg_mask[]) != is.na(r[[1]][])))))
  print("=========================")
}

summary(reg_mask) # summary(reg_mask[])
summary(as.data.frame(covariates))

rm(list=setdiff(ls(), c("rdata_path", "gsdms_data", "regSSP_data")))
gc()


# ## Reduce raster size ####
# ##  Remove NA values from covariate stack 
# ##  Mask retained with NA values
# regions <- c("vn", "aus", "til")
# rdata_path2 <- file.path(regSSP_data, "nonatables_lulcc") # change as per server: boab - "./nonatables_lulcc"
# if(!dir.exists(rdata_path2)){dir.create(rdata_path2)}
# 
# for(region in regions){
#   print(paste0("processing ", region, "... "))
#   reg_mask <- readRDS(file.path(rdata_path, paste0("mask_", region, ".rds")))
#   reg_layers <- list.files(rdata_path, full.names = TRUE, pattern = region)
#   ## files[grepl(paste0("(?=.*",region,")"), files, perl = TRUE)]
#   reg_layers <- reg_layers[!grepl("mask", reg_layers)]
#   
#   for(j in 1:length(reg_layers)){
#     print(paste0("processing file = ", j, " of ", length(reg_layers) ,":", reg_layers[j]))
#     reg_stack <- readRDS(reg_layers[j])
#     ind_nona <- which(!is.na(reg_mask[]))
#     
#     print(paste0("# layers in ", reg_layers[j], " = ", nlayers(reg_stack)))
#     new_dat <- data.table()
#     for(i in 1:nlayers(reg_stack)){
#       print(paste0("processing layer = ", i, " of ", nlayers(reg_stack)))
#       r <- reg_stack[[i]]
#       r <- getValues(r)
#       r <- r[ind_nona]
#       new_dat[,paste0("new_dat", "$", names(reg_stack[[i]])) := r]    
#     }
#     saveRDS(new_dat, file = paste0(rdata_path2, basename(reg_layers[j]), "_nona.", "rds"))
#     ## txt or csv files using fwrite much larger!
#   }
# }
# 
# rm(list=setdiff(ls(), c("rdata_path", "gsdms_data", "regSSP_data")))
# gc()

## Save raster stacks as data.table ####



## 3. Bioregions layer (for Australia only) ####
## Source: http://www.environment.gov.au/fed/catalog/search/resource/downloadData.page?uuid=%7B4A2321F0-DD57-454E-BE34-6FD4BDE64703%7D
bioreg <- readOGR(file.path(regSSP_data, "IBRA7_regions", "ibra7_regions.shp"))
reg_mask <- readRDS(file = file.path(rdata_path, "mask_aus.rds"))
bioreg_rast <- reg_mask
bioreg_rast[] <- NA

l <- length(bioreg@polygons)
for(i in 1:l){
  print(i)
  m <- raster::extract(reg_mask, bioreg[i,], cellnumbers = T)
  if(!is.null(m[[1]])){
    bioreg_rast[m[[1]][,1]] <- i
  }
}
bioreg_rast <- mask(bioreg_rast, reg_mask)
saveRDS(bioreg_rast, file.path(rdata_path, "bioregions_aus.rds"))


## 4. GTAP data ####
regions <- c("aus", "vn")
ssps <- paste0("ssp", 1:3)

## Area harvested for crops in 2016 - FAO data
## Data source:http://www.fao.org/faostat/en/#data/QC
## 'Item' in data is coded for GTAP sectors/commodities befor importing in R
harvested <- read.csv(file.path(regSSP_data,"fao_data/FAOSTAT_data_9-18-2018.csv"), header = T) 

for (region in regions){
  if (region == "aus"){
    harv <- harvested[which(harvested$Area == "Australia"),]
  }else if (region == "vn"){
    harv <- harvested[which(harvested$Area == "Viet Nam"),]
  }
  
  ## Estimate sum of area harvested (in ha) for each GTAP sector
  total_bysector <- aggregate(harv$Value, by = list(harv$GTAP.sector), FUN = function(x) sum(na.omit(x)))
  
  ## Estimate relative area harvested (per GTAP sector) relative to total area harvested for all GTAP sectors considered (7 sectors)
  total_bysector[,2] <- total_bysector[,2]/sum(total_bysector[,2])
  colnames(total_bysector) <- c("gtap_sector", "prop_harvested")
  saveRDS(total_bysector, file.path(rdata_path, paste0("harvested2016_", region, ".rds")))
  
  ## Specify additional gtap sectors (i.e. forestry, cattle)
  gtap_sectors <- c(as.character(total_bysector$gtap_sector), "frs", "ctl")
  
  ## Load GTAP land endowmenet by sector data for 2020 - 2071
  ## Source: Tom Kompas (Uni Melb), Van Ha Pham (ANU)
  for (ssp in ssps){
    temp <- as.data.table(read_xls(file.path(regSSP_data, "gtap_data", "ssps",
                                             paste0(region, "_", ssp, ".xls")), sheet = "land"))
    # temp <- temp[X__1 %in% gtap_sectors]
    temp <- temp[...1 %in% gtap_sectors] # on server X__1 is renamed as ...1
    temp[, c(2:6, ncol(temp)) := NULL]
    temp[, c(20:29) := NULL]
    colnames(temp)[1] <- "gtap_sector"
    
    ## Estimate land endowment values for each year from 2020 - 2070 using linear interpolation
    new_yrs <- c(2020:2070)
    new_temp <- matrix(NA, nrow(temp), length(new_yrs))
    for(k in 1:nrow(temp)){
      new_temp[k,] <- approx(y = temp[k,-1], x = as.numeric(colnames(temp)[-1]), xout = new_yrs)$y
    }
    rownames(new_temp) <- temp$gtap_sector
    colnames(new_temp) <- as.character(new_yrs)
    
    ## Add base year column for 2019
    new_temp <- cbind("2019" = 0, new_temp)
    
    ## Save output data
    saveRDS(new_temp, file = file.path(rdata_path, paste0("gtap_landendowments_", region,"_", ssp, ".rds")))
    
    ## print land endowments for 2070
    print(paste0("land endowments for 2070 for ", region, " and ", ssp, ":"))
    print(new_temp[,ncol(new_temp)])
    print("----------------------") 
  }
}


## 7. Biodiversity data ####
## Australia:  https://doi.org/10.15468/dl.khlzmu
## Tile 29: https://doi.org/10.15468/dl.yapqxq
## Viet Nam: https://doi.org/10.15468/dl.nt0ftl
source(file.path(".", "scripts", "0_filter.gbifcsv.R"))
gbif_backbone <- fread(file.path(regSSP_data, "gbif","Taxon.tsv"))
regions = c("aus", "til")
# For vn, we fit models to larger data from til (subset by vnm species) 
vnm_species <- unique(fread(file.path(regSSP_data, "gbif", 
                                      paste0("gbif_aves_vn.csv")), header = T, 
                            na.strings=c("NA", "", " "))$species)

for(region in regions){
  print(paste0("processing ", region, "... "))
  reg_mask <- readRDS(file.path(rdata_path, paste0("mask_", region, ".rds")))
  
  gbif_all <- fread(file.path(regSSP_data, 
                              "gbif", paste0("gbif_aves_", region, ".csv")), 
                    header = T, na.strings=c("NA", "", " "))
  
  filter.gbifcsv(gbif.downloaded.data = gbif_all,
                 gbif.nub.taxonomy = gbif_backbone,
                 subset.gbifnubtaxonomy.byclass = "Aves",
                 output_folder = rdata_path,
                 output_name = paste0("occ_", region, ".rds"),
                 domain.mask = reg_mask,
                 select_fields = c("decimallongitude", "decimallatitude", "species"))
  
  print(paste0("file saved for ", region, ": ", file.path(rdata_path, paste0("occ_", region, ".rds"))))
  
  ## Subset species for tile 29 as per species in Viet Nam
  if(region == "til"){
    dat_til <- readRDS(file.path(rdata_path, paste0("occ_til.rds")))
    til_species <- unique(dat_til$species)
    til_species <- til_species[which(til_species%in%vnm_species)]
    dat_til <- dat_til[dat_til$species%in%til_species,]
    saveRDS(dat_til, file = file.path(rdata_path, paste0("occ_til.rds")))
  }
}


## Log job duration
job_end <- Sys.time()
job_duration = job_end - job_start
cat(paste0(">> Job end = ", job_end, "\n\n"), file = masterlog, append = TRUE)
cat(paste0("Job duration for region ", region, " = ", job_duration, "\n\n"), file = masterlog, append = TRUE)
cat(paste0("Date:\n"), file = masterlog, append = TRUE)
sink(file = masterlog, append = TRUE)
Sys.Date()
sink(NULL)
cat("\n\nSession Info:\n", file = masterlog, append = TRUE)
sink(file = masterlog, append = TRUE)
sessionInfo()
sink(NULL)

rm(list=ls())
gc()


