## Functions

## Download bioclim data from url
download_from_url <- function (urls, zipdst, rasterdst) {
  
  for (url in urls){
    response = tryCatch({
      download.file(url, destfile = zipdst)
    }, 
    error = function(e){e}) # shouldn't get called
    
    print(response) # should be 0
    unzip(zipdst, exdir = rasterdst)
    file.remove(zipdst)
  }
}


## Update mask based on NAs in raster stack (to obtain min set NAs)
## Add to function to include masking raster stack based on updated mask
align.maskNA <- function(raster_stack, region_mask) {
  print(paste0("# NAs in input mask: ", summary(region_mask)[6]))
  
  for (i in names(raster_stack)){
    if (!sum(is.na(region_mask@data@values)) == sum(is.na(raster_stack[[i]]@data@values))) {
      nona <- which(!is.na(values(region_mask))) # non-na values in mas
      nas <- which(is.na(values(raster_stack[[i]])[nona])) # na values in covariate
      xys <- xyFromCell(region_mask, nona)
      xys <- xys[nas,]
      values(region_mask)[cellFromXY(region_mask, xys)] <- NA
    }
  }
  print(paste0("# NAs in output mask: ", summary(new_mask)[6]))
  
  return(region_mask)
}


## Download current bioclim data by specified tiles - by Simon Kapitza
get_wctiles <- function(tiles, var, path, ...){
  
  if(missing(path)){
    path <- getwd()
  }
  
  rs <- raster(nrows = 5, ncols = 12, xmn = -180, xmx = 180,
               ymn = -60, ymx = 90)
  rs[] <- 1:length(rs)
  
  tiles_names <- c(paste0(0, 0:11), paste0(1, 0:11), paste0(2, 0:11), paste0(3, 0:11), paste0(4, 0:11))
  
  points <- xyFromCell(rs, which(tiles_names%in%tiles))
  
  message("Loading required worldclim tiles...")
  biotiles <- list()
  for (i in 1:nrow(points)){
    biotiles[[i]] <- getData(name = "worldclim", var = var, path = path, res = 0.5, lon = points[i,1], lat = points[i,2], ...)
  }
  message(paste0("Data successfully downloaded to ", path))
  biotiles
}

## Merge downloaded bioclim tiles into single raster layer - by Simon Kapitza
merge_wctiles <- function(biotiles){
  
  if(!is.list(biotiles)) stop("Please provide list with stacks for each tile")
  
  if (length(biotiles) == 1){
    out <- biotiles[[1]]
    message("Only 1 tile detected, merging not necessary")
  }
  
  if(length(biotiles) > 1){
    bionames <- names(biotiles[[1]])
    out <- list()
    message("Merging tiles...")
    for (i in 1:nlayers(biotiles[[1]])){
      b <- biotiles[[1]][[i]]
      message(paste0(substr(bionames[i], 1, nchar(bionames[i])-3)," | " , i, "/", length(bionames)))
      
      for(j in 2:length(biotiles)){
        b <- raster::merge(b, biotiles[[j]][[i]])
      }
      
      out[[i]] <- b
    }
    out <- stack(out)
    names(out) <- bionames
  }
  out
}

## Reduce predictor set by dropping predictors with highest corrleation with another predictor - by Simon Kapitza
## Same function used in the landuse model
correlations <- function(covs, thresh = 0.7) {
  subs_cor <- sample(1:nrow(covs), size = 15000) # subset covs to calculate correlations
  cors <- cor(covs[subs_cor, ], method = "spearman")
  while (min(abs(cors[abs(cors) >= thresh])) != 1){
    values <- cors[which(abs(cors) > thresh)]
    # corellated <- which(abs(cors) > thresh)
    values[values ==1] <- NA
    # corellated[which(values== max(values, na.rm = T))]
    rows_highest_cor <- which(cors == max(values, na.rm = T), arr.ind = T)[,1]
    cors_cur <- abs(cors[rows_highest_cor,])
    '%ni%' <- Negate('%in%')
    m1 <- max(cors_cur[1,][cors_cur[1,]%ni%c(max(values, na.rm = T),1)])
    m2 <- max(cors_cur[2,][cors_cur[2,]%ni%c(max(values, na.rm = T),1)])
    out <- ifelse(m1 > m2, 1, 2)
    cors <- cors[-which(colnames(cors) == names(rows_highest_cor)[out]), 
                 -which(colnames(cors) == names(rows_highest_cor)[out])]
    nrow(cors)
  }
  return(cors)
}

## Write raster layers from objects to disk, so maxent can pull them from there rather than having them all loaded into memory - by Simon Kapitza
writeToDisk <- function(covariates, folder){
  dir.create(folder, recursive = T)
  writeRaster(covariates, filename = file.path(folder, paste0(names(covariates), ".tif")), bylayer = T, driver = "GTiff", overwrite = T)
}


## Calculate area-corrected richness for each raster in the stack
## NOTES Area-corrected richness is the ratio of value in a pixel vs the sum of values across the area for a species. It gives a relative value that can be compared across species (unlike the relative likelihoods). So we can sum this for species. This is still an index of area corrected richness and not a measure of the richness directly. 

area_corrted_richnes <- function (inputstack){
  output <- stack()
  if(class(inputstack) != "RasterStack"){
    return(NA)
  } else {
    for (i in 1:dim(inputstack)[3]){
      output <- stack(output, inputstack[[i]]/cellStats(inputstack[[i]], stat = "sum", na.rm = T))
      ## values come out the same as inputstack[[i]]@data@values/sum(inputstack[[i]]@data@values, na.rm = T)
      ## check by subsetting i=1 and comparing min and max with that of resultant raster
    }
    return(output)
  }
}