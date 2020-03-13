## Maxent SDMs using dismo
## Outputs...


## DPREPARE SDM DATA ####
## Set up work environment ####
setwd("./regSSP_birds/")

options(java.parameters = "-Xmx200g")
## 200g is size specified for RAM
## https://www.ibm.com/support/pages/error-watson-studio-local-r-studio-javalangoutofmemoryerror-java-heap-space
## https://stackoverflow.com/questions/41950018/rstudio-error-in-jarraym-java-lang-outofmemoryerror-java-heap-space

# ## Load older versions of doParallel, iterators annd foreach
# packageurl <- "https://cran.r-project.org/src/contrib/Archive/doParallel/doParallel_1.0.14.tar.gz"
# packageurl <- "https://cran.r-project.org/src/contrib/Archive/iterators/iterators_1.0.10.tar.gz"
# packageurl <- "https://cran.r-project.org/src/contrib/Archive/foreach/foreach_1.4.4.tar.gz"
# install.packages(packageurl, repos=NULL, type="source")

x <- c("sp", "raster", "data.table", "doParallel", "dismo", "rJava")
lapply(x, require, character.only = TRUE)
sessionInfo()

## If java doens't initialise correctly:
## https://support.apple.com/kb/DL1572?viewlocale=en_US&locale=en_US
## Error in .jnew("RectangularArrayBuilder", .jcast(array), dim) : 
##  java.lang.OutOfMemoryError: Java heap space
## to load rJava start RStudio from shell, using: LD_LIBRARY_PATH=$(/usr/libexec/java_home)/jre/lib/server: open -a 
## options("java.home"="/Library/Java/JavaVirtualMachines/jdk-11.0.1.jdk/Contents/Home/jre")

## Maxent trobleshooting
maxent()
## If maxent() shows error:
## Error in .getMeVersion() : file missing:
## /home/payalb/R/x86_64-pc-linux-gnu-library/3.6/dismo/java/maxent.jar.
## Please download it here: http://www.cs.princeton.edu/~schapire/maxent/
file.exists("/home/payalb/R/x86_64-pc-linux-gnu-library/3.6/dismo/java/maxent.jar") # maxent.3.3.3
list.files("/home/payalb/R/x86_64-pc-linux-gnu-library/3.6/dismo/java/")

# ## https://github.com/johnbaums/rmaxent
# devtools::install_github('johnbaums/rmaxent')
# library(rmaxent)
# get_maxent(version = "latest", quiet = FALSE)
# file.exists("/home/payalb/R/x86_64-pc-linux-gnu-library/3.6/dismo/java/maxent.jar") # maxent.3.3.3
# list.files("/home/payalb/R/x86_64-pc-linux-gnu-library/3.6/dismo/java/")
# 
source(file.path(".","scripts", "0_functions.R"))
regSSP_data <- "../regSSP_data" # '/Volumes/discovery_data/regSSP_data' 
rdata_path <- "./RData" # file.path(regSSP_data, "RData_lulcc") 
output_path <- "./output" # file.path(regSSP_data, "output_lulcc") 
files <- list.files(rdata_path, full.names = TRUE)

## Global parameters ####
regions <- c('til', 'aus') # models run for tile 29 which encompasses Vietnam, and Australia
covs_bias <- c("dibu", "diro", "pa", "popd", "roughness")
covs_exclude <- c("srtm", "popd", "dico", "carb", "nitro", "awco", "bulk", "landuse") # not used in SDM
# covs_exclude <- c(covs_exclude, "soilnitro") # might be used later? not in covariates.

## ***** loop through each region here ***** ####
for(region in regions){
  
  ## Log file
  job_start <- Sys.time()
  masterlog <- paste0("./sdm_data_run_", region, ".txt")
  writeLines(c(""), masterlog)
  cat(paste0(">> Job start = ", job_start, "\n"), file = masterlog, append = TRUE)
  
  ## Load data ####
  covariates <- files[grepl(paste0("covariates_", region), files)]
  covariates <- readRDS(covariates)
  
  ## Create bias layer to account for sampling bias ####
  ## Bias layer predicted based on covs_bias
  birds <- readRDS(file.path(rdata_path, paste0("occ_", region, ".rds")))
  birds <- as.data.frame(birds)
  colnames(birds)[1:2] <- c("long", "lat")
  subs <- sample(1:nrow(birds), size = 5000)
  obs <- birds[subs, c(1,2)]
  covariates_bias <- covariates[[which(names(covariates)%in%covs_bias)]]
  print(paste0("Modelling bias layer for region ", region, " ..."))
  mod <- dismo::maxent(x= covariates_bias, p= obs,nbg = 10000, factors = "pa", 
                       args= c("randomtestpoints=25","-J", "-P", "-p", "-h", "threshold=FALSE"))
  print(paste0("Predicting bias layer for region ", region, " ..."))
  bias <- dismo::predict(mod, covariates_bias)
  saveRDS(bias, file = file.path(output_path, paste0("bias_", region, ".rds")))
  rm(birds, subs, obs, covariates_bias, mod, bias)
  
  
  ## Discard correlated covariates ####
  covariates_subset <- covariates[[-which(names(covariates)%in%c(covs_exclude, covs_bias))]]
  cov_df <- getValues(covariates_subset)
  cov_df <- na.omit(cov_df)
  preds <- colnames(correlations(cov_df))
  saveRDS(preds, file = file.path(output_path, paste0("preds_", region, ".rds")))
  rm(covariates_subset, cov_df, preds)
  
  
  ## For Australia only: Bioregion layer ####
  ## Create list of adjacent bioregions for each unique bioregion (for background sampling)
  if(region == "aus"){
    bio_regs <- readRDS(file.path(rdata_path, "bioregions_aus.rds"))
    ecoregs <- na.omit(unique(bio_regs[]))
    ecoregs <- sort(ecoregs)
    adjs_list <- list()
    for (i in ecoregs){
      adjs <- adjacent(bio_regs, which(bio_regs[] == i), pairs = FALSE, include = TRUE)
      adjs_list[[i]] <- unique(bio_regs[adjs], na.rm = TRUE)
      print(i)
    }
    saveRDS(adjs_list, file = file.path(output_path, "adjs_list.rds"))
  }
  rm(bio_regs, ecoregs, adjs_list, i, adjs)
  
  
  ## Sample background points using bias layer ####
  if (region == "til"){
    bias_rast <- readRDS(file.path(output_path, paste0("bias_", region, ".rds")))
    inds <- which(!is.na(bias_rast[]))
    probs <- round(bias_rast[inds],3)
    inds_all <- sample(inds, p = probs, size = 10000)
    saveRDS(inds_all, file = file.path(output_path, paste0("bgp_", region, "_bias.rds")))
  }
  rm(bias_rast, inds, probs, inds_all)
  
  ## For Australia only: background points are constrained to the species bioregion plus one adjacent
  if(region == "aus"){
    bio_regs <- readRDS(file.path(rdata_path, "bioregions_aus.rds"))
    bias_rast <- readRDS(file.path(output_path, paste0("bias_", region, ".rds")))
    ## Min set non-NA values
    bias_rast <- mask(bias_rast, bio_regs)
    bio_regs <- mask(bio_regs, bias_rast)
    probs <- round(bias_rast[],3)
    aves_obs <- readRDS(file.path(rdata_path, paste0("occ_",region,".rds")))
    
    species <- sort(unique(aves_obs$species))
    biases <- c(TRUE, FALSE)
    adj_list <- readRDS(file.path(output_path, "adjs_list.rds"))
    for (bias in biases){
      samples_list <- list()
      for (i in 1:length(species)){
        subs <- aves_obs[as.character(aves_obs$species)%in%species[i],]
        spec_bio_regs <- na.omit(unique(extract(bio_regs, subs[,c(1,2)])))
        adj_bioregs <- na.omit(unique(unlist(adj_list[spec_bio_regs])))
        inds <- which(bio_regs[]%in%adj_bioregs)
        sample_size <- 10000
        if (length(inds) > 1000000){
          inds <- sample(inds, size = 1000000)
        }
        if (length(inds) < 10000){
          sample_size <- length(inds)
        }
        if(bias == TRUE){
          p <- round(probs[inds], 3)
          samples <- sample(x=inds, size = sample_size, p = p, replace = TRUE)
        }else if(bias == FALSE){
          samples <- sample(x=inds, size = 10000)
        }
        samples_list[[i]] <- samples
        print(i)
      }
      if(bias == TRUE){
        saveRDS(samples_list, file = file.path(output_path, "bgp_aus_bias_constrained.rds"))
      }else if(bias == FALSE){
        saveRDS(samples_list, file = file.path(output_path, "bgp_aus_random_constrained.rds"))
      }
    }
  }
  
  ## Log job duration
  masterlog <- paste0("./sdm_data_run_", region, ".txt")
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
}

rm(list=ls())
gc()



## RUN SDMs ####
## ***** loop through each region here ***** ####
regions <- c("til", "aus") # models run for tile 29 which encompasses Vietnam, and Australia
for (region in regions) {
  
  # system("ps")
  # system("pkill -f R")
  .rs.restartR()
  
  ## Set up work environment ####
  setwd("./regSSP_birds/")
  
  options(java.parameters = "-Xmx200g")
  x <- c("sp", "raster", "data.table", "doParallel", "dismo", "rJava")
  lapply(x, require, character.only = TRUE)
  
  source(file.path(".","scripts", "0_functions.R"))
  regSSP_data <- "../regSSP_data" # '/Volumes/discovery_data/regSSP_data' 
  rdata_path <- "./RData" # file.path(regSSP_data, "RData_lulcc") 
  output_path <- "./output" # file.path(regSSP_data, "output_lulcc") 
  files <- list.files(rdata_path, full.names = TRUE)
  
  
  ## Global parameters ####
  ssps <- paste0("ssp", 1:3)
  rcps <- c("45", "60", "85")
  quartiles <- c("q2", "q1", "q3")
  scens <- sort(apply(expand.grid(quartiles, ssps), 1, paste0, collapse="_"))
  scens_rcps <- sort(apply(expand.grid(quartiles, rcps), 1, paste0, collapse="_"))
  treatments <- c("pre","agg") # treatments <- c("pre", "ind", "dir", "agg")
  
  ## Log file
  job_start <- Sys.time()
  masterlog <- paste0("./sdm_model_run_", region, ".txt")
  writeLines(c(""), masterlog)
  cat(paste0(">> Job start = ", job_start, "\n"), file = masterlog, append = TRUE)
  
  
  ## Load data ####
  ## Observation data
  obs <- readRDS(file.path(rdata_path, paste0("occ_",region,".rds")))
  species_list <- sort(as.character(unique(obs$species)))
  preds <- readRDS(file.path(output_path, paste0("preds_",region,".rds")))
  
  ## Covariates for prediction, write into temp folder structure
  ##  Note: Covariates used for model building and predictions 
  ##    are different for Vietnam, see Methods
  temp <- file.path(output_path, paste0("temp_", region))
  if(region == "til"){
    region_preds <- "vn"
  }else{
    region_preds <- region
  }
  type <- c("covariates", "predlu", "bio")
  for(k in 1:3){
    if(type[k]%in%type[c(2,3)]){
      for(i in 1:length(scens)){
        if(type[k] == "predlu"){
          path <- output_path
          pattern <- paste0(type[k], "_", scens[i], "_", region_preds)
        }
        if(type[k] == "bio"){
          path <- rdata_path
          pattern <- paste0(type[k], scens_rcps[i],"_", region_preds)
        }
        covariates <- readRDS(list.files(path, pattern = pattern, full.names = TRUE))
        covariates <- covariates[[which(names(covariates)%in%c(preds, "landuse"))]]
        writeToDisk(covariates, folder = file.path(temp, pattern))
      }
    }else{
      covariates <- readRDS(list.files(rdata_path, 
                                       pattern = paste0("covariates_", region_preds), 
                                       full.names = TRUE))
      covariates <- covariates[[which(names(covariates)%in%c(preds, "landuse"))]]
      writeToDisk(covariates, folder = file.path(temp, "pres"))
    }
  }
  
  ## Covariates for modelling, write into temp folder structure
  covariates <- readRDS(list.files(rdata_path, pattern = paste0("covariates_", region), full.names = TRUE))
  covariates <- covariates[[which(names(covariates)%in%c(preds, "landuse"))]]
  writeToDisk(covariates, folder = file.path(temp, "models"))
  covariates <- stack(list.files(file.path(temp, "models"), full.names = TRUE))
  
  ## Bias layer and list of predictors to retain
  bias_rast <- readRDS(file.path(output_path, paste0("bias_", region, ".rds")))
  
  ## Background points
  if(region == "aus"){
    bio_region <- readRDS(file.path(rdata_path, "bioregions_aus.rds"))
    samples_list <- readRDS(file.path(output_path, "bgp_aus_bias_constrained.rds"))
  } else {
    inds <- readRDS(file.path(output_path, paste0("bgp_", region, "_bias.rds")))
    bgp <- SpatialPoints(xyFromCell(bias_rast, cell = inds))
    bg <- extract(covariates, bgp)
  }
  
  
  ## Build Maxent model ####
  ## Specify output folder and log file
  output <- file.path(output_path, paste0("models_", region))
  if(!dir.exists(output)){dir.create(output)}
  factors <- "landuse"
  l <- length(species_list)
  predictions <- TRUE
  
  
  ## Load Cluster
  cl <- makeCluster(70)
  registerDoParallel(cl)
  # clusterCall(cl, function() library(sp))
  # clusterCall(cl, function() library(raster))
  # clusterCall(cl, function() library(data.table))
  # clusterCall(cl, function() library(dismo))
  # clusterCall(cl, function() library(rJava))
  # clusterCall(cl, function() library(foreach))
  # clusterCall(cl, function() library(iterators))
  # clusterCall(cl, function() library(parallel))
  # clusterCall(cl, function() library(doParallel))
  
  
  results <- foreach(i = 1:l,
                     .packages = c('sp', 'raster', "data.table", 'doParallel', 
                                   "dismo", "rJava",'foreach', 'iterators')) %dopar% {
                                     
                                     logfile <- file.path(output, paste0(species_list[i], "_log.txt"))
                                     writeLines(c(""), logfile)
                                     cat(paste0("Log file for ", species_list[i], "\n\n"), file = logfile, append = TRUE)
                                     cat("#########################################  \n\n", file = logfile, append = TRUE)
                                     job_species_start <- Sys.time()
                                     cat(paste0(">> Species job start = ", job_species_start, "\n\n"), file = logfile, append = TRUE)
                                     
                                     ## Load covariate and species data
                                     cov <- covariates
                                     subs <- obs[obs$species == species_list[i],]
                                     
                                     ## Load background data
                                     if(region == "aus"){
                                       inds <- samples_list[[i]]
                                       bgp <- SpatialPoints(xyFromCell(bias_rast, cell = inds))
                                       bg <- extract(cov, bgp)
                                     }
                                     
                                     ## Prepare data for model building
                                     sp <- SpatialPoints(subs[,c(1,2)])
                                     pr <- extract(cov, sp)
                                     occ <- c(rep(0, 1, nrow(bg)), rep(1, 1, nrow(pr)))
                                     cov_df <- data.frame(rbind(bg, pr))
                                     
                                     ## Model 1: Determine covariate importance ####
                                     ## Initial model run
                                     cat(paste0("Start model fitting .... ", "\n"), file = logfile, append = TRUE)
                                     model_1_start <- Sys.time()
                                     cat(paste0("Covariate importance model start (model 1) = ", model_1_start, "\n"), file = logfile, append = TRUE)
                                     
                                     mod <- tryCatch(dismo::maxent(cov_df, 
                                                                   p= occ, 
                                                                   factors = factors, 
                                                                   args= c("randomtestpoints=25",
                                                                           "-J", "-p", "-h", 
                                                                           "threshold=FALSE")),
                                                     error = function(e) NA)
                                     
                                     ## Determine covariates with perm importance < 1 and remove from covariate set (see Methods)
                                     if(class(mod) != "MaxEnt"){
                                       cat(paste0("Model 1 failed to fit ... \n"), file = logfile, append = TRUE)
                                       return(NA)
                                     } else {
                                       
                                       res <- data.frame("names" = as.character(rownames(mod@results)), 
                                                         "results" = as.numeric(mod@results))
                                       perm <- res[grep("permutation.importance", res$names),]
                                       kickout <- perm$names[which(perm$results < 1)]
                                       kickout <- gsub('.permutation.importance', '', kickout)
                                       
                                       if(length(kickout) > 0){
                                         cov_df <- cov_df[,-which(colnames(cov_df)%in%kickout), drop = FALSE]
                                         names.new <- colnames(cov_df)
                                       }else{
                                         cov_df <- cov_df
                                         names.new <- colnames(cov_df)
                                       }
                                       
                                       factors.new <- factors[factors%in%names(cov_df)]
                                       cat(paste0("Landuse covariate discarded : ", 
                                                  length(factors.new) == 0, "\n"), file = logfile, append = TRUE)
                                       
                                       model_1_end <- Sys.time()
                                       cat(paste0("Covariate importance model end (model 3) = ", model_1_end, "\n\n"), file = logfile, append = TRUE)
                                       
                                       ## Model 2: Cross-validated model ####
                                       model_2_start <- Sys.time()
                                       cat(paste0("Cross-validation model start (model 2) = ", model_2_start, "\n"), file = logfile, append = TRUE)
                                       
                                       mod.new <- tryCatch(dismo::maxent(cov_df, 
                                                                         p= occ,
                                                                         factors = factors.new, 
                                                                         args= c("replicatetype=crossvalidate", 
                                                                                 "replicates=4", 
                                                                                 "-J","-P", "-p", "-h", 
                                                                                 "threshold=FALSE")), 
                                                           error = function(e) NA)
                                       
                                       if(class(mod.new) != "MaxEntReplicates"){
                                         cat(paste0("Model 2 failed to fit ... \n"), file = logfile, append = TRUE)
                                         return(NA)
                                       } else {
                                         
                                         modresults <- mod.new@results
                                         
                                         model_2_end <- Sys.time()
                                         cat(paste0("Cross-validation model end (model 2) = ", model_2_end, "\n\n"), file = logfile, append = TRUE)
                                         rm(mod, mod.new)
                                         gc()
                                         
                                         
                                         ## Model 3: Final model ####
                                         model_3_start <- Sys.time()
                                         cat(paste0("Final model start (model 3) = ", model_3_start, "\n"), file = logfile, append = TRUE)
                                         
                                         out_spec <- file.path(output, species_list[i])
                                         mod.final <- 
                                           tryCatch(dismo::maxent(cov_df, p= occ,
                                                                  factors = factors.new,
                                                                  args= c("randomtestpoints=25",
                                                                          "-J","-P", "-p", "-h", 
                                                                          "threshold=FALSE", "writeplotdata=TRUE"), 
                                                                  path = out_spec), 
                                                    error = function(e) NA)
                                         
                                         model_3_end <- Sys.time()
                                         cat(paste0("Final model end (model 3) = ", model_3_end, "\n\n"), file = logfile, append = TRUE)
                                         
                                         if(class(mod.final) != "MaxEnt"){
                                           cat(paste0("Model 3 failed to fit ... \n"), file = logfile, append = TRUE)
                                           return(NA)
                                         } else {
                                           if(predictions == FALSE) { next }
                                           
                                           
                                           ## Prediction ####
                                           prediction_start <- Sys.time()
                                           cat(paste0("Start prediction .... ", "\n"), file = logfile, append = TRUE)
                                           cat(paste0("Prediction start = ", prediction_start, "\n"), file = logfile, append = TRUE)
                                           
                                           ## Get maxSSS threshold to deliniate suitable habitat (Liu et al 2016) 
                                           thresh_maxSSS <- mod.final@results[names(mod.final@results[,1])%in%"Maximum.test.sensitivity.plus.specificity.Cloglog.threshold"]
                                           ## Note: Previously 'Maximum.test.sensitivity.plus.specificity.logistic.threshold' which is not listed anymore...???
                                           
                                           ## Get present cov set and static variables (don't change) for predictions
                                           cov <- stack(list.files(file.path(temp, "pres"), full.names = TRUE))
                                           sta <- cov[[-which(grepl(paste0(c("landuse", "bio"), 
                                                                           collapse = "|"), names(cov)))]]
                                           area_maxSSS <- numeric()
                                           cells_maxSSS <- list()
                                           m <- 0 # index counter for cells_maxSSS
                                           
                                           for (k in 1:length(treatments)) {
                                             
                                             for (j in 1:length(scens)) { #scens are a combination of SSP * Quartile
                                               m <- m + 1
                                               ##Load data to predict under different treatments and scenarios
                                               ## Present prediction (current land use and climate)
                                               if (treatments[k] == "pre") {
                                                 dyn <- cov[[which(grepl("bio", names(cov)))]]
                                                 luc <- cov[[which(grepl("landuse", names(cov)))]]
                                                 
                                                 ##'Aggregated' Land use predcitors as per ssp and climate predictors as per rcps
                                               } else if(treatments[k] == "agg") {
                                                 dyn <- stack(list.files(file.path(temp, 
                                                                                   paste0("bio", scens_rcps[j], "_", 
                                                                                          region_preds)), full.names = TRUE))
                                                 luc <- stack(stack(), list.files(file.path(temp, 
                                                                                            paste0("predlu_", scens[j], "_", 
                                                                                                   region_preds)), full.names = TRUE))
                                               }
                                               
                                               #   ## Direct: Land use stays the same (current), only clime predictors change as per rcps
                                               #  else if(treatments[k] == "dir") {
                                               #   dyn <- stack(list.files(file.path(temp, 
                                               #                                     paste0("bio", scens_rcps[j], "_", 
                                               #                                            region_preds)), full.names = TRUE))
                                               #   luc <- cov[[which(grepl("landuse", names(cov)))]]
                                               #   
                                               #   ## Indirect: Land use changes as per ssp, rcp predictors stay the same (current)
                                               # } else if(treatments[k] == "ind") {
                                               #   dyn <- cov[[which(grepl("bio", names(cov)))]]
                                               #   luc <- stack(stack(), list.files(file.path(temp, 
                                               #                                     paste0("predlu_", scens[j], "_", 
                                               #                                            region_preds)), full.names = TRUE))
                                               # }
                                               
                                               ## Combine static, dynamic and landuse covariates
                                               cov.new <- stack(sta, dyn, luc)
                                               cov.new <- stack(cov.new[[which(names(cov.new)%in%names.new)]])
                                               
                                               ## Predicted map as per region
                                               map <- dismo::predict(mod.final, cov.new)
                                               rm(cov.new)
                                               
                                               ## Apply bioregion constraint for aus (see Methods)
                                               if (region == "aus"){
                                                 map[which(!bio_region[]%in%na.omit(unique(bio_region[inds])))] <- NA
                                               }
                                               
                                               ## Get number of cells above or equal to maxSSS threshold: This is deemed suitable habitat
                                               area_maxSSS <- c(area_maxSSS, length(which(map[] >= thresh_maxSSS)))
                                               cells_maxSSS[[m]] <- which(map[] >= thresh_maxSSS)
                                               
                                               if(treatments[k] == "pre") { break }
                                             }
                                           }
                                           
                                           prediction_end <- Sys.time()
                                           cat(paste0("Prediction end = ", prediction_end, "\n\n"), file = logfile, append = TRUE)
                                         }
                                       }
                                     }
                                     
                                     rm(dyn, luc, mod.final)
                                     gc()
                                     
                                     job_species_end <- Sys.time()
                                     cat(paste0(">> Species job end = ", job_species_end, "\n\n"), file = logfile, append = TRUE)
                                     cat("#########################################  \n\n", file = logfile, append = TRUE)
                                     
                                     meta_data <- list(job_duration = job_species_end - job_species_start,
                                                       model_1_duration = model_1_end - model_1_start,
                                                       model_2_duration = model_2_end - model_2_start,
                                                       model_3_duration = model_3_end - model_3_start,
                                                       prediction_duration = prediction_end - prediction_start)
                                     
                                     list(modresults, area_maxSSS, i, species_list[i], cells_maxSSS, meta_data)
                                   }
  
  ## Store results
  saveRDS(results, file = file.path(output_path, paste0("results_", region, ".rds")))
  ## Structure of results
  ## results is a list of length(species)
  ## results[[x]][[1]] = maxent model output
  ## results[[x]][[2]] = vector with areas_maxSSS values (i.e. total habitat area); length = 1 pres treatment + ((length(treatments) - 1) * length(scens)) = 10
  ## results[[x]][[3]] = index of species
  ## results[[x]][[4]] = name of species
  ## results[[x]][[5]] = list of cells_maxSSS vectors (i.e. cells with suitable habitat); length = 1 pres treatment + ((length(treatments) - 1) * length(scens)) = 10
  ## results[[x]][[6]] = meta_data
  
  # unlink(temp, recursive = TRUE)
  stopCluster(cl)
  
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
  
  rm(list = ls())
  gc()
}


## Check files
spfiles <- list.files("./output/models_aus", pattern = ".txt")
spfiles <- gsubfn::strapplyc(spfiles, "(.*)_log", simplify = TRUE)
spdirs <- list.dirs("./output/models_aus", recursive = FALSE)
spdirs <- gsubfn::strapplyc(spdirs, "_aus/(.*)", simplify = TRUE)

'%!in%' <- function(x,y)!('%in%'(x,y))
splost_idx <- which(spfiles %!in% spdirs)
splost_names <- spfiles[splost_idx]
species_list[splost_idx]


## Species for which only one covariate is retained for the model
## show error: 
## Field Training gain with only bio16 not found in /tempdata/tmp/RtmpMIwsLt/raster/maxent/6793395443/maxentResults.csv
## Warning: Error calculating replicate summary: file /tempdata/tmp/RtmpMIwsLt/raster/maxent/6793395443/maxentResults.csv missing values for field "Training gain with only bio16"

# > species_list[splost_idx]
# [1] "Amaurornis cinerea"    "Anthus australis"      "Rhipidura rufiventris" "Tadorna radjah"       
# > splost_idx
# [1]  34  59 550 591

