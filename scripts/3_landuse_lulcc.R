## on mrc

## Landuse model - lulcc package, Moulds et al. 2015
## Using USGS - MODIS data: https://archive.usgs.gov/archive/sites/landcover.usgs.gov/global_climatology.html 
## 
## Outputs: 
## 1. lucovs - raster stack - covariates for lulcc model [2 files for vn, aus]
## 2. lulcc_dat - workspace - covariates and variables [2 files for vn, aus]
## 3. predlu - raster layer - predicted land use 
##      [18 files: vn, aus for 0.25, 0.5 and 0.75 quartiles 
##      under ssp 1, 2, and 3 + 2 text log files]
## 4. finaldmd - matrix - predicted demand - 6 years x 6 lu classes 
##      [18 files: vn, aus for 0.25, 0.5 and 0.75 quartiles 
##      under ssp 1, 2, and 3]
##
## TOTAL FILES CREATED IN OUTPOT FOLDER (output): 46


## Set up work environment ####
setwd("./regSSP_birds/")

rm(list = ls())
gc()
# system("ps")
# system("pkill -f R")

x <- c('sp', 'raster', 'glmnet', 'lulcc', 'doParallel')
lapply(x, require, character.only = TRUE)
# library(lattice)
# library(latticeExtra)
# library(RColorBrewer)
# library(rasterVis)
# library(Matrix)
# library(proto)
# library(gsubfn)

regSSP_data <- "../regSSP_data" # '/Volumes/discovery_data/regSSP_data' 
rdata_path <- "./RData" # file.path(regSSP_data, "RData_lulcc") 
output_path <- "./output" # file.path(regSSP_data, "output_lulcc") 

source(file.path(".", "scripts", "0_functions.R"))
files <- list.files(rdata_path, full.names = TRUE)

## PREPARE LANDUSE MODEL DATA ####
## ***** loop through each region here ***** ####
regions <- c('vn', 'aus')
for(region in regions){
  
  ## Log file
  job_start <- Sys.time()
  masterlog <- paste0("./landuse_data_run_", region, ".txt")
  writeLines(c(""), masterlog)
  cat(paste0(">> Job start = ", job_start, "\n"), file = masterlog, append = TRUE)
  
  print(paste0("processing ", region, "... "))
  
  ## Load data ####
  reg_mask <- readRDS(file.path(rdata_path, paste0("mask_", region, ".rds"))) 
  covariates_all <- readRDS(files[grepl(paste0("covariates_", region), files)])
  print("original covariate set: ")
  print(names(covariates_all))
  
  ## Landuse layer
  lu_original <- lu_reduced <- covariates_all$landuse
  print("original landuse classes: ")
  print(table(lu_original[]))
  
  print("Check: 0 values in land use layer: ")
  print(lu_original[lu_original == 0]) # check for zero values
  
  lu_classes <- c("urban", "crop", "grass", "shrub", "Oforest", "Cforest")
  lu_exclude <- c(7,8) 
  ## lu 7 & 8 are not modelled but added back to the predicted maps
  ## lu_7: 'herbaceous wetland, moss, lichen' 
  ## lu_8: 'bare, sparese, ice'
  lu_reduced[lu_reduced%in%lu_exclude] <- NA
  print("reduced landuse classes: ")
  print(table(lu_reduced[]))
  
  ## Protected areas layer
  pa <- covariates_all$pa
  print("original protected area layer: ")
  print(table(pa[])); plot(pa)
  ## mask PA classes that can't change
  pa[which(lu_original[]%in%lu_exclude)] <- 0
  print("reduced protected area layer: ")
  print(table(pa[])); plot(pa)
  
  
  ## Define observed landuse at t=0 ####
  obs <- ObsLulcRasterStack(x=lu_reduced,categories=c(1:6),
                            labels=paste0("lu", 1:6), t = 0)
  
  
  ## Prepare explanatory variables ####
  covs_exclude <- c("pa", "popd", "landuse")
  covariates <- covariates_all[[-which(grepl(paste(covs_exclude, collapse = "|"), 
                                             names(covariates_all)))]] 
  print("reduced covariate set: ")
  print(names(covariates))
  
  ## Discard correlated covariates
  cov_df <- getValues(covariates)
  print(paste0("# data points in covariates: ", nrow(cov_df)))
  cov_df <- na.omit(cov_df)
  print(paste0("# data points in covariates without NAs: ", nrow(cov_df)))
  
  cov_df <- scale(cov_df) # scale values for comparison
  preds <- colnames(correlations(cov_df)) # see functions 
  covariates <- covariates[[which(names(covariates)%in%preds)]]
  cov_names <- names(covariates)
  bio_inds <- which(grepl(cov_names, pattern = "bio"))
  
  ## Centre and scale covariate data for lulcc
  cov_cent <- attr(cov_df,c("scaled:center"))[which(colnames(cov_df)%in%cov_names)]
  cov_scal <- attr(cov_df,c("scaled:scale"))[which(colnames(cov_df)%in%cov_names)]
  covariates <- stack(scale(covariates, center = cov_cent, scale = cov_scal))
  ## Rename covariates as per lulcc
  names(covariates) <- paste0("ef_", sprintf("%02d", 1:(nlayers(covariates))))
  
  
  ## Models ####
  print(paste0("defining models for ", region, " ..."))
  
  ## Partition observed data into training and test set for GLM
  train.dat.partition <- 20000 # can change to 10% later
  part <- sample(which(!is.na(obs[[1]][])), train.dat.partition)
  
  ## Extract covariate data for the training partition
  ef <- ExpVarRasterList(covariates)
  train.data <- getPredictiveModelInputData(obs=obs, ef=ef, cells=part, t=0)
  train.data <- as.matrix(train.data)
  
  ## Define models for each lu class 
  lu_cats <- obs@labels
  forms <- list()
  
  for(i in 1:length(lu_cats)){
    
    ## Get coefficients by fitting glmnet on training data 
    test <- cv.glmnet(y = train.data[,i], 
                      x = train.data[,-c(1:6)], family = "binomial", alpha = 1)
    
    ## Identify non-significant covariates
    not_sig <- names(which(coef(test)[,1] == 0))
    
    ## Specify model formula for each lu class using significant variables only
    if(length(not_sig)==0){
      forms[[i]] <- 
        as.formula(paste(lu_cats[i], "~", paste(names(ef), collapse="+")))
    }else{
      forms[[i]] <- 
        as.formula(paste(lu_cats[i], "~", 
                         paste(names(ef)[-which(names(ef)%in%not_sig)], 
                               collapse="+")))
    }
    print(paste0(lu_cats[i], "..."))
  }
  
  ## Fit models (glms) on training data
  train.data <- as.data.frame(train.data)
  glm.models <- glmModels(formula=forms, family=binomial, 
                          data=train.data, obs=obs)
  
  # see link for error: https://stats.stackexchange.com/questions/11109/how-to-deal-with-perfect-separation-in-logistic-regression
  
  ## Save variables for region ####
  print(paste0("saving lulcc data for ", region, " ..."))
  save(
    lu_original,
    lu_classes,
    lu_exclude,
    pa,
    obs,
    preds,
    cov_names,
    bio_inds,
    cov_cent,
    cov_scal,
    ef,
    glm.models,
    file = paste0(output_path, "/lulcc_dat_", region, ".RData")
  )
  saveRDS(readAll(covariates), file = paste0(output_path, "/lucovs_", region, ".rds"))
  print(paste0("lulcc data for ", region, " saved as: ", 
               paste0(output_path, "/lulcc_dat_", region, ".RData")))
  print(paste0("lulcc covariates for ", region, " saved as: ", 
               paste0(output_path, "/lucovs_", region, ".rds")))
  print("-----------------------------------")
  
  rm(list=setdiff(ls(), c("regSSP_data",
                          "rdata_path",
                          "output_path",
                          lsf.str(pattern = "correlations"))))
  gc()
  
  ## Log job duration
  masterlog <- paste0("./landuse_data_run_", region, ".txt")
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

rm(list = ls())
gc()


## RUN LANDUSE MODEL ####
## ***** loop through each region here ***** ####
regions <- c("vn", "aus") # models run for tile 29 which encompasses Vietnam, and Australia
for (region in regions) {
  
  # system("ps")
  # system("pkill -f R")
  .rs.restartR()
  
  ## Set up work environment ####  
  setwd("./regSSP_birds/")
  x <- c('sp', 'raster', 'glmnet', 'lulcc', 'doParallel', 'foreach')
  lapply(x, require, character.only = TRUE)
  
  regSSP_data <- "../regSSP_data" # '/Volumes/discovery_data/regSSP_data' 
  rdata_path <- "./RData" # file.path(regSSP_data, "RData_lulcc") 
  output_path <- "./output" # file.path(regSSP_data, "output_lulcc") 
  
  source(file.path(".", "scripts", "0_functions.R"))
  files <- list.files(rdata_path, full.names = TRUE)
  
  ## Log file
  job_start <- Sys.time()
  masterlog <- paste0("./lulcc_model_run_", region, ".txt")
  writeLines(c(""), masterlog)
  cat(paste0(">> Job start = ", job_start, "\n"), file = masterlog, append = TRUE)
  
  
  ## Predictions ####
  ## Load workspace for region
  load(paste0(output_path,"/lulcc_dat_", region, ".RData"))
  covariates <- readRDS(paste0(output_path, "/lucovs_", region, ".rds"))
  reg_mask <- readRDS(file.path(rdata_path, paste0("mask_", region, ".rds"))) 
  
  ## Define clues parameters ####
  w <- matrix(1, 6, 6) 
  w[1,2:6] <- 0 # because lu_1 'urban' cannot change to any other land use
  elas = c(1, 0.6, 0.5, 0.5, 0.6, 0.9)
  clues.parms <- list(jitter.f=0.002,
                      scale.f=0.000001,
                      max.iter=5000,
                      max.diff= 200,
                      ave.diff= 200)
  
  # elas = c(1, 0.6, 0.6, 0.6, 0.7, 0.9)
  # clues.parms <- list(jitter.f = 0.002,
  #                     scale.f = 0.000001,
  #                     max.iter = 5000,
  #                     max.diff = 500,
  #                     ave.diff = 500)
  
  ## Notes: 
  ## Model does not converge WITH warning at max.iter + warning = 5000/10000/20000
  ## Model does not converge WITHOUT warning at max.iter = 20000, elas = c(1, 0.6, 0.7, 0.7, 0.7, 0.9)
  
  
  ## Specify scenarios ####
  ssps <- paste0("ssp", 1:3)
  rcps <- c("45", "60", "85")
  
  ## Quartiles under which to predict (main results are under q2, the medians)
  quartiles <- c("q2", "q1", "q3")
  
  ## Specify time steps
  yrs <- 2070 - 2019
  timesteps <- c(11, 21, 31, 41, 51)
  
  ## Initiate log file
  obs_n <- length(which(!is.na(obs@layers[[1]][]))) # non-NA values in observed lu data
  logfile <- file.path(output_path, paste0("landuse_", region, "_log.txt"))
  writeLines(c(""), logfile)
  
  ## Run Model ####
  ## Prepare final demand trajectories (including grass/shrub and forest estimations)
  cl <- makeCluster(3, outfile="")
  registerDoParallel(cl)
  clusterCall(cl, function() library('sp', 'raster', 'lulcc'))
  
  for (j in 1: length(ssps)){
    
    ## Load observed land-use data
    obs2 <- obs
    
    ## Load land-use demand for ssp-rcp scenario
    dmd <- as.matrix(readRDS(file.path(output_path, paste0("landdemand_", 
                                                           region, "_", ssps[j], 
                                                           ".rds"))))
    
    foreach(k = 1:length(quartiles), .packages = c("sp", "raster", "lulcc")) %dopar% {
      
      msg1 = paste0("Preparing quantile for region ", 
                    region, " ... Quartile: ", quartiles[k], " | Scenario: ", 
                    ssps[j],  "-rcp", rcps[j], "\n")
      msg2 = paste0(Sys.time(), "\n")
      cat(paste(msg1, msg2, sep = '\n'),  file = logfile, append = T)
      print(msg1)
      print(Sys.time())
      
      ## Load dynamic covariates (i.e. future bioclim layers) for ssp-rcp scenario and quartile
      fut_bioclim <-
        readRDS(files[grepl(paste0("bio", paste0(c(
          quartiles[k], rcps[j], region), collapse = "_")), files)])
      
      ## Discard correlated covariates
      fut_bioclim <-
        fut_bioclim[[which(grepl(paste0(preds, collapse = "|"), names(fut_bioclim)))]]
      
      ## Name covariates as per lulcc
      names_dyn_cov <- names(covariates)[bio_inds]
      names(fut_bioclim) <- names_dyn_cov
      
      ## Crop raster by mask
      fut_bioclim <- crop(fut_bioclim, reg_mask)
      
      ## Scale covariates
      fut_bioclim <-
        stack(scale(readAll(fut_bioclim), center = cov_cent[bio_inds], 
                    scale = cov_scal[bio_inds]))
      
      ## Estimate annual increment of dynamic covariates (i.e. future bioclim layers)
      annual_incr_dyn <- (fut_bioclim - covariates[[bio_inds]])/yrs
      
      ## Load static and dynamic covariates into separate objects
      cov_stat <- covariates[[-bio_inds]]
      # cov_names[-bio_inds]
      
      cov_dyn <- covariates[[bio_inds]]
      # cov_names[bio_inds]
      rm(fut_bioclim)
      
      
      for(t in 1:length(timesteps)){
        
        msg1 = paste0("Preparing model for region ", 
                      region, " ... Quartile: ", quartiles[k], " | Scenario: ", 
                      ssps[j],  "-rcp", rcps[j], " | Timestep: ", t, "\n")
        msg2 = paste0(Sys.time(), "\n")
        cat(paste(msg1, msg2, sep = '\n'),  file = logfile, append = T)
        print(msg1)
        print(Sys.time())
        
        ## Prepare timestep specific data stack to predict lu suitability
        covariates_final <- stack(cov_stat, cov_dyn + 
                                    (annual_incr_dyn * timesteps[t]))
        gc()
        
        ## Estimate predicted suitability for grass, shrub and Oforest and Cforest lu classes
        subs2 <- sample(which(!is.na(covariates_final[[1]][])), 10000)
        ## Predict for lu_classes[3]: "grass"
        pgrass <- 
          mean(predict(glm.models@models[[c(3)]], 
                       data = covariates_final[subs2], type = "response"))
        ## Predict for lu_classes[4]: "shrub"
        pshrub <- 
          mean(predict(glm.models@models[[c(4)]], 
                       data = covariates_final[subs2], type = "response"))
        ## Predict for lu_classes[5]: "Oforest"
        pforesto <- 
          mean(predict(glm.models@models[[c(5)]], 
                       data = covariates_final[subs2], type = "response"))
        ## Predict for lu_classes[6]: "Cforest"
        pforestc <- 
          mean(predict(glm.models@models[[c(6)]], 
                       data = covariates_final[subs2], type = "response"))
        
        ## Proportional suitability for each lu class
        sp <- c(pgrass, pshrub, pforesto, pforestc)
        prop_suit <- sp / sum(sp)
        rm(sp)
        
        ## Area avialable for change: total area - (ar(urban) + ar(crop))
        area_avail <- obs_n - sum(dmd[t+1,], na.rm = TRUE)
        
        ## Estimate lu demand based on predicted suitability * area available
        dmd[t+1,3] <- round(area_avail * prop_suit[1])
        dmd[t+1,4] <- round(area_avail * prop_suit[2])
        dmd[t+1,5] <- round(area_avail * prop_suit[3])
        dmd[t+1,6] <- round(area_avail * prop_suit[4])
        
        ## Estimate unallocated land area at t+1 compared to t
        unalloc_land <- rowSums(dmd)[1] - rowSums(dmd)[t + 1]
        
        ## Allocate unallocated land area to a randomly selected landuse class in t+1
        col_adj <- sample(1:ncol(dmd), 1)
        dmd[t+1, col_adj] <- dmd[t+1, col_adj] + unalloc_land
        
        ## Check if adjusted allocation sums to total area
        if(rowSums(dmd)[1] - rowSums(dmd)[t + 1] != 0) stop("sum of demand not equal to total area")
        
        ## Reassign original predictor names so they can be found by predict function
        names(covariates_final)[bio_inds] <- names_dyn_cov
        ef <- ExpVarRasterList(readAll(covariates_final), pattern = "ef")
        ## gives warning: 
        # Warning message:
        #   In fun(libname, pkgname) : couldn't connect to display ":0"
        
        ## Sync NAs otherwise clues won't converge...
        inds <- unique(c(which(is.na(obs2@layers[[1]][])), 
                         which(is.na(ef@maps[[1]][]))))
        
        for (efs in 1:length(ef)){
          ef@maps[[efs]][inds] <- NA
        }
        pa[inds] <- NA
        obs2@layers[[1]][inds] <- NA
        dmd[t,] <- table(obs2[])
        
        msg1 = paste0("Starting CLUE for region ", 
                      region, " ... Quartile: ", quartiles[k], " | Scenario: ", 
                      ssps[j],  "-rcp", rcps[j], " | Timestep: ", t, "\n")
        msg2 = paste0(Sys.time(), "\n")
        cat(paste(msg1, msg2, sep = '\n'),  file = logfile, append = T)
        print(msg1)
        print(Sys.time())
        
        ## Specify clues model
        clues.model <- CluesModel(obs = obs2, 
                                  ef = ef, 
                                  models = glm.models, 
                                  time = 0:1,
                                  demand=dmd[c(t:(t+1)),], 
                                  rules = w, 
                                  params = clues.parms, 
                                  elas = elas,
                                  mask = pa)
        
        ## Allocate predicted change in lu classes
        st <- system.time(ordered_out <- allocate(clues.model))
        cat(paste(st, "\n"), file = logfile, append = T)
        
        ## Compare input and output
        print("Compare input (demand) and output (allocated demand) from lulcc: \n")
        dmd
        summary(ordered_out)$output
        
        # ## Allocate changes in line with demand, suitability model, and allocation order (see Methods)
        # ## Make model object
        # ordered.model <- OrderedModel(obs=obs2,
        #                               ef=ef,
        #                               models=glm.models,
        #                               time=0:1,
        #                               demand=dmd[c(i:(i+1)),],
        #                               order=c(1:6),
        #                               rules = w,
        #                               mask = pa)
        # ## Allocate predicted change in lu classes
        # ordered_out <- allocate(ordered.model, stochastic = F)
        
        ## New land-use map becomes input for next iteration of t
        obs2 <- ObsLulcRasterStack(ordered_out@output[[2]])
        print(paste0("Quartile: ", quartiles[k], " | Scenario: ", 
                     ssps[j],  "-rcp", rcps[j], " | Timestep: ", t))
      }
      rm(ef)
      
      ## Save outputs
      out <- ordered_out@output[[2]]
      
      ## Check if total area in output == total area in observed
      print("Is total area of predicted lu map = total area of observed lu map: ")
      sum(table(out[])) == obs_n
      
      ## Add lu classes 7 & 8 to predicte lu
      out[which(lu_original[]%in%lu_exclude[1])] <- 7
      out[which(lu_original[]%in%lu_exclude[2])] <- 8
      names(out) <- "landuse"
      
      saveRDS(readAll(out), file = paste0(output_path, "/predlu_", quartiles[k], 
                                          "_", ssps[j], "_", region,  ".rds"))
      
      saveRDS(dmd, file = paste0(output_path, "/finaldmd_", quartiles[k], 
                                 "_", ssps[j], "_", region,  ".rds"))
    }
  }
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
}

rm(list=ls())
gc()



## Check output: Compare output of Clues model (finaldmd) with projected land damand from GTAP (predlu)
region <- "aus"
x <- list.files(output_path, pattern = region,  full.names = TRUE)
xr <- x[grepl("predlu", x)]
xd <- x[grepl("finaldmd", x)]
landuse <- readRDS(paste0("./RData/covariates_", region, ".rds"))
landuse <- landuse$landuse
lu_tabled <- table(landuse[])
rs <- list()
print("Land use classes considered: \n")
print(lu_classes)
for(i in 1:(length(xr))){
  r <- readRDS(xr[i])
  rs[[i]] <- r
  d <- readRDS(xd[i])
  check <- round(table(r[])[1:6]/d[6,], 0) #compare values
  print(paste0("Is CLUEs output = projected GTAP demand for: ", 
               stringr::str_extract(xr[i], "q.{0,10}")))
  print(all(check == 1))
  print("Final predicted demand (CLUEs model):")
  print(d[6,])
}

## ?? ...
# rs_df <- as.data.frame(stack(rs))
# for(j in 1:6){
#   for(i in 1:6){
#     print(identical(rs_df[,j], rs_df[,i]))
#   }
# }
plot(obs2@layers[[1]])


## LU model not converging....
## Warning implies that the model did not converge. To check compare 
## table(ordered_out@obs[])
## with dmd in t+1
## when run in parallel, error: error in unserialize(socklist[[n]]) : error reading from connection
## seems to be a problem in parallelisation, run in serial or ask for more memory on boab/run on MRC
