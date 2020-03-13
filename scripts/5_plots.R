## Set work environment ####
setwd("./regSSP_birds/")
x <- c('sp', 'raster', 'RColorBrewer', 'data.table', 'readxl', 'lattice', 'latticeExtra', 'rasterVis', 'pBrackets')
lapply(x, require, character.only = TRUE)

regSSP_data <- "../regSSP_data" # '/Volumes/discovery_data/regSSP_data'
rdata_path <- "./RData" # file.path(regSSP_data, "RData_lulcc")
output_path <- "./output" # file.path(regSSP_data, "output_lulcc") 
fig_path <- "./figures"
if(!dir.exists(fig_path)){dir.create(fig_path)}

## Global parameters ####
ssps <- paste0("ssp", 1:3)
rcps <- c("45", "60", "85")
quartiles <- c("q2", "q1", "q3")
scens <- sort(apply(expand.grid(quartiles, ssps), 1, paste0, collapse="_"))


## COMMODITY & ASSOCIATED LANDUSE CHANGE PLOTS ####

## >> Specify region ####
regions <- c('vn', 'aus')
region = 'aus'

## Barplots for gtap_endowments ####
dmd1 <- readRDS(paste0(rdata_path, "/gtap_landendowments_", region, "_ssp1.rds"))
dmd2 <- readRDS(paste0(rdata_path, "/gtap_landendowments_", region, "_ssp2.rds"))
dmd3 <- readRDS(paste0(rdata_path, "/gtap_landendowments_", region, "_ssp3.rds"))
dmd <- cbind(dmd1[,52], dmd2[,52], dmd3[,52]) # demand in 2070
dmd <- (dmd - mean(dmd))/(max(dmd) - min(dmd)) #mean normalization
dmd <- dmd*100 # convert to percentages
colnames(dmd) <- c("ssp1","ssp2","ssp3")
rownames(dmd) <- c("paddy", "wheat","cereal grains", "vegetables,fruits,nuts", "oil seeds", "sugarcane, sugarbeet", "plant-based fibres", "other crops", "cattle,sheep,goats,horses", "forestry")
cols <- brewer.pal(dim(dmd)[1], "RdBu")
par(mar = c(5,4,9,6), oma = c(0,0,0,0))
barplot(dmd, horiz = T, beside = T, col = cols, xlab = "% change in land harvested\n (2019-2070)", ylab = "scenarios")
par(xpd=TRUE)
legend(20,47, rownames(dmd), fill = cols)

## Plot and legend separate
par(mar = c(5,5,2,2), oma = c(0,0,0,0))
barplot(dmd, horiz = T, beside = T, col = cols, xlab = "% change in land harvested (2019-2070)", ylab = "scenarios", cex.axis = 1.5, cex.names = 1.5, cex.lab = 1.5)
plot(NULL)
legend(0.2,1, rownames(dmd), fill = cols, cex = 1.2)

## SSP 3 (dmd3) only
par(mar = c(5,15,1,1), oma = c(0,0,0,0))
x <- barplot(dmd[,1], horiz = T, beside = T, col = cols, yaxt="n", xlab = "% change in land harvested (2019-2070)", cex.axis = 1.2, cex.names = 1.2, cex.lab = 1.2)
labs <- rownames(dmd)
text(x=-28, y=x, labels = labs, xpd=TRUE, srt=0, cex=1.2, pos=2)


## GTAP data plot ####
ssp = "ssp3"
temp <- as.data.table(read_xls(file.path(regSSP_data, "gtap_data", "ssps",
                                         paste0(region, "_", ssp, ".xls")), sheet = "land"))
temp <- temp[X__1 %in% rownames(dmd1)]
temp[, c(2:6, ncol(temp)) := NULL]
temp[, c(20:29) := NULL] ## retain data upto 2070 only
temp[,1 := NULL]
temp<- as.matrix(temp)
rownames(temp) <- rownames(dmd1)
temp2070 <- as.matrix(temp[,c(18)])
rownames(temp2070) <- rownames(dmd1)
par(mar = c(5,4,9,6), oma = c(0,0,0,0))
barplot(temp2070, horiz = T, beside = T, col = cols, xlab = "% change in land harvested\n (2019-2070)", ylab = "scenarios")
par(xpd=TRUE)
legend(20,47, rownames(temp2070), fill = cols)


## Total land use change table for vn and aus ####
luchange <- data.frame()
lu_classes <- c("Urban", "Cropland", "Herbaceous", "Shrubs", "Open forest", "Closed forest")
scensq2 <- scens[grep("q2_", scens)]
ctr <- 1

for (region in regions){
  m <- matrix(NA, ncol = 4, nrow = 6)
  ## Pull out data for year0 from any scenario
  m[, 1] <- as.numeric(readRDS(paste0(output_path, "/finaldmd_", scensq2[1], "_", region, ".rds"))[1,])
  for (i in 1:3){
    ## Pull out data for last year for each scenario
    m[, i + 1] <- as.numeric(readRDS(paste0(output_path, "/finaldmd_", scensq2[i], "_", region, ".rds"))[6,])
  }
  colnames(m) <- paste0(c("yr0_", "ssp1_", "ssp2_", "ssp3_"), region)
  if (ctr == 1) {
    luchange <- rbind(luchange, as.data.frame(m))
  } else {
    luchange <- cbind(luchange, as.data.frame(m))
  }
  ctr <- ctr +1
}

# ## Calculate change relative to entire landscape
# rs <- colSums(luchange) .... ???
# mn <- cbind(round((luchange[,1:4] - luchange$yr0_vn), 2), round((luchange[,5:8] - luchange$yr0_aus), 2))
# mn <- mn[, -(grep("yr0", colnames(mn)))]
# 
# luchange <- round(t(t(mn)/rs) * 100, 2)
# cols1 <- colorRampPalette(c("darkred", "white"))(50)
# cols2 <- colorRampPalette(c("white", "darkgreen"))(50)
# breaks <- seq(-3, 3, length.out = 101)
# new_breaks <- scales::rescale(breaks^3, to = c(-4, 4))
# colnames(luchange) <- rep(c("SSP-1", "SSP-2","SSP-3"), 2)
# rownames(luchange) <- lu_classes
# 
# ## Plot
# x <- 1:ncol(luchange)
# y <- 1:nrow(luchange)
# centers <- expand.grid(y,x)
# par(mar = c(2, 8, 5, 2))
# image(x, y, t(luchange),
#       col = c(cols1, cols2),
#       breaks = new_breaks,
#       xaxt = 'n',
#       yaxt = 'n',
#       xlab = '',
#       ylab = '',
#       ylim = c(max(y) + 0.5, min(y) - 0.5)
# )
# 
# text(centers[,2], centers[,1], c(luchange), col= "black", cex = 1.5)
# mtext(paste(attributes(luchange)$dimnames[[2]]), at=1:ncol(luchange), padj = -1, cex = 1.5)
# mtext(c("Australia", "Vietnam") , at= c(1.5, 3.5), padj = -3, cex = 1.5)
# mtext(attributes(luchange)$dimnames[[1]], at=1:nrow(luchange), side = 2, las = 1, adj = 1.05, cex = 1.5)


## PREDICTED LANDUSE MAPS ####
## >> Specify region ####
region  <- 'aus'

## Current land use ####
landuse_all <- landuse <- readRDS(paste0(rdata_path, "/covariates_", region, ".rds"))$landuse
lu_classes <- c("urban", "crop", "grass", "shrub", "Oforest", "Cforest", "wetland", "bare")
# lu_exclude <- c(7,8) 
#   ## lu 7 & 8 are not modelled but added back to the predicted maps
#   ## lu_7: 'herbaceous wetland, moss, lichen' 
#   ## lu_8: 'bare, sparese, ice'
# landuse[landuse%in%lu_exclude] <- NA

par(mfrow=c(1,1), mar = c(0.1, 0.1, 0.1, 0.1))
plot(landuse_all, legend = F, axes=F, box=F, col = ochRe::ochre_palettes$namatjira_qual)
plot(NULL)
legend("center", legend = lu_classes, fill = ochRe::ochre_palettes$namatjira_qual, bty = "n", cex = 1.7)

## Predicted landuse ####
## ssp1
lu11 <- readRDS(paste0(output_path, "/predlu_q1_ssp1_", region, ".rds"))
lu21 <- readRDS(paste0(output_path, "/predlu_q2_ssp1_", region, ".rds"))
lu31 <- readRDS(paste0(output_path, "/predlu_q3_ssp1_", region, ".rds"))
## ssp2
lu12 <- readRDS(paste0(output_path, "/predlu_q1_ssp2_", region, ".rds"))
lu22 <- readRDS(paste0(output_path, "/predlu_q2_ssp2_", region, ".rds"))
lu32 <- readRDS(paste0(output_path, "/predlu_q3_ssp2_", region, ".rds"))
## ssp3
lu13 <- readRDS(paste0(output_path, "/predlu_q1_ssp3_", region, ".rds"))
lu23 <- readRDS(paste0(output_path, "/predlu_q2_ssp3_", region, ".rds"))
lu33 <- readRDS(paste0(output_path, "/predlu_q3_ssp3_", region, ".rds"))

plot(stack(lu11,lu21,lu31,lu12,lu22,lu32,lu13,lu23,lu33), legend = F, axes=F, box=F, col = ochRe::ochre_palettes$namatjira_qual)

plot(lu21, legend = F, axes=F, box=F, col = ochRe::ochre_palettes$namatjira_qual)
plot(lu22, legend = F, axes=F, box=F, col = ochRe::ochre_palettes$namatjira_qual)
plot(lu23, legend = F, axes=F, box=F, col = ochRe::ochre_palettes$namatjira_qual)

plot(NULL)
legend("center", legend = lu_classes, fill = ochRe::ochre_palettes$namatjira_qual, bty = "n", cex = 1.7)


## Plot difference in land use (current - scenario) ####
par(mfrow=c(1,1), mar = c(0.1, 0.1, 0.1, 0.1), oma = c(0, 0, 0, 0))
plot(overlay(landuse, lu21, fun=function(x,y) as.logical(x==y)), 
     col= c("salmon", "steelblue"), box=FALSE, axes=FALSE, legend = F)
plot(overlay(landuse, lu22, fun=function(x,y) as.logical(x==y)), 
     col= c("salmon", "steelblue"), box=FALSE, axes=FALSE, legend = F)
plot(overlay(landuse, lu23, fun=function(x,y) as.logical(x==y)), 
     col= c("salmon", "steelblue"), box=FALSE, axes=FALSE, legend = F)
plot(overlay(lu11, lu23, fun=function(x,y) as.logical(x==y)), 
     col= c("salmon", "steelblue"), box=FALSE, axes=FALSE, legend = F)
plot(NULL)
legend("topright", legend = c("change", "no change"), fill=c("salmon", "steelblue"), bty = "n", cex = 3)


## ********* Landuse change (levelplot) maps ####
cols <- viridis::cividis(101, alpha = 1, begin = 0.2, end = 1, direction = 1)
breaks <- seq(0, 1, length.out = length(cols)) #do we need the sqrt of break points here to?

f0 <- readRDS(file.path(rdata_path, paste0("covariates_", region, ".rds")))$landuse
r_stack <- stack()
# for(k in 1:length(quartiles)){ # for quartile, k[1] = q2
k = 1
for(i in 1:length(ssps)){ # for ssp
  name <- paste0("predlu_", quartiles[k], "_", ssps[i], "_", region)
  print(name)
  f1 <- readRDS(file.path(output_path, paste0(name, ".rds")))
  dif <- f0 - f1
  dif[which(dif[] != 0 & !is.na(dif[]))] <- 1
  #r <- aggregate(dif, fact = 5)
  r <- focal(dif, w = matrix(1, 3, 3), fun = function(x) {mean(x, na.rm = TRUE)})
  names(r) <- name
  r_stack <- stack(r_stack, r)
  par(mar = c(0,0,0,0), oma = c(0,0,0,0))
  png(file.path(fig_path, paste0(name, ".png")), bg = "transparent")
  #plot(r, col = cols, legend = FALSE, axes=FALSE, box=FALSE)
  
  levelplot(sqrt(r), col.regions = cols, colorkey = NULL, margin = F, par.settings = list(
    axis.line = list(col = "transparent"),
    strip.background = list(col = 'transparent'),
    strip.border = list(col = 'transparent')),
    scales = list(draw = FALSE),
    xlab = NULL,
    ylab = NULL,
    at= breaks
  )
  dev.off()
}

## Plot legend
plot(NULL)
cols <- viridis::cividis(101, alpha = 1, begin = 0.2, end = 1, direction = 1)
legend_image <- as.raster(matrix(cols, nrow=1))
par(mar = c(0,0,0,0), oma = c(0,0,0,0))
plot(c(0,1),c(0,1), type = "n", axes = FALSE,xlab = '', ylab = '')
text(x=sqrt(seq(0,1,l=5)), y = 0.25, labels = seq(0,1,l=5), cex = 0.8)
rasterImage(legend_image, 0, 0.5, 1,1)



## HABITAT SUITABILITY MAPS ####
regions <- c("til", "aus")
treatments <- c("pres", paste0("agg_", scens))

## Prepare (& save) species data for plotting ####
## List of species to exclude based on AUC values
auc_list <- list()
exclude_list <- list()
for(j in 1:2){ # for regions
  res <- readRDS(file.path(output_path, paste0("results_", regions[j], ".rds")))
  # res <- readRDS("/Volumes/discovery_data/birds_ccimpacts/output/results_til.rds")
  auc_list[[j]] <- numeric()
  for(i in 1:length(res)){
    auc_list[[j]][i] <- res[[i]][[1]][,1][which(names(res[[i]][[1]][,1]) == "Test.AUC")]
  }
  exclude_list[[j]] <- which(auc_list[[j]] < 0.7)
}
names(auc_list) <- names(exclude_list) <- c("vn", "aus")
saveRDS(auc_list, file.path(output_path, "auclist.rds"))
saveRDS(exclude_list, file.path(output_path, "excludelist.rds"))

## Retain results for significant AUC values for habitat change
final_data <- list()
for(j in 1:2){ # for regions
  results <- res <- readRDS(file.path(output_path, paste0("results_", regions[j], ".rds")))
  species <- sapply(res, FUN = function (x) {cbind(x[[4]])})
  res <- t(sapply(res, FUN = function (x) {cbind(x[[2]])}))
  res <- data.frame(species, res)
  res <- res[-exclude_list[[j]],]
  if(length(which(res[,2] == 0)) != 0) {res <- res[-which(res[,2] == 0),]}
  names(res)[-1] <- treatments
  res[,-1] <- log(res[,-1]/res[,2]) # scenario values/ pres values
  res <- res[-which(!is.finite(rowSums(res[,-1]))),]
  final_data[[j]] <- res
}
names(final_data) <- c("vn", "aus")
saveRDS(final_data, file.path(output_path, "final_sdmdat.rds"))


## Decrease in habitat plot ####
final_data <- readRDS(file.path(output_path, "final_sdmdat.rds"))

lims <- c(log(0.1), log(8))
padj <- - 0.05
wth <- 0.35
jit <- 0.3
tf <- 1
of <- 0.1
region_labels <- c("Vietnam", "Australia")

pdf(file.path(fig_path, "habitat_change.pdf"), pointsize = 12)
for (j in 1:2){
  par(mar = c(0.2,4,0.2,0.2))
  par(mgp = c(3, 0.5, 0))
  plot(1, ylim = lims, xaxt = "n", type = "n", bty = "n", ylab = "", yaxt = "n")
  axis(2, at = c(log(0.125), log(0.25), log(0.5), 
                 log(1), log(2), log(4), log(8)), 
       labels = c(expression(paste(over(1,8), " x")),
                  expression(paste(over(1,4), " x")),
                  expression(paste(over(1,2), " x")),
                  expression(paste(1, " x")),
                  expression(paste(2, " x")),
                  expression(paste(4, " x")),
                  expression(paste(8, " x"))),
       las = 2, cex.axis = 0.8)
  
  title(ylab = "decrease", cex.lab = 1,
        line = 2.5, adj = 0.3)
  
  title(ylab = "increase", cex.lab = 1,
        line = 2.5, adj = 0.7)
  
  ## Plot data points with means, by country (j) and treatment * SSP (i)
  dp <- final_data[[j]][,which(grepl("agg_q2", colnames(final_data[[j]])))]
  dp <- as.matrix(dp)
  alpha <- 0.1
  means <- colMeans(dp)
  par(new=TRUE)
  par(mar = c(0,0,0.2,0))
  plot(1:7, ylim = lims, axes = F, type = "n", xlim = c(0.5, 7.5))
  
  for(i in 1:ncol(dp)){
    pts <- dp[,i]
    
    ## Only plot points within certain bounds
    pts <- pts[which(pts > log(0.125) & pts < log(8))]
    points(rep(i+1, length(pts)) + 
             runif(length(pts), -jit, jit), 
           pts, pch = 20, col = grey(0.7), cex = 1)
    
    ## Means
    segments((i+1) - wth, means[i], (i+1) + wth, means[i],
             col= grey(0), border=par("fg"), lty= 1, lwd= 2, xpd=FALSE)
  }
  
  text(x = 2, y = lims[1] - of - 0.5, "SSP scenarios", cex = tf, xpd = NA)
  title(xlab = region_labels[j], cex.lab = 1.5, line = -32, adj = 0.8)
  scen_labels <- paste0("SSP", c(1,2,3))
  for (i in 1:3){
    text(x = i+1, y = log(0.125) - 0.2, scen_labels[i], cex = tf)
  }
}
dev.off()


## Boxplots for AUC ###
auc_list <- readRDS(file.path(output_path, "auclist.rds"))
exclude_list <- readRDS(file.path(output_path, "excludelist.rds"))
par(mar = c(1.5,8,2,6))
f <- boxplot(auc_list[[1]][-exclude_list[[1]]], auc_list[[2]][-exclude_list[[2]]],
             boxlty = 0,
             staplecol = "black",
             col = "grey",
             boxwex = 0.8,
             boxlty = 0, 
             bty= "n",
             axes = F
)

axis(2, tck = -0.025, at = c(0.7, 0.8, 0.9, 1), labels = c(0.7, 0.8, 0.9, 1))
title(xlab = paste0("n = ", length(auc_list[[1]]) - length(exclude_list[[1]])), adj = 0.1, line = 0)
title(xlab = "Vietnam", adj = 0.05, line =- 6)
title(xlab = paste0("n = ", length(auc_list[[2]]) - length(exclude_list[[2]])), adj = 0.9, line = 0)
title(xlab = "Australia", adj = 0.95, line = - 6)
title(ylab = "AUC", line = 2)


## Plot varibale contributions ####
## Estimate variable contributions
contr_list <- list()
contr_num <- list()
for(j in 1:2){
  preds <- c(readRDS(file.path(output_path, paste0("preds_", regions[j], ".rds"))), "landuse")
  res <- readRDS(file.path(output_path, paste0("results_", regions[j], ".rds")))
  contr <- matrix(ncol = length(preds), nrow = length(res))
  colnames(contr) <- preds
  
  for (i in 1:length(res)){
    pred_values <- res[[i]][[1]][which(grepl("contribution", rownames(res[[i]][[1]])))]
    pred_names <- rownames(res[[i]][[1]])[which(grepl("contribution", rownames(res[[i]][[1]])))]
    n <- unlist(strsplit(pred_names, split = ".", fixed = T))
    pred_names <- n[seq(1,length(n), 2)]
    contr[i, na.omit(match(pred_names, colnames(contr)))] <- pred_values
  }
  contr_list[[j]] <- contr
  contr_num[[j]] <- apply(contr, 2, FUN = function (x) {length(which(is.na(x)))})/length(res) * 100
}
saveRDS(contr_list, file.path(output_path, "var_contributions.rds"))


## Barplots
for(j in 1:2){
  means <- colMeans(contr_list[[j]], na.rm = T)
  ns <- apply(contr_list[[j]], 2, FUN = function (x) {round(length(which(!is.na(x)))/length(x) * 100)})
  means <- means/sum(means) * 100
  names <- names(means)
  names[which(names == "dila")] <- "dist lakes"
  names[which(names == "diri")] <- "dist rivers"
  names(means) <- names
  names(ns) <- c("slope", "dist to lakes", "dist to rivers", 
                 "prec wettest month", "prec warm quarter", 
                 "prec cold quarter", "mean diurnal range", 
                 "max temp warm month", "temp annual range", "landuse")
  #breaks <- barplot(sort(means), horiz = T, axes = F, las = 2, border = NA, xlim = c(0, 25))
  par(mar = c(0,16,0,1))
  breaks <- barplot(sort(ns), horiz = T, axes = F, las = 2, border = NA, xlim = c(0, 100), cex.names = 1.4)
  text(sort(ns)-5, breaks, paste0(sort(ns), "%"), adj = 1, cex = 1.5)
  title(xlab = region_labels[j], line = -2, adj = 1, cex.lab = 1.3)
}



## Species plots ####
## >> Specify region ####
regions <- c("til", "aus")
region <- 'aus'

final_data <- readRDS(file.path(output_path, "final_sdmdat.rds"))
res <- readRDS(file.path(output_path, paste0("results_", region, ".rds")))

## Total habitat area per species
spar <- t(sapply(res, "[[", 2))
colnames(spar) <- treatments
rownames(spar) <- sapply(res, "[[", 4)
spar <- spar[,which(grepl("agg_q2_", colnames(spar)))]
spar <- spar[which(rownames(spar) %in% final_data[[c]]$species),]
spar <- na.omit(spar)
spar <- as.data.frame(spar)

## Create summary rasters from stacked SDMs ####
final_maps <- list()
if(region == "til") region <- "vn"
reg_mask <- readRDS(file.path(rdata_path, paste0("mask_", region, ".rds")))
final_maps <- stack(reg_mask, reg_mask, reg_mask, reg_mask)
names(final_maps) <- c("pres", "ssp1", "ssp2", "ssp3")

sp_preds <-  sapply(res, "[", 5)
names(sp_preds) <- sapply(res, "[[", 4)
sp_preds <- sp_preds[which(names(sp_preds)%in%final_data[[c]]$species)]

treatments <- c("pres", paste0("agg_", scens))
treat_idx <- c(1, which(grepl("agg_q2_", treatments)))

for(k in 1:4){
  ## pull our results for specific treatments across all species
  layer <- sapply(sp_preds, "[", treat_idx[k])
  
  for(j in 1:length(sp_preds)){
    if(!names(sp_preds)[j]%in%final_data[[c]]$species) next
    if(length(layer[[j]]) != 0) final_maps[[k]][layer[[j]]] <- final_maps[[k]][layer[[j]]] + 1 # +1 to value in final_maps, for all cells noted to have suitable habitat
    print(j)
  }
  print(paste0("Processing treatment ", treatments[k], " ..."))
}
saveRDS(final_maps, file.path(output_path, "species_maps.rds"))
final_maps_scaled <- final_maps/nrow(final_data[[c]])

for(k in 1:nlayers(final_maps_scaled)){
  print(paste0("Resampling output map ", k))
  final_maps_scaled[[k]] <- focal(final_maps_scaled[[k]], w = matrix(1, 3, 3), fun = function(x) {mean(x, na.rm = TRUE)})
}
names(final_maps_scaled) <- c("pres", "ssp1", "ssp2", "ssp3")
saveRDS(final_maps_scaled, file.path(output_path, "species_maps_resampled.rds"))

## Plot stacked rasters ####
final_maps_scaled <- readRDS(file.path(output_path, "species_maps_resampled.rds"))
sum_max <- max(values(final_maps_scaled), na.rm = TRUE)
sum_min <- min(values(final_maps_scaled), na.rm = TRUE)
par(mar = c(0.2,0,0,0), oma = c(0,0,0,0))
breaks <- seq(0, 0.5, length.out = 20)
cols <- colorRampPalette(c("lightgoldenrod1", "rosybrown2", "royalblue4"))(length(breaks) - 1)
plot(final_maps_scaled, legend = F, axes=F, box=F, col = cols)
plot(NULL)
plot(final_maps_scaled[[1]], legend.only = T, axes=F, box=F, cex = 1.5, col= cols)

## Level plot ####
f <- final_maps
levelplot(sqrt(f[[1]]), col.regions = cols, margin = F, colorkey= NULL, par.settings = list(
  axis.line = list(col = "transparent"),
  strip.background = list(col = 'transparent'), 
  strip.border = list(col = 'transparent')),
  scales = list(draw = FALSE),
  xlab = NULL,
  ylab = NULL,
  at = breaks
)

## Legend
cols <- cividis(101, alpha = 1, begin = 0.2, end = 1, direction = 1)
legend_image <- as.raster(matrix(cols, nrow=1))
par(mar = c(0,0,0,0), oma = c(0,0,0,0))
plot(c(0,1),c(0,1), type = "n", axes = FALSE,xlab = '', ylab = '')
text(x=sqrt(seq(0,1,l=5)), y = 0.25, labels = seq(0,0.02,l=0.005), cex = 3)
rasterImage(legend_image, 0, 0.5, 1,1)

## Plot differences ####
plot(overlay(final_maps_raw[[4]], final_maps_raw[[1]], fun=function(x,y) as.logical(x==y)), 
     col= c("salmon", "steelblue"), box=FALSE, axes=FALSE, legend = F)
cols <- cividis(101, alpha = 1, begin = 0.2, end = 1, direction = 1)
tm <- c('ind', "dir")
region_labels <- c("Vietnam", "Australia")
dev.off()
maximum <- sqrt(max(c(stack(final_maps[[1]])[]), stack(final_maps[[2]])[], na.rm = TRUE))


## Identify species with declines ####
if(region == 'aus') { c <- 2 } else { c <- 1}
dat <- final_data[[c]]
dat <- dat[,which(grepl("agg_q2", colnames(dat)))]
dat <- cbind(final_data[[c]][,1], dat)

apply(dat[,-1], 2, FUN = function(x) {length(which(exp(x) < 0.5))})
# > 50% decline (median predictions)
apply(dat[,-1], 2, FUN = function(x) {length(which(exp(x) < 0.1))})
# > 90% decline (median predictions)

ssp3_90p_sp <- as.character(dat[,1])[which(exp(dat$agg_q2_ssp3) < 0.1)]
ssp3_90p_idx <- which(dat[,1]%in%ssp3_90p_sp) 

## Extract species
sp_preds <-  sapply(res, "[", 5)
names(sp_preds) <- sapply(res, "[[", 4)
sp_preds <- sp_preds[which(names(sp_preds)%in%ssp3_90p_sp)]

treatments <- c("pres", paste0("agg_", scens))
treat_idx <- c(1, which(grepl("agg_q2_", treatments)))

if(region == "til") region <- "vn"
reg_mask <- readRDS(file.path(rdata_path, paste0("mask_", region, ".rds")))
cols <- colorRampPalette(c("lightgoldenrod1", "rosybrown2", "royalblue4"))(length(breaks) - 1)

for(j in 1: length(sp_preds)){
  sp_maps <- list()
  sp_maps <- stack(reg_mask, reg_mask, reg_mask, reg_mask)
  names(sp_maps) <- c("current", "ssp1", "ssp2", "ssp3")
  
  for(k in 1:4){
    layer <- sp_preds[[j]]
    l <- treat_idx[k]
    if(length(layer[[l]]) != 0) sp_maps[[k]][layer[[l]]] <- sp_maps[[k]][layer[[l]]] + 1
  }
  par(mar = c(0,0,0,0), oma = c(1,0,0,0))
  plot(sp_maps, legend = F, axes=F, box=F, col = cols)
  title(names(sp_preds)[j], cex.main = 1, line=-26, font=2)
}

