# -*- coding: utf-8 -*-
# ---
# This script requires a functions script titled "SDHMfunctions"
# This script is designed to create ensemble species distribution models from presence point data and predictor layers
# Written by William Wiskes 
# Last update 8/29/2021
# ---
# -*- coding: utf-8 -*-
# +

#Load Library
library(sp)
library(sf)
library(rgdal)
library(terra)
library(tidyverse)
library(dplyr)

library(tinytex)
set.seed(1234)

# -

#UTM Nad 83 (for point data)
prj.utmN83z12 <- "+proj=utm +zone=12 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
# https://epsg.io/42303  (5070)use this one
prj.aeaN83 <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"  # epsg:42303



dat <-st_read("Heritage_So_Leatherside_Chub_ce.shp")
dat <- sf::st_transform(dat, 
                        crs = prj.aeaN83)

# +

###Import Ancillary Visualization/Reference Data
utah <- st_read("/vsicurl/https://storage.googleapis.com/predictors_public/aquatics/utah.json") #for pretty map
drainages_ut <- st_read("/vsicurl/https://storage.googleapis.com/predictors_public/aquatics/watersheds.json") 


# +

#Import predictor layer
streams_att <- st_read("/vsicurl/https://storage.googleapis.com/predictors_public/aquatics/aquatics_cleaned.json")  
#streams_att <- st_read("data/gis_layers/streams_pred.shp")



# +

utah <- sf::st_transform(utah, 
                         crs = prj.aeaN83)
streams_sf <- sf::st_transform(streams_att, 
                               crs = prj.aeaN83)
drainages_sf <- sf::st_transform(drainages_ut, 
                                 crs = prj.aeaN83)
st_transformdrainages_sf <- sf::st_transform(drainages_ut, 
                                             crs = prj.aeaN83)

# -

####################################################################
####################################################################
#data preparations
####################################################################
# Rename column by name: change "beta" to "two"
names(drainages_sf)[names(drainages_sf)=="HUC_12"] <- "huc12"
names(drainages_sf)
#replace streams_sf with the merged database that includes huc names
streams_sf <- st_join(streams_sf,drainages_sf,by = "huc12")
####################################################################
#End Data merge
####################################################################

# +
#WAY TO TEST EQUALITY FOR PROJ4 IN FUTURE WILL PUT IN EPSG
##all.equal(st_crs(drainages_ut)$proj4string,st_crs(streams_sf)$proj4string)
#verify CRS of transform


# +

plot(sf::st_geometry(streams_sf), col= "dodgerblue", main = "Species Presence locations")
plot(sf::st_geometry(utah) , border = "black", lwd =2, add = TRUE)
plot(sf::st_geometry(dat), pch=17,col= "green", add = TRUE)


# +
### 2.1.5.a Point correction
# -



# +

###2.1.6 Point Correction
##########################################################
##########################################################
#  Point Correction
#   moves points to respective stream polyline
#   ensure all layers are in common crs
##########################################################
library(maptools)
library(rgeos)

# +
# streams_sp <- as_Spatial(streams_sf) #convert to SpatialPoint
# pts.outSP <- as_Spatial(dat) 

#################start testing
streams_sf <- streams_sf[ ! st_is_empty( streams_sf ) , ]
streams_sp <- as_Spatial(streams_sf) #convert to SpatialPoint
dat <- dat[ ! st_is_empty( dat ) , ]
pts.outSP <- as_Spatial(dat) 
#####################end testing
# -

#crs(streams_sp)
#crs(pts.outSP)
pts_snap<-snapPointsToLines(pts.outSP,    #The object(points) to be 
                            #corrected
                            streams_sp,   #The nhd streams sp object
                            maxDist=2000) #The maxmimum distance a point 
#can be moved.  
pts_snap_sf <-sf::st_as_sf(pts_snap)              
#The snapped point object being converted back to an sf object.
#                 coords = station@coords, 
#                 crs = prj.utmN83z12, 
#                remove = F)
pts_snap_sf<- sf::st_transform(pts_snap_sf, 
                               crs = prj.aeaN83)
##2.3 Species Data Organization
##2.3.1 Establish SDHM output resolution 
streams_pts_sgbp <-st_intersects(streams_sf,            
                                 #The polyline(stream layer in this case)
                                 st_buffer(pts_snap_sf, 
                                           #applying a buffer to the snapped points.
                                           dist = 200), 
                                 #distance of buffer around point
                                 sparse = TRUE)
streams_pts_sgbp_1 <- lengths(streams_pts_sgbp)> 0   
#Takes the sgbp list and removes segments that do not intersect
stream_presence <-streams_sf[streams_pts_sgbp_1,]    
#Uses the previous sgbp to select those segments that do intersect.  

########################################
stream_presence <- stream_presence %>% 
  distinct(comid, .keep_all = TRUE)
########################################

head(stream_presence)

# +

plot(stream_presence$geometry,col = "dodgerblue", main = "Corrected/Transformed Species Presence Segments")
plot(utah$geometry, add = TRUE)

# +

###2.3.3
##########MODELING DOMAIN AND BOUNDING BOX DEVELOPEMENT
### Create MCP
pts_snap.mcpsf <- sf::st_convex_hull(sf::st_union(pts_snap_sf))
# -

plot(utah$geometry)
plot(pts_snap.mcpsf, add = TRUE)
plot(pts_snap_sf$geometry, add = TRUE, pch = 18, col = "red")
##############################
#add drainage to polygon
#dissolve into huc8 sub drainages
huc_8_sample <-st_transformdrainages_sf %>%
  group_by(HUC_8) %>%
  summarize()
#giggle plot
plot(huc_8_sample$geometry)
plot(pts_snap.mcpsf, add = TRUE, col = "blue")
plot(pts_snap_sf$geometry, add = TRUE, pch = 18, col = "red")

sample_huc <- st_overlaps(huc_8_sample, pts_snap.mcpsf)
huc_frame <- lengths(sample_huc)> 0   
#Takes the sgbp list and removes segments that do not intersect
stream_sample_huc <-huc_8_sample[huc_frame,] 

stream_sample_huc <-stream_sample_huc %>%
  mutate(together = 123)%>%
  group_by(together) %>%
  summarize()
library(smoothr)
stream_sample_huc <- fill_holes(stream_sample_huc,area_thresh <- units::set_units(100000, km^2))
#giggle plot
plot(stream_sample_huc$geometry)


# +

#THIS IS THE MODELING DOMAIN WE CHOSE 
plot(utah$geometry, main = "MCP Vs HUC 10 Modeling Domain")
plot(streams_sf$geometry)
plot(stream_sample_huc$geometry, add = TRUE)
#plot(pts_snap.mcpsf, add = TRUE, col = "blue") #MCP plot if wanted
plot(pts_snap_sf, add = TRUE, pch = 18, col = "red")

# +
#######################################################
#End frame development


# +

###2.3.4
## Generate Presence/Pseud-Absences data frame
##################################################################
# -

#new code with drainage huc selection by mcp
sample <- st_covered_by( streams_sf, stream_sample_huc)
streams_frame <- lengths(sample)> 0   
#Takes the sgbp list and removes segments that do not intersect
stream_sample <-streams_sf[streams_frame,] 

plot(utah$geometry)
#plot(pts_snap.mcpsf, 
#     add = TRUE, 
#    col = "green")
plot(stream_sample$geometry, col = "dodgerblue", add = TRUE)
plot(pts_snap_sf$geometry, add = TRUE, pch = 18, col = "red")
plot(stream_sample_huc$geometry, add = TRUE)

#This selects the points that do not contain the species. denoted by the [-x]
sample2 <- st_covered_by( stream_sample, stream_presence)
streams_frame2 <- lengths(sample2)>0   
#Takes the sgbp list and removes segments that do not intersect
stream_sample2 <-stream_sample[!streams_frame2,]
#select the points that are absences.

#giggle plot
#plot(utah$geometry, main = " Presence/Absence Stream Segments")
plot(stream_sample_huc$geometry,main = " Presence/Absence Stream Segments", 
     col=rgb(red=.30, green=0.2, blue=1.00, alpha=0.4))
#plot(pts_snap.mcpsf, col=rgb(red=1.0, green=0.2, blue=.50, alpha=0.2), add = TRUE)
plot(stream_sample2$geometry, col = "dodgerblue",add= TRUE, lwd= .5 )
plot(stream_presence$geometry, col = "yellow", add = TRUE, lwd = 1)
#MAY ADD A INDIVIDUAL DRAINAGE TRANSPARENT LAYER
#rm(stream_sample)
#gc()
#data frame merge
#assigns a 1 to segments that are present
tryme <-stream_presence %>%
  mutate(presabs = 1)
#assigns points of absence in frame to 0 or not known presence
tryme2 <- stream_sample2 %>%
  mutate(presabs = 0)
complete_frame_streams <- rbind(tryme,tryme2)
#balance checks
nrow(complete_frame_streams)#should be same
nrow(stream_sample)

#giggle me plot plot ugly though
plot(utah$geometry)
plot(tryme2, add = TRUE, col = "dodgerblue") #absence points
plot(tryme, add = TRUE, col = "red")  #presence points

###Selecting Pseudo-absences from mcp non present confirmed points
#The set seed function must be set here to repeat results for selection
#The number of pseudo absence points that will be drawn from the background 
#will be equal to the number of grid cells with occurrence points for RF and 
#BRT.  For GLM, GAM, and MAXENT, at least 10x the number of grid cells with 
#occurrence points will be drawn.  This approach has been shown to maximize 
#sensitivity or to minimize the possibility of an omission error in the model 
#output (Barbet-Massin et al. 2012). 
set.seed(123)
sample.pseudo <- tryme2[sample(1:nrow(tryme2), 4 * nrow(tryme), replace = F), ] #the 2 denotes the multiple of samples compared to presence points

plot(sample.pseudo$geometry)
plot(tryme$geometry, col = "dodgerblue", add = TRUE)
nrow(sample.pseudo)
nrow(tryme)


# +

###2.3.5 Construct the SDHM training data frame
#build database of presence/absences
tryme3 <- rbind(tryme, sample.pseudo)
nrow(tryme3) #ensure number of rows equals that from above

# +
######################################################################
#end database construction
##########

# +

######################################################################
######################################################################
#   Data Exploration and covariate developement
######################################################################
# -


sample_subset <-tryme3[sapply(tryme3, function(x) is.numeric(x))] #pull only numeric columns of data for correlat
nrow(sample_subset)


# +

presabs <- sample_subset[,"presabs"] #reorder columns so presabs is first item
tester <- dplyr::select(sample_subset, -presabs)
# -



sample_subset_2 <-cbind(st_drop_geometry(presabs),tester)
sample_subset_22<-cbind(st_drop_geometry(presabs),st_drop_geometry(tester))



options (warn = - 1)
###############################################################################
######## START CORRELATIONS AMONG TOPO PREDICTORS
# numeric correlations
cut.point <- 0.7 # set cutpoint for correlation
c1 <- cor(sample_subset_2[, c(1:ncol(sample_subset_2)-1)], use = "pairwise.complete.obs",  #note the geometry must be dropped to run correlation
          method = "spearman") # est. correlation
#c1 # examine
c2 <- subset(c1 > cut.point | c1 < -cut.point) # matrix of cor>cutpoint
#c2 # examine; FALSE indicates cor<cutpoint 

# START MODIFIED panel.cor CORRELATION FUNCTION		   
#   determine correlations among predictor variables: modified from 
#   http://addictedtor.free.fr/graphiques/graphcode.php?graph=137
panel.cor <- function(x, y, digits=2, prefix="", cex.cor) 
{ usr <- par("usr"); on.exit(par(usr)) 
par(usr = c(0, 1, 0, 1)) 
r <- abs(cor(x, y, use = "pairwise.complete.obs")) 
txt <- format(c(r, 0.123456789), digits=digits)[1] 
txt <- paste(prefix, txt, sep="") 
if(missing(cex.cor)) cex <- 0.8/strwidth(txt) 
test <- cor.test(x,y) 
# borrowed from printCoefmat
Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                 cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                 symbols = c("***", "**", "*", ".", " ")) 
text(0.5, 0.5, txt, cex = cex * r)
text(.8, .8, Signif, cex=cex, col=2) 
}
# END MODIFIED panel.cor CORRELATION FUNCTION	

# plot correlations using modified panel.cor function
#png(filename = "Variable Assmt", width = 1000, height = 1000) 
pairs(st_drop_geometry(tryme3[, c(20:30)]), lower.panel = panel.smooth, 
      upper.panel = panel.cor, main = "Variable Correlation Assessment") 
#dev.off()
#######################################################################
#exploration
#######################################################################


#### START function variable importance
varimp.glm <- function(tr.spp, tr.var, pres, pf, pl) {
  tmp.mat <- matrix(ncol = 2, nrow = (pl - pf + 1))
  for (i in pf:pl) {
    # option: linear+quadratic; linear only
    tmp <- glm(tr.spp[, 1] ~ tr.var[, i] + I((tr.var[, i])*(tr.var[, i])), na.action = na.omit, 
               family = binomial)
    # linear only glm
    #tmp <- glm(tr.spp[, pres] ~ tr.var[, i], na.action = na.omit, family = binomial)
    tmp.mat[(i - pf + 1), 1] <- tmp$aic
    tmp.mat[(i - pf + 1), 2] <- (1 - (tmp$deviance/tmp$null.deviance))
  }
  return(tmp.mat)
} 

options(scipen = 999)

tr.vip <- as.data.frame(sample_subset_22[,1:ncol(sample_subset_22)]) # keep only P/A & predictors
#tr.vip<- tr.vip[which(rowMeans(!is.na(tr.vip[,1:ncol(tr.vip)])) > 0), which(colMeans(!is.na(tr.vip)[,1:ncol(tr.vip)]) > 0)]
tr.vip <- mutate_all(tr.vip, function(x) as.numeric(as.character(x)))
head(tr.vip)

tr.spp <- sample_subset_2[,1:ncol(sample_subset_2)]
pres <- 1 # column for presence:absence
v.start <- 2 # column start predictor variables
v.stop <- ncol(tr.vip) # last column predictor variables
v.num <- v.stop - 8 # number predictor variables


# +
#######################################################################
#######################################################################
########################################################################
#this is my problem line! help me
#############################################
#############################################
##############################################
###############################################

dev.fit <- varimp.glm(tr.vip, tr.vip, pres, v.start, v.stop)# call VIP function

# +
#dev.fit <- dev.fit[-146,]#this line fixes my plot but why
###############################################################
################################################################
################################################################
# -

#names wont work because 146 was removed
# built basic barplot if desired
d.max <- ceiling(signif(max(dev.fit[, 2]), 2) * 10)/10 # max of y-axis
ylim.r <- range(0, d.max) # range y-axis

x.labs <- names(tr.vip[,2:ncol(tr.vip)]) # x-axis labels THe minus 1 removes the second geometry!



# +

#png(filename = "allvip.png", width = 550, height = 400)
barplot(dev.fit[, 2], col = "grey", ylim = ylim.r, main = "VIPs", 
        ylab = "adj.D2", names = x.labs) # barplot
abline(h = 0) # add horizontal line
abline(mean(dev.fit[, 2]), 0, lt = 3) # ref lines; dash=mean adj.dev 
#dev.off()



# +

options(scipen = 999)
#filter and build a barplot with only top 26
library(dplyr)
predictors_vimp <-cbind(names(tr.vip[2:v.stop]),dev.fit)
predictors_vimp <-as.data.frame(predictors_vimp)
colnames(predictors_vimp) = c("preds", "col2", "col3")
predictors_vimp <- predictors_vimp %>%
  arrange(desc(col3))
predictors_vimp$col3 <- as.numeric(as.character(predictors_vimp$col3))
head(predictors_vimp, 40)
# -

selected_preds <-predictors_vimp[,1]#adds names to draw from later for directed correlation

predictors_vimp <- predictors_vimp[1:35,]
predictors_vimp <- predictors_vimp[2:35,]
d.max <- ceiling(signif(max(predictors_vimp$col3), 2) * 10)/10 # max of y-axis
ylim.r <- range(0, d.max) # range y-axis
x.labs <- predictors_vimp$preds # x-axis labels


barplot(predictors_vimp$col3, col = "dodgerblue", ylim = ylim.r, main = "VIPs", 
        ylab = "adj.D2", 
        names = x.labs, # barplot
        cex.names=.7, las = 2.2)
#abline(h = 0) # add horizontal line
abline(mean(predictors_vimp[, 3]), 0, lt = 3) # ref lines; dash=mean adj.dev



# +

#png(filename = "vip.png", width = 550, height = 400)
barplot(predictors_vimp$col3, col = "dodgerblue", ylim = ylim.r, main = "VIPs", 
        ylab = "adj.D2", 
        names = x.labs, # barplot
        cex.names=.70, las = 2.2)
abline(h = 0) # add horizontal line
abline(mean(predictors_vimp[, 3]), 0, lt = 3) # ref lines; dash=mean adj.dev
#dev.off()  


# +

cov_names <- unique(predictors_vimp$preds)# use highest importance variables
subset_covariates_1 <-sample_subset_2[,c("presabs",cov_names)]
# -

st_write(subset_covariates_1,  "covariates1.shp", 
         delete_dsn=FALSE, 
         update = TRUE, 
         driver = "ESRI Shapefile")

# +
#writeOGR(as_Spatial(subset_covariates_1), ".", "data/gis_layers/subsetcovarR", driver="ESRI Shapefile")
# -

#database without presence absence
stream_sample3 <-complete_frame_streams[,cov_names]
geom_stream_sample <- complete_frame_streams$geometry

stream_sample3 <-cbind(stream_sample3, geom_stream_sample)

stream_sample3.ng <- stream_sample3

# load libraries now if desired; loaded below when needed
library(PresenceAbsence)  # fxns: optimal.thresholds, presence.absence.accuracy, 
#       auc.roc.plot
library(DAAG) 


library(pander)
Type <- c("Parametric", "Semi-parametric", "Machine Learning", "Machine Learning", "Machine Learning")
Model <- c("Logistic Regression", "Generalizard Additive", "Maximum Entropy", "Boosted Regression Trees", "Random Forest")
pander::pander(cbind.data.frame(Model, Type))


mod.form <- function(dat, r.col, p.col) {
  # generic glm formula construction function; inputs as: 
  #  resp =>col 1 in dataframe; preds=>col 2 thru ncol such that p.col=2 
  #  NOTE: vars as factors must be coerced PRIOR to formula construction 
  # example call: mod.form(dat1, 1, 2)
  n.col <- ncol(dat) # No. columns in dataframe
  resp <- colnames(dat[r.col]) # assign resp column name
  resp <- paste("as.factor(", colnames(dat[r.col]), ")", sep = "") # assign resp column name
  pred <- colnames(dat[c(p.col:n.col)]) # assign preds column names
  mod.formula <- as.formula(paste(resp, "~", paste(pred, collapse = "+"))) # build formula
} 

easy_one <-mod.form(subset_covariates_1,1,2)

# +

################################################################################
######## START MODEL #1 => FULL MODEL, ALL VARS INCLUDED
# -

dat1 <-sample_subset_2

mod1.LR <- glm(easy_one,family = binomial, data = sample_subset_2)

#### build initial model all variables: mod1.LR
# hard code version
#cov_list <- colnames(stream_sample3.ng)
cov_list <- colnames(stream_sample2)

# +

#mod1.LR <- glm(mod.form(dat1, 1, 2), family = binomial, data = dat1)
# -

# model 1 summary
summary(mod1.LR) # full model summary stats 

# model 1 fit
mod1.fit <- 100 * (1 - mod1.LR$deviance/mod1.LR$null.deviance) # model fit
mod1.fit  # examine fit
mod1.pred <- predict(mod1.LR,stream_sample2, type = "response") # model prediction

# +
#Giggle Plot
# -

mod1.LR.pred <- cbind(stream_sample2,mod1.pred)
#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('grey','green'))  #This adds a column of color values
mod1.LR.pred$Col <- rbPal(10)[as.numeric(cut(mod1.LR.pred$mod1.pred,breaks = 10))]


plot(mod1.LR.pred$geometry,col = mod1.LR.pred$Col) #predicted presence
plot(tryme$geometry, add = TRUE, col = "yellow")  #verified presense

# +
####### END MODEL #1 => FULL MODEL, ALL VARS INCLUDED
###############################################################################
# -


################################################################################
######## START MODEL #2 => REDUCED VARIBALE MODEL
# build parsimonious model w/var reduction techniques: mod2.LR
#   variable reduction: backwards
mod2.LR <- step(mod1.LR, trace = F) # backwards stepwise variable reduction
mod2.fit <- 100 * (1 - mod2.LR$deviance/mod2.LR$null.deviance)
mod2.fit # model fit

# model 1 v. model 2 fit
100 * (1 - mod1.LR$deviance/mod1.LR$null.deviance)  # fit model 1
100 * (1 - mod2.LR$deviance/mod2.LR$null.deviance)  # fit model 2

# model 2 prediction
mod2.pred <- predict(mod2.LR,type = "response") # model prediction
head(mod2.pred)

# model 1 v. model 2 prediction
head(mod1.pred) # mod 1 prediction
head(mod2.pred) # mod 2 prediction

# model 2 summary
summary(mod2.LR) # reduced model summary 



######## END MODEL #2 
################################################################################
################################################################################
######## START RESUBSTITUTION ACCURACY CALCULATIONS, MODEL=LOGISTIC GLM
# requires pkg PresenceAbsence
#   build testing dataframe using mod2 predictions
modl <- "mod2.LR" # add var to keep track of model
dat2 <- cbind(modl, dat1[1], mod2.pred) # build dataframe w/mod2 predictions
head(dat2, 20) # examine prediction dataframe

# determine best threshold using PresenceAbsence package Sec7.1
#library(PresenceAbsence) # PresenceAbsence for accuracy metrics
# help(optimal.thresholds) # options for optimizing threshold
mod.cut <- optimal.thresholds(dat2, opt.methods = c("ObsPrev")) # default threshold=0.5
mod.cut # examine threshold=DEFAULT of 0.5
modF.cut <- numeric(0)
modF.cut$LR.cut <- mod.cut


# generate confusion matrix
mod2.cfmat <- table(dat2[[2]], factor(as.numeric(dat2$mod2.pred >= mod.cut$mod2.pred)))
mod2.cfmat # examine

# calculate model accuracies with standard deviation=F
mod2.acc <- presence.absence.accuracy(dat2, threshold = mod.cut$mod2.pred, st.dev = F)
tss <- mod2.acc$sensitivity + mod2.acc$specificity - 1 # code TSS metric
mod2.acc <- cbind(mod2.acc[1:7], tss) # bind all metrics
mod2.acc[c(1, 4:5, 7:8)] # examine accuracies

# plotting AUC
auc.roc.plot(dat2, color = T) # basic AUC plot; pkg PresenceAbsence 
gc()


mod2.pred <- predict(mod2.LR,stream_sample2, type = "response")

#Bind the predictions to the dataframe then produce map product
mod2.LR.pred <- cbind(stream_sample2,mod2.pred)

#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('blue','green'))  #This adds a column of color values
mod2.LR.pred$Col <- rbPal(10)[as.numeric(cut(mod2.LR.pred$mod2.pred,breaks = 10))]

plot(utah$geometry, lwd = 2, col=rgb(red=0, green=0.5, blue=1, alpha=0.2 ), main = "Logistic Regression Final Prediction") #plot(huc_8_sample, add = TRUE, col=rgb(red=.5, green=0.5, blue=.50, alpha=0.1 ))
plot(mod2.LR.pred$geometry,col = mod2.LR.pred$Col, add = TRUE)
plot(tryme$geometry, add = TRUE, col = "orange")

save(mod2.LR, file="lr.Rdata")
# save plot if desired 
#setwd(path.figs)
#savePlot(filename = "LR_pred.pdf", type = "pdf")
######## END RESUBSTITUTION ACCURACY CALCULATIONS, MODEL=LOGISTIC GLM
################################################################################

mod.cut #view probability cut for model

slsc_layer <- stream_sample2 #create feauture layer format
slsc_layer$lr.prob <- mod2.pred #add probability values to table
#head(slsc_layer,25) #view layer 

# +

##Generalized Linear Model
# -

library(gam)

################################################################################
######## START SOME PRELIMINARY FIGURES
# some scatterplot of remote sensing-based variables
plot(subset_covariates_1[2:20]) # simple data scatterplots

# +
# NOTE: car package has great graphs BUT conflicts w/plots of gam models
#  if load car MUST start new R session to use plots of gam models
#   conversely, run library car in new R window and build plots there
#  even detach(package:car) does not solve this; a true R bug ...
# library(car)  # NOT RUN; fancy scatterplot w/lines if desired
# scatterplotMatrix(dat1[2:6], main = "Scatterplots of RS data")  # scatterplot

# +

####### END SOME PRELIMINARY FIGURES
###############################################################################0
# -

mod.form(subset_covariates_1, 1, 2)


# GAM model 0: all linear
#   can use glm for linear model 
mod0.GAM <- glm(easy_one, family = binomial, data = sample_subset_2)


summary(mod0.GAM) # model #0 summary



# +

mod.form2 <- function(dat, r.col, p.col) {
  # generic glm formula construction function; inputs as: 
  #  resp =>col 1 in dataframe; preds=>col 2 thru ncol such that p.col=2 
  #  NOTE: vars as factors must be coerced PRIOR to formula construction 
  # example call: mod.form(dat1, 1, 2)
  n.col <- ncol(dat) # No. columns in dataframe
  resp <- colnames(dat[r.col]) # assign resp column name
  resp <- paste("as.factor(", colnames(dat[r.col]), ")", sep = "") # assign resp column name
  pred <- colnames(dat[c(p.col:n.col)]) # assign preds column names
  mod.formula <- as.formula(paste(resp, "~",  paste("s(",pred,")", collapse = "+"))) # build formula
} 
# -

mod.form3 <- function(dat, r.col, p.col) {
  # generic glm formula construction function; inputs as: 
  #  resp =>col 1 in dataframe; preds=>col 2 thru ncol such that p.col=2 
  #  NOTE: vars as factors must be coerced PRIOR to formula construction 
  # example call: mod.form(dat1, 1, 2)
  n.col <- ncol(dat) # No. columns in dataframe
  resp <- colnames(dat[r.col]) # assign resp column name
  resp <- paste("as.factor(", colnames(dat[r.col]), ")", sep = "") # assign resp column name
  pred <- colnames(dat[c(p.col:n.col)]) # assign preds column names
  mod.formula <- as.formula(paste(resp, "~",  paste("s(",pred,",5)", collapse = "+"))) # build formula
} 

easy_one2 <-mod.form2(subset_covariates_1,1,2)
easy_one3 <-mod.form3(subset_covariates_1,1,2)

# +

easy_one2 <-(as.factor(presabs) ~ s(OBJECTID) + s(LevelPathI) + s(FromNode) + 
               s(Pathlength) + s(TerminalPa) + s(UpLevelPat) + s(UpHydroseq) + 
               s(ToNode) + s(DnLevel) + s(StreamLeve) + s(tmin8110ws) + 
               s(qe_ma) + s(StreamOrde) + s(StreamCalc) + s(tmean8110w) + 
               s(qa_ma) +  s(DnLevelPat) + s(DnHydroseq) + s(tmin8110ca) + 
               s(qc_ma) + s(pcturbaniz) + s(va_ma) + s(SLOPE_1) + s(vc_ma) + 
               s(ve_ma) + s(mgows) + s(DivDASqKM) + s(ArbolateSu) + s(TotDASqKM) + 
               s(tmean8110c) + s(tmax8110ws) + s(rddenscat) + s(rddenscatr))
# -


easy_one3 <-(as.factor(presabs) ~ s(OBJECTID, 5) + s(LevelPathI, 5) + s(FromNode, 
                                                                        5) + s(Pathlength, 5) + s(TerminalPa, 5) + s(UpLevelPat, 
                                                                                                                     5) + s(UpHydroseq, 5) + s(ToNode, 5) + s(DnLevel, 5) + s(StreamLeve, 
                                                                                                                                                                              5) + s(tmin8110ws, 5) + s(qe_ma, 5) + s(StreamOrde, 5) + 
               s(StreamCalc, 5) + s(tmean8110w, 5) + s(qa_ma, 5) 
             + s(DnLevelPat, 5) + s(DnHydroseq, 5) + s(tmin8110ca, 
                                                       5) + s(qc_ma, 5) + s(pcturbaniz, 5) + s(va_ma, 5) + s(SLOPE_1, 
                                                                                                             5) + s(vc_ma, 5) + s(ve_ma, 5) + s(mgows, 5) + s(DivDASqKM, 
                                                                                                                                                              5) + s(ArbolateSu, 5) + s(TotDASqKM, 5) + s(tmean8110c, 5) + 
               s(tmax8110ws, 5) + s(rddenscat, 5) + s(rddenscatr, 5))

formatGam <- function(dat) {
  resp <- colnames(dat[1]) # assign resp column name
  pred <- colnames(dat[c(4:ncol(dat)-1)]) # assign preds column names # -1 to remove geom
  mod.form <- as.formula(paste(as.factor(resp), "~", paste(paste("lo(",pred,",5)"), collapse = "+"))) # formula
}
head(sample_subset_2)

# +

# GAM model 1: all smoothers using defaults df (smoothers=4)
mod1.GAM <- gam(formatGam(sample_subset_2), family = binomial, data = sample_subset_2)
#summary(mod1.GAM) # model #1 summary
summary(mod2.LR)



# +
# # GAM model 2: all smoothers using specified df (smoothers=5)
# date()
# mod2.GAM <- gam(easy_one3, family = binomial, data = sample_subset_2)
# date()
# summary(mod2.GAM) # model #2 summary 

# +
# ###############################################################################
# ####### START FINAL GAM MODEL OUTPUT
# #build GAM parameter table output; mod2.GAM chosen based on fit / accuracies
# names(mod2.GAM$coefficients) # find the model terms
# -

# build model
modfin1.gam <- gam(formatGam(sample_subset_2), family = binomial, data = sample_subset_2)

# +
# final model summary
table <- summary(modfin1.gam)
names(modfin1.gam$coefficients) # determine which vars were smoothed vs. linear

names(table)
table
summary.glm(modfin1.gam)$coefficients[c(4:5), ] # access linear terms

# -

test<- na.omit(table$anova)
row.names(test)
new <- str_sub(row.names(test),1,nchar(row.names(test))-4)
substring(new, 4)

# final model fit and accuracy
modfinl.fit <- 100 * (1 - modfin1.gam$deviance/modfin1.gam$null.deviance) # model fit
modfinl.fit

# final model resubstituion accuracy
modl <- "modfinl.GAM" # add var to keep track of model
modfinl.pred <- predict(modfin1.gam, type = "response") # predict by model
dat2 <- cbind(modl,
              #sample_subset_2[1], 
              modfinl.pred)  # build dataframe w/modfinl predictions
dat2 <- cbind(dat2, subset_covariates_1[1]) #requires the cbind to be done in two parts
dat2 <- dat2[,c("modl","presabs","modfinl.pred" )] #reorder the features


mod.cut <- optimal.thresholds(dat2, opt.methods = c("MaxSens+Spec")) # threshold
mod.cut # examine threshold
modF.cut$gam.cut <- mod.cut
# generate resubstitution confusion matrix
modfinl.cfmat <- table(dat2[[2]], 
                       factor(as.numeric(dat2$modfinl.pred >= mod.cut$modfinl.pred)))
modfinl.cfmat # examine

# calculate model resubstitution accuracies with standard deviation=F
modfinl.acc <- presence.absence.accuracy(dat2, threshold = mod.cut$modfinl.pred, st.dev = F)
tss <- modfinl.acc$sensitivity + modfinl.acc$specificity - 1 # code TSS metric
modfinl.acc <- cbind(modfinl.acc[1:7], tss) # bind all metrics
modfinl.acc[c(1, 4:5, 7:8)] # examine resubstituion accuracies

# final model 5-fold cross-validation accuracy
cv.modfinl <- CVbinary(modfin1.gam, nfolds = 5, print.details = F) # crossval predict
cv.predfinl <- cv.modfinl$cvhat # assign new name to jackknife estimates
cv.datfinl <- cbind(modl, dat1[1], cv.predfinl) # build obs and prediction dataframe
cv.cutfinl <- optimal.thresholds(cv.datfinl, opt.methods = c("MaxSens+Spec")) # threshold

# generate 5-fold cross-validation  confusion matrix
cv.cfmatfinl <- table(cv.datfinl[[2]], 
                      factor(as.numeric(cv.datfinl$cv.predfinl > cv.cutfinl$cv.predfinl))) # confusion 

# calculate model 5-fold cross-validation accuracies with standard deviation=F
cv.accfinl <- presence.absence.accuracy(cv.datfinl, 
                                        threshold = cv.cutfinl$cv.predfinl, st.dev = F) # calculate accuracies
tss <- cv.accfinl$sensitivity + cv.accfinl$specificity - 1 # code TSS metric
cv.accfinl <- cbind(cv.accfinl[1:7], tss) # bind all metrics
cv.accfinl$model <- modl # variable substitution
cv.accfinl[c(1, 4:5, 7:8)] # examine accuracies

# GAM plots using interactive mode
#  vars to plot here are: [2:3,5:8]
par(mfrow = c(2, 3)) # set par

# +
#plot(mod.finl, ask = T) # interactive graphs

# +
# plot categorical var geonutr ; plot number=3
# par(mfrow = c(1, 1)) # reset par


# +
# plot(modfin1.gam, ask = T) # interactive graphs


# +

####### START FINAL GAM MODEL OUTPUT
###############################################################################

# +
################################################################################
######## START PREDICTION & RECLASSIFICATION OF FINAL GAM MODEL
# NOT RUN: basic prediction code to build probability map
# predict for final GAM model
#   ASSUMEs stream_sample3.ng is a table 

#remove null columns from df
stream_sample3.ng <- stream_sample3
stream_sample3.ng <- stream_sample3.ng[, (colnames(stream_sample3.ng) %in% substring(new, 4))]
##
stream_sample3.ng <- stream_sample2 
stream_sample3.ng 
modfin1.gam
modfinl.GAM <- predict( modfin1.gam,stream_sample3.ng,
                        type = "response")
# -

#Bind the predictions to the dataframe then produce map product
modfinl.GAM.pred <- cbind(stream_sample2,modfinl.GAM)
#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('blue','green'))  #This adds a column of color values
modfinl.GAM.pred$Col <- rbPal(10)[as.numeric(cut(modfinl.GAM.pred$modfinl.GAM,breaks = 10))]

plot(utah$geometry, lwd = 2, col=rgb(red=0, green=0.5, blue=1, alpha=0.2 ), main = "GAM Final Prediction") 
#plot(huc_8_sample, add = TRUE, col=rgb(red=.5, green=0.5, blue=.50, alpha=0.1 ))
plot(modfinl.GAM.pred$geometry,col = modfinl.GAM.pred$Col, add = TRUE)
plot(tryme, add = TRUE, col = "yellow")

save(modfin1.gam, file="gam.Rdata")
######## END PREDICTION & RECLASSIFICATION OF FINAL GAM MODEL
################################################################################


# +
##Construct Feature Table
# -


mod.cut #view probability cut for model

slsc_layer <- stream_sample2 #create feauture layer format
slsc_layer$gam.prob <- modfinl.GAM.pred$modfinl.GAM #add probability values to table
#head(slsc_layer,25) #view layer 

# +
### 3.2.2 Entropy models (MAXENT)

# +

mod.form4 <- function(dat, s.col) {
  # generic glm formula construction function; inputs as: 
  #  resp =>col 1 in dataframe; preds=>col 2 thru ncol such that p.col=2 
  #  NOTE: vars as factors must be coerced PRIOR to formula construction 
  # example call: mod.form(dat1, 1, 2)
  n.col <- ncol(dat) # No. columns in dataframe
  pred <- colnames(dat[c(s.col:n.col)]) # assign preds column names
  pred <- colnames(dat[c(s.col:n.col)]) # assign preds column names
}
# -

easy_one4<-mod.form4(subset_covariates_1,1)

library(dismo)
library(raster)
library(rJava)
#the presence only dataframe tryme


max.preds.1 <- st_drop_geometry(tryme3[,c(easy_one4)])

# save 20% presence for testing, 80% training; 5-fold x-val
set.seed(1234) # set.seed if desire repeatability
#library(dismo) # dismo needed for kfold and maxent
fold <- kfold(max.preds.1, k = 5) # data k-folds == 5
pres.tst <- max.preds.1[fold == 1, ] # build test data; the 20%
pres.tr <- max.preds.1[fold != 1, ] # build training data; the 80%


pa <-pres.tr[,"presabs"]

easy_one5<-mod.form4(subset_covariates_1,2) #grab preds only(start on column 2)

max.preds <- pres.tr[,c(easy_one5)]


mod1.MAX <- maxent( max.preds, #training data (pred.dom)
                    pa)        #presabs vector (pres.tr)

# DEPENDING ON SYSTEM:  
#  typing "mod1.MAX" in R redirects output to browser
# OR
#  browser automatically opens
mod1.MAX

plot(mod1.MAX) # var importance plot
# save plot if desired 
#setwd(path.figs)
#savePlot(filename = "mod7fig01.pdf", type = "pdf")

response(mod1.MAX) # prediction vs. var plot
# save plot if desired 
# setwd(path.figs)
# savePlot(filename = "mod7fig02.pdf", type = "pdf")

#### NOT RUN; assessment if desire maxent assessment only
mod2.MAX <- mod1.MAX # copy for maxent evaluation
# assessment #1
#mod2.bak <- randomPoints(pred.dom, 1000)
#mod2.val <- evaluate(mod2.MAX, p = pres.tst, a = mod2.bak, x = pred.dom)
#mod2.val # examine
# assessment #2
#pts.tst <- data.frame(extract(pred.dom, pres.tst))
#pts.bak <- data.frame(extract(pred.dom, mod2.bak))
#mod2.val <- evaluate(mod2.MAX, p = pts.tst, a = pts.bak)
####

# +

library(PresenceAbsence)
# -


#filter presence and absence for model evaluation  
pres.tst <- subset_covariates_1 %>% 
  filter(presabs == 1)

bak.xy <- subset_covariates_1 %>% 
  filter(presabs == 0) 


# evaluate model (an x-fold process: see help(evaluate))
mod1.val <- evaluate(mod1.MAX, p = pres.tst, a = bak.xy)
#x = pred.dom) # x-fold cross-val
mod1.val # examine
threshold(mod1.val) # some validation information

# +

# modl <- "mod1.MAX"
# mod1.pred <- pres.pred # tmp assignment
# tmp.p <- cbind(modl, subset_covariates_1[,1], mod1.pred) # pres dataframe !! spp in col=3
# mod1.pred <- bak.pred # tmp assignment
# tmp.b <- cbind(modl, bak.xy[3], mod1.pred) # bak dataframe
# dat2 <- data.frame(rbind(tmp.p, tmp.b)) # final dataframe
# head(dat2) # examine 
# -

# extract thresholds; see help(threshold) for options
#   eg, kappa is "max kappa"
threshold(mod1.val)[[1]] # extract max kappa
threshold(mod1.val)[[2]] # returns threshold spec_sens
mod.cut <- threshold(mod1.val) # extract max kappa
mod.cut # view maxent thresholds
#modF.cut$max.cut <- mod.cut$spec_sens
modF.cut$max.cut <- threshold(mod1.val)[[2]]
#optimal.thresholds(subset_covariates_1, opt.methods = 1:6) # PresenceAbsence thresholds 

mod.cut

# closest comparion thresholds
mod.cut[c(1:2)] # maxent thresholds via dismo
#optimal.thresholds(subset_covariates_1, opt.method = 3) # PresenceAbsence thresholds   

# +
#################
#################this is the problem code section
# -

mod1.cfmat <- table(subset_covariates_1[[1]], factor(as.numeric(dat2$modfinl.pred >= mod.cut$spec_sens)))
mod1.cfmat # examine



# calculate model accuracies with standard deviation=F; assume PresenceAbsence
mod1.acc <- presence.absence.accuracy(dat2, threshold = mod.cut$spec_sens, st.dev = F)
tss <- mod1.acc$sensitivity + mod1.acc$specificity - 1 # code TSS metric
mod1.acc <- cbind(mod1.acc[1:7], tss) # bind all metrics
mod1.acc[c(1, 4:5, 7:8)] # examine accuracies 

# plotting AUC
auc.roc.plot(dat2, color = T) # basic AUC plot; pkg PresenceAbsence
# save plot if desired 
#setwd(path.figs)
#savePlot(filename = "mod7fig03.pdf", type = "pdf")
######## END  MAXENT MODEL #1: TRUE PRESENCE ONLY DATA
################################################################################


################################################################################
######## START SPATIAL PREDICTION & CLASSIFICATION
# NOT RUN:  mod1 spatial prediction
mod1.MAXprob = predict(mod1.MAX, st_drop_geometry(stream_sample3.ng)) # predict entire model domain
#writeRaster(mod1.MAXprob, filename = "mod1.MAXprob.img", format = "HFA")

#Bind the predictions to the dataframe then produce map product
mod1.MAXprob <- cbind(stream_sample2,mod1.MAXprob)
#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('blue','red'))  #This adds a column of color values
mod1.MAXprob$Col <- rbPal(10)[as.numeric(cut(mod1.MAXprob$mod1.MAX,breaks = 5))]

plot(utah$geometry, lwd = 2, col=rgb(red=0, green=0.5, blue=1, alpha=0.2 ), main = "Maxent Final Prediction") 
#plot(huc_8_sample, add = TRUE, col=rgb(red=.5, green=0.5, blue=.50, alpha=0.1 ))
plot(mod1.MAXprob$geometry,col = mod1.MAXprob$Col, add = TRUE)
plot(tryme, add = TRUE, col = "yellow")
plot(mod1.MAXprob$geometry,col = mod1.MAXprob$Col)
plot(tryme, add = TRUE, col = "yellow")


save(mod1.MAX, file="maxent.Rdata")

###Construct Feature Table
mod.cut #view probability cut for model

slsc_layer$max.prob <- mod1.MAXprob$mod1.MAXprob #add probability values to table
#####################end maxent


# +
### 3.2.3 Classification-like models (CT, MDA)
### 3.2.4 Machine learning models (NNET, GARP)

# +
### 3.2.5 Averaging models (Random Forest, BRT)


# +
##############################################################
# Random Forest Model
# -

library(randomForest)

# build model formula as alternative hard code
mod.form <- function(dat ,r.col, p.col) {
  # generic formula construction function; inputs as:
  #  resp =>col 1 in dataframe such that r.col=1, 
  #  preds=>col 2 thru ncol in dataframe such that p.col=2
  #  NOTE: predictor vars as factors; coerce PRIOR to formula construction
  # example call: mod.form(dat1,1,2)
  n.col <- ncol(dat) # No. columns in dataframe
  resp <- colnames(dat[r.col]) # assign resp column name
  resp <- paste("as.factor(", colnames(dat[r.col]), ")", sep = "") # assign resp column name
  pred <- colnames(dat[c(p.col:n.col)]) # assign preds column names
  mod.formula <- as.formula(paste(resp, "~", paste(pred, collapse = "+"))) # build formula 
}
######## END INITIALIZATION
################################################################################

# +
###############################################################################
####### START RF MODEL #1
# RF model 1
# -

mod1.RF <- randomForest(mod.form(subset_covariates_1, 1, 2), importance = T, 
                        keep.forest = T, data = subset_covariates_1) # RF model w/mod.form fxn

mod1.pred <- predict(mod1.RF, type = "prob")[, 2] # predict from model
#head(mod1.pred) # examine 
######## END RF MODEL #1
################################################################################

modl <- "mod1.RF" # add var to keep track of model
dat2 <- cbind(modl, dat1[1], mod1.pred) # build dataframe w/mod1 predictions
head(dat2, 2) # examine prediction dataframe 

# determine best threshold using PresenceAbsence package Sec7.1
#   see help(optimal.thresholds) for more info
#library(PresenceAbsence) # PresenceAbsence for accuracy metrics
#help(optimal.thresholds) # options for optimizing threshold
mod.cut <- optimal.thresholds(dat2, opt.methods = c("PredPrev=Obs"), req.sens = 0.95)
mod.cut # sensitivity set at 0.95
#mod.cutK=optimal.thresholds(dat2,opt.methods=c('MaxKappa')); mod.cutK # MaxKappa option

# +
modF.cut$rf.cut <- mod.cut


# -

################################################################################
######## START RESUBSTITUTION ACCURACY CALCULATIONS, MODEL=RF
# build testing dataframe using model #1 predictions
modl <- "mod1.RF" # add var to keep track of model
dat2 <- cbind(modl, dat1[1], mod1.pred) # build dataframe w/mod1 predictions
head(dat2, 2) # examine prediction dataframe 


# determine best threshold using PresenceAbsence package Sec7.1
#   see help(optimal.thresholds) for more info
#library(PresenceAbsence) # PresenceAbsence for accuracy metrics
#help(optimal.thresholds) # options for optimizing threshold
mod.cut <- optimal.thresholds(dat2, opt.methods = c("PredPrev=Obs"), req.sens = 0.95)
mod.cut$pred # sensitivity set at 0.95
#mod.cutK=optimal.thresholds(dat2,opt.methods=c('MaxKappa')); mod.cutK # MaxKappa option

mod1.cfmat <- table(dat2[[2]], factor(dat2$mod1.pred >= mod.cut$mod1.pred))
mod1.cfmat # examine

# calculate model accuracies with standard deviation=F
mod1.acc <- presence.absence.accuracy(dat2, threshold = mod.cut$mod1.pred, st.dev = F)
tss <- mod1.acc$sensitivity + mod1.acc$specificity - 1 # code TSS metric
mod1.acc <- cbind(mod1.acc[1:7], tss) # bind all metrics
mod1.acc[c(1, 4:5, 7:8)] # examine accuracies

# plotting AUC
auc.roc.plot(dat2, color = T) # basic AUC plot; pkg PresenceAbsence 

# +
# save plot if desired 
# setwd(path.figs)
# savePlot(filename = "mod8fig04.pdf", type = "pdf")
####### END RESUBSTITUTION ACCURACY CALCULATIONS, MODEL=RF
###############################################################################


# +
##PLot RF model
# -


gc()

mod1.RF = predict(mod1.RF, stream_sample3.ng, type = "prob")[, 2]


mod1.rfprob <- cbind(stream_sample2,mod1.RF)
#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('blue','green'))  #This adds a column of color values
mod1.rfprob$Col <- rbPal(10)[as.numeric(cut(mod1.rfprob$mod1.RF,breaks = 10))]

plot(utah$geometry, lwd = 2, col=rgb(red=0, green=0.5, blue=1, alpha=0.2 ), main = "Random Forest Final Prediction") 
#plot(huc_8_sample, add = TRUE, col=rgb(red=.5, green=0.5, blue=.50, alpha=0.1 ))
plot(mod1.rfprob$geometry,col = mod1.rfprob$Col, add = TRUE)
plot(tryme, add = TRUE, col = "yellow")

save(mod1.RF, file="rf.Rdata")

# +
##Construct Feature Table
# -


mod.cut #view probability cut for model

#NLSC_feature_lyr <- stream_sample3 #create feauture layer format
slsc_layer$rf.prob <- mod1.RF #add probability values to table
head(slsc_layer,5) #view layer 
############end random Forest


# +
#####################################################################
# Boosted Regression Trees Model
# -


# load libraries now if desired; loaded below when needed
library(gbm)
#library(dismo)
#library(PresenceAbsence)

# +
###############################################################################
####### START INITIALIZATION
# -

dat3 <- as.data.frame(subset_covariates_1)


# BRT model formulation !! WARNING !! resp for datasets in col=4
dat3 <- na.omit(dat3) # remove NAs - BRT not picky but drop 'em anyway
resp <- paste("as.factor(", colnames(dat3[1]), ")", sep = "") # assign resp to col number
n.col <- ncol(dat3) # number of columns
pred <- (2:n.col) # assign predictors to column numbers 


# +
####### END INITIALIZATION
###############################################################################
# -

###############################################################################
####### START BRT MODELS
library(gbm) # load gbm package for BRT
library(dismo) # for BRT calls per Elith et al (2008) JAnimalEcol 77:802-813

# basic BRT model - example with warning
#   CPU time contingent on LR & TC;  ~30-35 min runtime w/LR=0.0001
date() # start time stamp
mod1.BRT <- gbm.step(data = dat3, gbm.x = pred, gbm.y = 1, family = "bernoulli",
                     tree.complexity = 15, learning.rate = 1e-04, bag.fraction = 0.75, n.folds = 10, 
                     n.trees = 5000, plot.main = TRUE, keep.fold.fit = TRUE) 
# save NT plot if desired 
#setwd(path.figs)
#savePlot(filename = "mod9fig01.pdf", type = "pdf")
date() # stop time stamp

date() # start time stamp
# basic BRT model - LR adjusted up (ie faster)
#   CPU ~20 min runtime w/LR=0.01
#   CPU ~5 min runtime w/LR=0.1
mod2.BRT <- gbm.step(data = dat3, gbm.x = pred, gbm.y = 1, family = "bernoulli", 
                     tree.complexity = 10, learning.rate = 0.001, bag.fraction = 0.75, n.folds = 5, 
                     plot.main = TRUE, keep.fold.fit = TRUE,n.trees = 100)
date() # stop time stop

# +
# save NT plot if desired 
# setwd(path.figs)
# savePlot(filename = "mod9fig02.pdf", type = "pdf")
# -

# examine BRT output
ls(mod2.BRT) # examine BRT objects
head(mod2.BRT$fitted) # model fit values
head(dat1$spp65) # observed values
mod2.BRT$contributions # relative variable importance

# examine response:predictor plots
par(mfrow = c(3, 4))
gbm.plot(mod2.BRT, n.plots = 10) # response:predictor plots 
par(mfrow = c(1, 1))
# save plot if desired 
#setwd(path.figs)
#savePlot(filename = "mod9fig03.pdf", type = "pdf")

# search for & examine interactions
mod2.int <- gbm.interactions(mod2.BRT) # examine pairwise interactions
mod2.int$rank.list # matrix of 5 top interactions 
#mod2.int$interactions # NOT RUN: matrix all pairwise interactions

# plot 3 top pairwise interactions
par(mfrow = c(1, 3))
gbm.perspec(mod2.BRT, mod2.int$rank.list[1, 1], mod2.int$rank.list[1, 3], theta = 30)
gbm.perspec(mod2.BRT, mod2.int$rank.list[2, 1], mod2.int$rank.list[2, 3], theta = 30)
gbm.perspec(mod2.BRT, mod2.int$rank.list[3, 1], mod2.int$rank.list[3, 3], theta = 30) 

# +
# save plot if desired 
# setwd(path.figs)
# savePlot(filename = "mod9fig04.pdf", type = "pdf")
# ####### END BRT MODELS
###############################################################################
# -


################################################################################
######## START ACCURACY CALCULATIONS, MODEL=BRT
# build testing dataframe using model predictions
modl <- "mod2.BRT" # add var to keep track of model
dat2 <- cbind(modl, dat1[1], mod2.BRT$fitted, mod2.BRT$fold.fit) # build dataframe
names(dat2)[3:4] <- c("pred", "cvpred") # rename vars
head(dat2, 2) # just to see logit scale
dat2$cvpred <- exp(dat2$cvpred)/(1 + exp(dat2$cvpred)) # convert from logit
head(dat2, 2) # examine prediction dataframe 


# determine best threshold using PresenceAbsence package
#   see help(optimal.thresholds) for more info
#library(PresenceAbsence)  # PresenceAbsence for accuracy metrics
mod.cut <- optimal.thresholds(dat2, opt.methods = c("ObsPrev")) # threshold=PREVALENCE
mod.cut # examine
#mod.cut2 <- optimal.thresholds(dat2, opt.methods = c("MaxKappa")) # MaxKappa option
#mod.cut # MaxKappa threshold
modF.cut$BRT.cut <- mod.cut

# generate confusion matrix
mod2.cfmatR <- table(dat2[[2]], factor(as.numeric(dat2$pred >= mod.cut$pred)))
mod2.cfmatX <- table(dat2[[2]], factor(as.numeric(dat2$cvpred >= mod.cut$cvpred)))
mod2.cfmatR # examine
mod2.cfmatX # examine 


# calculate model accuracies with standard deviation=F
mod2.acc <- presence.absence.accuracy(dat2, threshold = mod.cut$pred, st.dev = F)
tss <- mod2.acc$sensitivity + mod2.acc$specificity - 1 # code TSS metric
mod2.acc <- cbind(mod2.acc[1:7], tss) # bind all metrics
mod2.acc[c(1, 4:5, 7:8)] # examine accuracies

# plotting AUC
auc.roc.plot(dat2, color = T) # basic AUC plot; pkg PresenceAbsence 

brt.predvals <- dat2 %>%
  mutate(presabs.BRT = if_else(pred >= mod.cut$pred,1,0))


# +
# save plot if desired 
# setwd(path.figs)
# savePlot(filename="mod9fig05.tiff",type="tiff")
####### END ACCURACY CALCULATIONS, MODEL=RF
###############################################################################
# -

brt.predvals

# +
###############################################################################
####### START BRT SPATIAL PREDICTION

# mod2.BRTprob  =predict(pred.dom, mod2.BRT, 
#  n.trees = mod2.BRT$gbm.call$best.trees, type = "response", 
#  filename = "mod2.BRTprob.img")

# +

# mod2.BRTclas = reclassify(mod2.BRTprob, c(0, mod.cut[[2]], 0, 
#  mod.cut[[2]], 1, 1))
# writeRaster(mod2.BRTclas, filename = "mod2.BRTclas.img", format = "HFA")
####### END BRT SPATIAL PREDICTION
###############################################################################

# +
# BRT prediction and classification
# -


modFprob.BRT <- predict(mod1.BRT, stream_sample2, n.trees = 100, 
                        type = "response", filename = "modFprob.BRT.img", overwrite = T) # prob map

gbm.step(data = dat1, gbm.x = pred, gbm.y = 1, family = "bernoulli", 
         tree.complexity = 10, learning.rate = 0.0001, bag.fraction = 0.75, n.folds = 5, 
         plot.main = TRUE, keep.fold.fit = TRUE,n.trees = 100)  

save(mod1.BRT, file="brt.Rdata")

saveRDS(modF.cut,"modF.cut.rds")
save(modF.cut, file="modF.cut..RData")

slsc_layer$brt.prob <- modFprob.BRT #add probability values to table


head(slsc_layer,25) #view layer 

# +
##################################################end of script that works correctly
# -
































































# +
# the issue with the code below is:
#slsc_layer$lr.prob <- mod2.pred #add probability values to table

# +
#compile allvalues
slsc_layer <- stream_sample2 #create feauture layer format
slsc_layer$lr.prob <- mod2.pred #add probability values to table
slsc_layer$gam.prob <- modfinl.GAM #add probability values to table
slsc_layer$max.prob <- mod1.MAXprob$mod1.MAXprob #add probability values to table
slsc_layer$rf.prob <- mod1.RF #add probability values to table
slsc_layer$brt.prob <- modFprob.BRT #add probability values to table

head(slsc_layer,25) #view layer 



# +
### 3.2.6 Naïve Bayes (NB)


# +
### Subsection Outputs

# +
## 3.3 Build Ensemble Models

# +

################################################################################
######## START ENSEMBLE PROCESS
# load probability & classified maps; create stacks
#library(raster)
# -

prob.dom <- slsc_layer[,c("lr.prob","gam.prob","max.prob","rf.prob","brt.prob")]
prob.dom # examine
prob.dom <- st_drop_geometry(prob.dom)
# standardize all prediction maps 0-1.0
layers <- {} # initialize (empty) list of raster layers
for (i in 1:length(prob.dom)) {
  m1 <- prob.dom[[i]] # get a prob map
  m2 <- 1/max(prob.dom[[i]]) * prob.dom[[i]] # standardize all probs to max=1
  m2.5 <- names(prob.dom) # split prob layer name apart
  assign(paste(m2.5[i],"STD", sep = "."), m2)
  layers <- cbind(layers, get(paste(m2.5[i],"STD", sep = ".")))
}

# +

m2.5 <-paste(m2.5, "std", sep = ".")   #add abbreviation to names "std"
colnames(layers) <- m2.5
# -


colnames(layers) <- c("lr.prob.std","gam.prob.std","max.prob.std","rf.prob.std","brt.prob.std")
#m2.5  #assign column names to scaled prob values
slsc_layer <-cbind(slsc_layer,layers)   #bind scaled probs to feature layer


slsc_layer$Ensemble_ave <-  ((slsc_layer$lr.prob.std+
                                slsc_layer$gam.prob.std+
                                slsc_layer$max.prob.std+
                                slsc_layer$rf.prob.std+
                                slsc_layer$brt.prob.std)/5)


slsc_layer$lr.class <- ifelse(slsc_layer$lr.prob >modF.cut$LR.cut$mod2.pred, 1, 0)
modF.cut$LR.cut$mod2.pred


slsc_layer$gam.class <- ifelse(slsc_layer$gam.prob >modF.cut$gam.cut$modfinl.pred, 1, 0)
modF.cut$gam.cut$modfinl.pred

slsc_layer$max.class <- ifelse(slsc_layer$max.prob >modF.cut$max.cut, 1, 0)
modF.cut$max.cut

slsc_layer$rf.class <- ifelse(slsc_layer$rf.prob >modF.cut$rf.cut$mod1.pred, 1, 0)

slsc_layer$brt.class <- ifelse(slsc_layer$brt.prob >modF.cut$BRT.cut$cvpred, 1, 0)



slsc_layer <- slsc_layer %>%
  mutate(concordance = lr.class+gam.class + max.class + rf.class + brt.class)



#Write out Shp for stream_sample3 and stream_sample3.ng
library(sf)

# +
st_write(pts_snap.mcpsf,  "mcp_polygon.shp", 
         delete_dsn=FALSE, 
         update = TRUE )
st_write(stream_sample3$geometry,  "stream_sample3.shp", 
         delete_dsn=FALSE, 
         update = TRUE, layer_options = "OVERWRITE=yes" )
st_write(slsc_layer,  "ybcfinal.shp", delete_dsn=FALSE, update = TRUE, , layer_options = "OVERWRITE=true")
st_write(tryme3[,"presabs"],  "presabs.shp", 
         delete_dsn=FALSE, 
         update = TRUE )
st_write(stream_sample_huc,  "stream_sample_huc.shp", 
         delete_dsn=FALSE, 
         update = TRUE )
st_write(pts_snap_sf,  "presence_points.shp", 
         delete_dsn=FALSE, 
         update = TRUE )


library(mapview)
# -

library("mapview")
mapview(slsc_layer, 
        zcol = c("Ensemble_ave"), map.types = c("Esri.WorldImagery", "OpenTopoMap"))

mapview(slsc_layer, zcol = "concordance", 
        map.types = c("Esri.WorldImagery", "OpenTopoMap")) 

mapview(slsc_layer, zcol = c("lr.class", "gam.class","max.class", "rf.class", "brt.class"), map.types = c("Esri.WorldImagery", "OpenTopoMap"))+pts_snap_sf$geometry


# +
### Subsection Outputs


# +
#  4 SDHM IMPLEMENTATION

# +
## 4.2 DWR and cooperator general management and conservation goals
### 4.2.1 DWR-initiated; clueless here

# +
## 4.3 Public outreach
### 4.3.1 DWR-initiated; clueless here



