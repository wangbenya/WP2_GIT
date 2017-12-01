
## load the library 
rm(list=ls())

library(raster)
library(sp)
library(randomForest)
library(magrittr)
library(rgdal)
library(gstat)
library(ggplot2)
library(mlr)
library(SemiPar)
library(Hmisc)
library(foreign)
library(maptools)
library(prettymapr)
library(mlrMBO)
library(parallelMap)
library(caret)
library(automap)
library(reshape2)

## start the parallel 
parallelStartSocket(16)

WGS84 <- CRS("+proj=utm +zone=50 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

study_area <- shapefile("~/WP2/data/study_area.shp")
water <- shapefile("~/WP2/data/water.shp")

## load the veg, soil, land use, groundwater subarea,
## surface water subare, catchment
Soil <- raster("~/WP2/data/soil1.ovr")
Veg <- raster("~/WP2/data/vegetation.ovr")
Land_use <- raster("~/WP2/data/landuse1.ovr")
Cat <- raster("~/WP2/data/catch_name2.tif.ovr")
## define the function 
## preprocess 
study_area <- spTransform(study_area, WGS84)
extent <- c(study_area@bbox[1, 1:2], study_area@bbox[2, 1:2])

water <- spTransform(water, WGS84)

pre <- function(x) {
  projection(x) <- WGS84
  extent(x) <- extent
  x <- raster::mask(x, study_area)
  return(x)
}

read_points <- function(read_data) {
  SP <- SpatialPoints(read_data[, 2:3], proj4string = CRS("+proj=longlat +ellps=GRS80 +no_defs"))
  SP <- spTransform(SP, WGS84)
  SP@bbox <- study_area@bbox
  if (length(zerodist(SP)) >= 1) {
    SP <- SP[-(zerodist(SP)[, 1]),]
  }
  #   plot(study_area_withW)
  #  points(SP@coords)
  return(SP)
}

read_pointDataframes <- function(read_data) {
  SP <- SpatialPoints(read_data[, 2:3], proj4string = CRS("+proj=longlat +ellps=GRS80 +no_defs"))
  SPD <- SpatialPointsDataFrame(SP, read_data)
  SPD <- spTransform(SPD, WGS84)
  SPD@bbox <- study_area@bbox
  if (length(zerodist(SPD)) >= 1) {
    SPD <- SPD[-(zerodist(SPD)[, 1]),]
  }
  # plot(study_area_withW)
  #points(SPD@coords)
  return(SPD)
}

# Add X and Y to training 
add_S1S2 <- function(dataset) {
  dataset$s1 <- coordinates(dataset)[, 1]
  dataset$s2 <- coordinates(dataset)[, 2]
  return(dataset)
}

get_landscape<-function(df){
  landscape_all<-data.frame()
  for (ii in seq(1,length(df))){
    aa<-as.data.frame(df[[ii]])
    aa<-subset(aa,aa$Soil!="NA")
    Soil=tail(names(sort(table(aa[,1]))),1)
    Veg=tail(names(sort(table(aa[,2]))),1)
    Landuse=tail(names(sort(table(aa[,3]))),1)
    Catchment=tail(names(sort(table(aa[,4]))),1)
    GW_depth=mean(aa[,5])
    Distance=mean(aa[,6])
    Distance_LP=mean(aa[,7])
    Distance_GWC=mean(aa[,8])
    sing_land<-data.frame(Soil,Veg,Landuse,Catchment,GW_depth,Distance,Distance_LP,Distance_GWC)
    landscape_all<-rbind(landscape_all,sing_land)
  }
  return(landscape_all)
}

## preprocess the landscape raster
Soil <- pre(Soil)
Veg <- pre(Veg)
Land_use <- pre(Land_use)
Cat <- pre(Cat)

v_Veg<-values(Veg)
v_Veg[v_Veg %in% c(2,3,4)]=1
v_Veg[v_Veg %in% c(8,9)]=8
v_Veg[v_Veg %in% c(12,13)]=12
v_Veg[v_Veg %in% c(18,19,20)]=18
values(Veg)<-v_Veg

v_land<-values(Land_use)
v_land[v_land %in% c(1,2,13)]=1
v_land[v_land %in% c(5,7,12,6,11)]=5
v_land[v_land %in% c(8,10)]=8
values(Land_use)<-v_land

v_soil<-values(Soil)
v_soil[v_soil %in% c(1,2)]=1
v_soil[v_soil %in% c(4,5)]=4
v_soil[v_soil %in% c(6,7)]=6
v_soil[v_soil %in% c(11,12)]=11
v_soil[v_soil %in% c(13,14)]=13
values(Soil)<-v_soil

# Create an empty grid where n is the total number of cells
r <- raster(study_area)
res(r) <- res(Soil) # 10 km if your CRS's units are in km
base_grid <- as(r, 'SpatialGrid')
#plot(base_grid)

## M2, using RF to predict the DON
depth <- read.csv("~/WP2/data/sampling_depth.csv",header=T) %>% read_pointDataframes(.)

# Define the 1st order polynomial equation
f_depth <- as.formula(sampling_d ~ 1)
# Add X and Y to training 
depth<-add_S1S2(depth)
# variogram on the de-trended data.
var.depth <- variogram(f_depth, depth)
#plot(var.depth)
dat.fit_depth <- fit.variogram(var.depth,vgm(c("Sph","Exp")))
# created in the earlier step)
depth_k <- krige(f_depth, depth, base_grid, dat.fit_depth) %>% raster(.) %>% raster::mask(., study_area)
#plot(depth_k)
depth_k@data@names<-"GW_depth"

#Now make the map
### distance
water <- raster::rasterize(water, depth_k)
water_distance <- raster::mask(distance(water),study_area)
water_distance@data@names<-"Distance_to_water"


left_up<-water 
values(left_up)<-NA
values(left_up)[[1]]<-1
distance_LP <-raster::mask(distance(left_up),study_area)
distance_LP@data@names<-"Distance_LP"


GW_center<-data.frame(Latitude=c(6495000,6475000,6460000,6448000,6403000),Longitude=rep(402000,5),values=1)
GW_center <- SpatialPoints(GW_center[, c(2:1)], proj4string = WGS84)
GW_center@bbox <- study_area@bbox
base_GWC<-water 
values(base_GWC)<-1
Distance_GWC<-distanceFromPoints(base_GWC,GW_center)
Distance_GWC@data@names<-"Distance_GWC"

## load the data 
landscapes<-stack(Soil,Veg,Land_use,Cat,depth_k,water_distance,distance_LP,Distance_GWC)
names(landscapes) <- c("Soil", "Veg", "Landuse","Catchment", "GW_depth", "Distance","Distance_LP","Distance_GWC")

## load the data 
set.seed(666)
seed.list<-sample(1:1000,300,replace =F)


all_points<-read.csv("~/WP2/data/all_data1127.csv",header = T)
extra_n<-read.csv("~/WP2/data/extra_n.csv",header = T)
extra_n<-subset(extra_n,!(extra_n$WIN_Site_ID %in% all_points$WIN_Site_ID))
DON_GW4<-read.csv("~/WP2_GIT/DON_GW4.csv",header = T)
DOC_GW4<-read.csv("~/WP2_GIT/DOC_GW4.csv",header = T)
NOx_GW4<-read.csv("~/WP2_GIT/NOx_GW4.csv",header = T)
NH4_GW4<-read.csv("~/WP2_GIT/NH4_GW4.csv",header = T)
TN_GW4<-read.csv("~/WP2_GIT/TN_GW4.csv",header = T)
extra_n<-subset(extra_n,!(extra_n$WIN_Site_ID %in% all_points$WIN_Site_ID))

#Make a distance matrix
all_points2<-read_pointDataframes(all_points)
all_points2<-add_S1S2(all_points2)
a=data.frame(all_points2)[,c('s1','s2')]
d <- pointDistance(a, lonlat=F)
n1.d <- apply(d, 1, function(x) order(x, decreasing=F)[2])
n2.d <- apply(d, 1, function(x) order(x, decreasing=F)[3])
n3.d <- apply(d, 1, function(x) order(x, decreasing=F)[4])
n4.d <- apply(d, 1, function(x) order(x, decreasing=F)[5])
n5.d <- apply(d, 1, function(x) order(x, decreasing=F)[6])

newdata <- cbind(all_points2, all_points2[n1.d,"DON"],all_points2[n2.d,"DON"],all_points2[n3.d,"DON"],all_points2[n4.d,"DON"],all_points2[n5.d,"DON"])
newdata$DON_m3<-(newdata$DON.1+newdata$DON.2+newdata$DON.3)/3

newdata$variance<-(newdata$DON-newdata$DON_m3)^2
newdata$dev<-abs(newdata$DON-newdata$DON_m3)/newdata$DON_m3

newdata[newdata$dev<=1,"type"]=1
newdata[newdata$dev>1,"type"]=0

all_points<-data.frame(newdata)

all_points<-subset(all_points,all_points$type==1)

## set the parameters for mlr
seed=35
set.seed(seed)
reg_rf = makeLearner("regr.randomForest")
#class_rf$par.vals<-list(importance=T)
ctrl = makeTuneControlIrace(maxExperiments = 200L)
rdesc = makeResampleDesc("CV", iters = 5)

## define the parameter spaces for RF      
para_rf = makeParamSet(
  makeDiscreteParam("ntree", values=seq(50,200,20)),
  makeIntegerParam("nodesize", lower = 10, upper = 15),
  makeIntegerParam("mtry", lower = 2, upper =6)
  #  makeDiscreteParam("coefReg", values=seq(0.05,0.2,0.05))
)

model_build <- function(dataset, n_target) {
  #set.seed(719)
  ## define the regression task for DON 
  WP3_target = makeRegrTask(id = "WP3_target", data = dataset, target = n_target)
  ## cross validation
  ## 10-fold cross-validation
  rin = makeResampleInstance(rdesc, task = WP3_target)
  res_rf = mlr::tuneParams(reg_rf, WP3_target, resampling = rdesc, par.set = para_rf, control = ctrl,
                           show.info = FALSE,measures=rsq)
  lrn_rf = setHyperPars(reg_rf, par.vals = res_rf$x)
  
  ## train the final model 
  #set.seed(719)
  rf <- mlr::train(lrn_rf, WP3_target)
  return(rf)
}

all_results<-data.frame()

for (tt in c(1)){
  print(tt)
  seeds<-seed.list[tt]
  set.seed(seeds)
  trainIndex <- createDataPartition(all_points$DON, p = 0.8, list = FALSE)
  
  training <- all_points[trainIndex,]
  testing <- all_points[-trainIndex,]
  
  ## load the point data 
  training_df <- read_pointDataframes(training)
  testing_df <-  read_pointDataframes(testing) 
  
  training_points<- read_points(training)
  testing_points <- read_points(testing)
  
  ## map1, using kringing for DON interpolation
  f.1 <- as.formula(log10(DON) ~ 1)
  # Add X and Y to training 
  training_df<-add_S1S2(training_df)
  testing_df<-add_S1S2(testing_df)
  
  # Compute the sample variogram; note that the f.1 trend model is one of the
  var.smpl1 <- variogram(f.1, training_df)
  plot(var.smpl1)
  # Compute the variogram model by passing the nugget, sill and range value
  dat.fit1 <- fit.variogram(var.smpl1,vgm(c("Sph")))
  
  plot(var.smpl1,dat.fit1)
  # Perform the krige interpolation (note the use of the variogram model
  kriging_DON_m1 <- krige(f.1, training_df, base_grid, dat.fit1) %>% raster(.) %>% raster::mask(., study_area)
  values(kriging_DON_m1) <- 10 ^ (values(kriging_DON_m1))
  dat.krg_DON<-kriging_DON_m1
  
  map1_predict <- data.frame(observed_DON=testing_df@data$DON,predicted_DON=raster::extract(kriging_DON_m1, testing_points))
  
  M1_rmse<-postResample(map1_predict[,2],map1_predict[,1])[1]
  M1_r2<-postResample(map1_predict[,2],map1_predict[,1])[2]
  
  map1_train <- data.frame(observed_DON=training_df@data$DON,predicted_DON=raster::extract(kriging_DON_m1, training_points))
  
  M1_rmse_train<-postResample(map1_train[,2],map1_train[,1])[1]
  M1_r2_train<-postResample(map1_train[,2],map1_train[,1])[2]
  
  ## M2, using RF to predict the DON
  for (a in seq(100,1500,100)){
    for (b in seq(100,3000,100)){
  capture_zone_land<-function(df){
    num<-nrow(df)
    landscape_data<-data.frame()
    for (r in seq(1,num)){
      p1_long<-df@coords[r,1]
      p1_lat<-df@coords[r,2]
      pg<-spPolygons(rbind(c(p1_long,p1_lat),c(p1_long+a,p1_lat+b),c(p1_long+2*a,p1_lat+b),
                           c(p1_long+2*a,p1_lat-b),c(p1_long+a,p1_lat-b),c(p1_long,p1_lat)))  
      projection(pg)<- WGS84
      p1_landscape<-raster::extract(landscapes,pg)
      p1_landscape<-get_landscape(p1_landscape)
      landscape_data<-rbind(landscape_data,p1_landscape)
    }
    return(landscape_data)
  }
  
  landscape_train <- capture_zone_land(training_df)
  landscape_test <- capture_zone_land(testing_df)
  
  M2_train <- cbind(as.data.frame(landscape_train), training_df@data[c("DON","Collect_Month","date_","s1","s2")])
  M2_test <- cbind(as.data.frame(landscape_test), testing_df@data[c("DON","Collect_Month","date_","s1","s2")])
  
  names(M2_train) <- colnames(M2_test)
  
  common_landscape<-function(land){
    land_dataset<-data.frame(table(M2_train[,land]))
    land_common<-subset(land_dataset,land_dataset[,2]==max(land_dataset[,2]))[1]
    return(as.matrix(land_common))
  }
  
  soil_max = common_landscape("Soil")[1]
  veg_max=common_landscape("Veg")[1]
  landuse_max = common_landscape("Landuse")[1]
  cat_max = common_landscape("Catchment")[1]
  
  max_list<-list(soil_max,veg_max,landuse_max,cat_max)
  
  for (ii in seq(1,4)){
    M2_train[,ii]<-factor(M2_train[,ii],levels = unique(values(landscapes[[ii]]))[-1])
    M2_test[,ii]<-factor(M2_test[,ii],levels=unique(values(landscapes[[ii]]))[-1])
    M2_test [(which(!(M2_test[,ii] %in% M2_train[,ii]))),ii]<-as.numeric(max_list[[ii]])
    
    M2_train[,ii]<-droplevels(M2_train[,ii])
    M2_test[,ii]<-factor(M2_test[,ii],levels = levels(M2_train[,ii]))
    
  }
  
  ## build the model for map2
  #names(M2_train)<-c("Soil", "Veg", "Landuse","Catchment", "GW_depth", "Distance", "DON","s1","s2")
  #names(M2_test)<-c("Soil",  "Veg", "Landuse","Catchment", "GW_depth", "Distance", "DON","s1","s2")
  
  WP2Train<-M2_train
  WP2Test<-M2_test
  
  WP2Train$date_<-(WP2Train$date_)^2
  WP2Test$date_<-(WP2Test$date_)^2
  
  WP2Train$Distance_GWC<-(WP2Train$Distance_GWC)^2
  WP2Test$Distance_GWC<-(WP2Test$Distance_GWC)^2
  
  WP2Train$GW_depth<-(WP2Train$GW_depth)^2
  WP2Test$GW_depth<-(WP2Test$GW_depth)^2
  
  WP2Test$Distance<-(WP2Test$Distance)^2
  WP2Test$Distance<-(WP2Test$Distance)^2
  
  #WP2Train$Latitude<--WP2Train$Latitude
  #M2_test$Latitude<--M2_test$Latitude
  
  #WP2Train$Distance<-log10(WP2Train$Distance+0.01)
  #WP2Test$Distance<-log10(WP2Test$Distance+0.01)
  
  for(i in c(5:8,11:13)){
    
    min_train<-min(WP2Train[,i])
    max_train<-max(WP2Train[,i])
    
    WP2Train[,i]<-(WP2Train[,i]-min_train)/(max_train-min_train)
    WP2Test[,i]<-(WP2Test[,i]-min_train)/(max_train-min_train)
    
    sd_train<-sd(WP2Train[,i])
    mean_train<-mean(WP2Train[,i])
    
    WP2Train[,i]<-(WP2Train[,i]-mean_train)/sd_train
    WP2Test[,i]<-(WP2Test[,i]-mean_train)/sd_train
    
  }
  
  set.seed(seeds)
  WP2Train<-WP2Train[,-c(7,12,13)]
  WP2Test<-WP2Test[,-c(7,12,13)]
  
  rf_DON_m2 <- model_build(WP2Train,"DON")
  
  map2_predict<-predict(rf_DON_m2,newdata=WP2Test)
  map2_train<-predict(rf_DON_m2,newdata=WP2Train)
  
  M2_rmse<-postResample(map2_predict$data$response, map2_predict$data$truth)[1]
  M2_r2<-postResample(map2_predict$data$response, map2_predict$data$truth)[2]
  
  M2_rmse_train<-postResample(map2_train$data$response, map2_train$data$truth)[1]
  M2_r2_train<-postResample(map2_train$data$response, map2_train$data$truth)[2]
  
# create the training and testing sets
  ## map3 predict accuracy
  sing_acc<-data.frame(a,b,M2_r2,M2_r2_train)
  all_results<-rbind(all_results,sing_acc)
  
  print(all_results)
}
}
}


