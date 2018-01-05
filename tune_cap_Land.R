
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
DEM<-raster("~/WP2/data/topo_ProjectRaster3.tif.ovr")

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


reclass <- function(df, i, j) {
  df[, "DON"][df[, "DON"] <= i] <- "Low"
  df[, "DON"][df[, "DON"] < j] <- "Medium"
  df[, "DON"][(df[, "DON"] != "Low") & (df[, "DON"] != "Medium")] <- "High"
  df[, "DON"] <- factor(df[, "DON"], levels = c("Low", "Medium", "High"))
  return(df)
}


reclass4<-function(df,i,j){
  for (t in c(1,2)){
    df[, t][df[, t] <=i] <- "Low"
    df[, t][df[, t] < j] <- "Medium"
    df[, t][(df[, t] != "Low") & (df[, t] != "Medium")] <- "High"
    df[, t] <- factor(df[, t], levels = c("Low", "Medium", "High"))
  }
  return(df)
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
    slope=mean(aa[,9])
    aspect=mean(aa[,10])
    sing_land<-data.frame(Soil,Veg,Landuse,Catchment,GW_depth,Distance,Distance_LP,Distance_GWC,slope,aspect)
    landscape_all<-rbind(landscape_all,sing_land)
  }
  return(landscape_all)
}

## preprocess the landscape raster
Soil <- pre(Soil)
Veg <- pre(Veg)
Land_use <- pre(Land_use)
Cat <- pre(Cat)
DEM <- pre(DEM)


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

x<-terrain(DEM,opt=c("slope","aspect"),neighbors = 8,unit = "degrees")
slope=x[[1]]
aspect=x[[2]]

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
landscapes<-stack(Soil,Veg,Land_use,Cat,depth_k,water_distance,distance_LP,Distance_GWC,slope,aspect)
names(landscapes) <- c("Soil", "Veg", "Landuse","Catchment", "GW_depth", "Distance","Distance_LP","Distance_GWC","slope","aspect")

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

newdata <- cbind(all_points2, all_points2[n1.d,"DON"],all_points2[n2.d,"DON"],all_points2[n3.d,"DON"])
newdata$DON_m3<-(newdata$DON.1+newdata$DON.2+newdata$DON.3)/3

newdata$dev<-abs(newdata$DON-newdata$DON_m3)/newdata$DON_m3

newdata[newdata$dev<=5,"type"]=1
newdata[newdata$dev>5,"type"]=0

all_points<-data.frame(newdata)
all_points<-subset(all_points,all_points$type==1)

## set the parameters for mlr
seed=35
set.seed(seed)
reg_rf = makeLearner("regr.randomForest")
class_rf = makeLearner("classif.randomForest")

#class_rf$par.vals<-list(importance=T)
ctrl = makeTuneControlIrace(maxExperiments = 500L)
rdesc = makeResampleDesc("CV", iters = 5)

## define the parameter spaces for RF      
para_rf = makeParamSet(
  makeDiscreteParam("ntree", values=seq(200,500,50)),
  makeIntegerParam("nodesize", lower = 50, upper = 55),
  makeIntegerParam("mtry", lower = 2, upper =3)
  #  makeDiscreteParam("coefReg", values=seq(0.05,0.2,0.05))
)

model_build <- function(dataset, n_target) {
  #set.seed(719)
  ## define the regression task for DON 
  WP3_target = makeRegrTask(id = "WP3_target", data = dataset, target = n_target)
  rin = makeResampleInstance(rdesc, task = WP3_target)
  res_rf = mlr::tuneParams(reg_rf, WP3_target, resampling = rdesc, par.set = para_rf, control = ctrl,
                           show.info = FALSE)
  lrn_rf = setHyperPars(reg_rf, par.vals = res_rf$x)
  ## train the final model 
  #set.seed(719)
  rf <- mlr::train(lrn_rf, WP3_target)
  return(rf)
}


model_build2 <- function(dataset, n_target) {
  #set.seed(719)
  ## define the regression task for DON 
  WP3_target = makeClassifTask(id = "WP3_target", data = dataset, target = n_target)
  rin = makeResampleInstance(rdesc, task = WP3_target)
  res_rf = mlr::tuneParams(class_rf, WP3_target, resampling = rdesc, par.set = para_rf, control = ctrl,
                           show.info = FALSE)
  lrn_rf = setHyperPars(class_rf, par.vals = res_rf$x)
  ## train the final model 
  #set.seed(719)
  rf <- mlr::train(lrn_rf, WP3_target)
  return(rf)
}

a1=1.0
a2=2.0
all_results<-data.frame()
  training <- all_points
  ## load the point data 
  training_df <- read_pointDataframes(training)  
  training_points<- read_points(training)
  
  ## map1, using kringing for DON interpolation
  f.1 <- as.formula(log10(DON) ~ 1)
  # Add X and Y to training 
  training_df<-add_S1S2(training_df)  
  ## M2, using RF to predict the DON
 for (a in seq(100,1500,100)){
    for (b in seq(100,2000,100)){

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
  M2_train <- cbind(as.data.frame(landscape_train), training_df@data[c("DON","s1","s2")])

   M2_train<-reclass(M2_train,a1,a2)
#  M2_train$DON<-log10(M2_train$DON)
#  M2_test$DON<-log10(M2_test$DON)
  
  for(i in c("GW_depth","Distance","Distance_GWC","slope","aspect","s1","s2")){
    
    min_train<-min(M2_train[,i])
    max_train<-max(M2_train[,i])
    
    M2_train[,i]<-(M2_train[,i]-min_train)/(max_train-min_train)
    
    sd_train<-sd(M2_train[,i])
    mean_train<-mean(M2_train[,i])
    
    M2_train[,i]<-(M2_train[,i]-mean_train)/sd_train
    
  }
  
  WP2Train<-M2_train[,-c(4,7,9,10)]
  
  rf_DON_m2 <- model_build2(WP2Train,"DON")
  map2_train<-predict(rf_DON_m2,newdata=WP2Train)
  map2_train_cla <- data.frame(observed_DON=map2_train$data$truth,predicted_DON=map2_train$data$response)
  #map2_predict_cla<-reclass4(map2_predict_cla,a1,a2
  M2_ACC_train<-postResample(map2_train_cla[,2],map2_train_cla[,1])[1]
  ## create the training and testing sets 
  sing_acc<-data.frame(a1,a2,M2_ACC_train)
  all_results<-rbind(all_results,sing_acc)
  print(all_results)
 
}}
