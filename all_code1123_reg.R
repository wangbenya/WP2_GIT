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
    sing_land<-data.frame(Soil,Veg,Landuse,Catchment,GW_depth,Distance)
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
# Run the regression model
lm.depth <- lm(f_depth, data = depth)
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
## load the data 
landscapes<-stack(Soil,Veg,Land_use,Cat,depth_k,water_distance)
names(landscapes) <- c("Soil", "Veg", "Landuse","Catchment", "GW_depth", "Distance")

## load the data 
set.seed(666)
seed.list<-sample(1:1000,50,replace =F)

all_points<-read.csv("~/WP2/data/all_data1210.csv",header = T)

extra_n<-read.csv("~/WP2/data/extra_n.csv",header = T)
extra_n<-subset(extra_n,!(extra_n$WIN_Site_ID %in% all_points$WIN_Site_ID))
DON_GW4<-read.csv("~/WP2_GIT/DON_GW4.csv",header = T)
DOC_GW4<-read.csv("~/WP2_GIT/DOC_GW4.csv",header = T)
NOx_GW4<-read.csv("~/WP2_GIT/NOx_GW4.csv",header = T)
NH4_GW4<-read.csv("~/WP2_GIT/NH4_GW4.csv",header = T)
TN_GW4<-read.csv("~/WP2_GIT/TN_GW4.csv",header = T)


      ## set the parameters for mlr
      seed=35
      set.seed(seed)
      reg_rf = makeLearner("regr.RRF")
      #class_rf$par.vals<-list(importance=T)
      ctrl = makeTuneControlIrace(maxExperiments = 200L)
      rdesc = makeResampleDesc("CV", iters = 5)
      
      ## define the parameter spaces for RF      
      para_rf = makeParamSet(
        makeDiscreteParam("ntree", values=seq(100,200,20)),
        makeIntegerParam("nodesize", lower = 4, upper = 10),
        makeIntegerParam("mtry", lower = 4, upper =12),
        makeDiscreteParam("coefReg", values=seq(0.1,0.3,0.05))
      )

      model_build <- function(dataset, n_target) {
        set.seed(719)
          ## define the regression task for DON 
          WP3_target = makeRegrTask(id = "WP3_target", data = dataset, target = n_target)
          ## cross validation
          ## 10-fold cross-validation
          rin = makeResampleInstance(rdesc, task = WP3_target)
          res_rf = mlr::tuneParams(reg_rf, WP3_target, resampling = rdesc, par.set = para_rf, control = ctrl,
                                   show.info = FALSE,measures=rsq)
          lrn_rf = setHyperPars(reg_rf, par.vals = res_rf$x)
        
        ## train the final model 
        set.seed(719)
        rf <- mlr::train(lrn_rf, WP3_target)
        return(rf)
      }

all_results<-data.frame()

 for (tt in seq(1,5)){
      print(tt)
      seeds<-seed.list[tt]
      set.seed(seeds)
      trainIndex <- createDataPartition(all_points$DON, p = .8, list = FALSE)
  
      training <- all_points[trainIndex,]
      testing <- all_points[-trainIndex,]
  
  ## load the point data 
  training_df <- training[,c(1,2,3,5)] %>% read_pointDataframes(.)
  testing_df <-  read_pointDataframes(testing) 
  
  training_points<-training[,c(1,2,3,5)] %>% read_points(.)
  
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
  a=150
  b=750
  capture_zone_land<-function(df){
  num<-nrow(df)
  landscape_data<-data.frame()
  for (r in seq(1,num)){
  p1_long<-df@coords[r,1]
  p1_lat<-df@coords[r,2]
  pg<-spPolygons(rbind(c(p1_long,p1_lat),c(p1_long,p1_lat+a),c(p1_long+sqrt(2)/2*b,p1_lat+a+sqrt(2)/2*b),
                       c(p1_long+a+sqrt(2)/2*b,p1_lat+sqrt(2)/2*b),c(p1_long+a,p1_lat),c(p1_long,p1_lat)))  
  projection(pg)<- WGS84
  p1_landscape<-raster::extract(landscapes,pg)
  p1_landscape<-get_landscape(p1_landscape)
  landscape_data<-rbind(landscape_data,p1_landscape)
  }
  return(landscape_data)
  }
  
  landscape_train <- capture_zone_land(training_df)
  landscape_test <- capture_zone_land(testing_df)
  
  M2_train <- cbind(as.data.frame(landscape_train), training_df@data[c("DON","Longitude","Latitude")])
  M2_test <- cbind(as.data.frame(landscape_test), testing_df@data[c("DON","Longitude","Latitude")])
  
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
      M2_train[,ii]<-as.factor(M2_train[,ii])
      M2_test [(which(!(M2_test[,ii] %in% M2_train[,ii]))),ii]<-as.numeric(max_list[[ii]])
      M2_test[,ii]<-factor(M2_test[,ii],levels=levels(M2_train[,ii]))
  }
  
  ## build the model for map2
  names(M2_train)<-c("Soil", "Veg", "Landuse","Catchment", "GW_depth", "Distance", "DON","Longitude","Latitude")
  names(M2_test)<-c("Soil",  "Veg", "Landuse","Catchment", "GW_depth", "Distance", "DON","Longitude","Latitude")
  
  WP2Train<-M2_train
  WP2Test<-M2_test
  
  WP2Train$Distance<-log10(WP2Train$Distance+0.01)
  
  WP2Test$Distance<-log10(WP2Test$Distance+0.01)
  
  for(i in c(5,6,8,9)){

    min_train<-min(WP2Train[,i])
    max_train<-max(WP2Train[,i])
    
    WP2Train[,i]<-(WP2Train[,i]-min_train)/(max_train-min_train)
    WP2Test[,i]<-(WP2Test[,i]-min_train)/(max_train-min_train)
    
    sd_train<-sd(WP2Train[,i])
    mean_train<-mean(WP2Train[,i])
    
    WP2Train[,i]<-(WP2Train[,i]-mean_train)/sd_train
    WP2Test[,i]<-(WP2Test[,i]-mean_train)/sd_train
    
  }
  
  #WP2Train$DON<-log10(WP2Train$DON)
  #WP2Test$DON<-log10(WP2Test$DON)
  
  set.seed(seeds)
  WP2Train<-createDummyFeatures(WP2Train,target = "DON")
  WP2Test<-createDummyFeatures(WP2Test,target = "DON")
  
  rf_DON_m2 <- model_build(WP2Train,"DON")
  
  map2_predict<-predict(rf_DON_m2,newdata=WP2Test)
  map2_train<-predict(rf_DON_m2,newdata=WP2Train)
  
  M2_rmse<-postResample(map2_predict$data$response, map2_predict$data$truth)[1]
  M2_r2<-postResample(map2_predict$data$response, map2_predict$data$truth)[2]
  
  M2_rmse_train<-postResample(map2_train$data$response, map2_train$data$truth)[1]
  M2_r2_train<-postResample(map2_train$data$response, map2_train$data$truth)[2]
  
  ## map4, kriging first and then rf
  # kriging for DOC
  f.DOC <- as.formula(log10(DOC) ~ 1)
  
  training_DOC <- training[,c(1,2,3,4)] %>% rbind(.,extra_n[,c(1,2,3,4)]) %>%
    rbind(.,DOC_GW4) %>% subset(.,.[,"DOC"]!="NA") %>% read_pointDataframes(.)
  
  training_DOC<-add_S1S2(training_DOC)
  var.smpl_DOC <- variogram(f.DOC, training_DOC)
  plot(var.smpl_DOC)
  
  dat.fit_DOC <- fit.variogram(var.smpl_DOC,vgm(c("Sph","Exp")))
  plot(var.smpl_DOC,dat.fit_DOC)
  # Perform the krige interpolation (note the use of the variogram model
  dat.krg_DOC <- krige(f.DOC, training_DOC, base_grid, dat.fit_DOC) %>% raster(.) %>% raster::mask(., study_area)
  values(dat.krg_DOC) <- 10 ^ (values(dat.krg_DOC))
  
  ## create rasterstack with kriging data
  kriging_nutrietn<-stack(dat.krg_DON,dat.krg_DOC)
  names(kriging_nutrietn) <- c("DON_k","DOC_k")
  
  ## extract the data from landscapes_withN

  capture_zone_nutrient<-function(df){
    num<-nrow(df)
    landscape_data<-data.frame()
    for (r in seq(1,num)){
      p1_long<-df@coords[r,1]
      p1_lat<-df@coords[r,2]
      pg<-spPolygons(rbind(c(p1_long,p1_lat),c(p1_long,p1_lat+a),c(p1_long+sqrt(2)/2*b,p1_lat+a+sqrt(2)/2*b),
                       c(p1_long+a+sqrt(2)/2*b,p1_lat+sqrt(2)/2*b),c(p1_long+a,p1_lat),c(p1_long,p1_lat)))  
      projection(pg)<- WGS84
      p1_landscape<-raster::extract(kriging_nutrietn,pg)
      land_data<-data.frame(DON_k=mean(p1_landscape[[1]][,1]),DOC_k=mean(p1_landscape[[1]][,2]))
      landscape_data<-rbind(landscape_data,land_data)
    }
    return(landscape_data)
  }
  
  landscape_train_withKN <- capture_zone_nutrient(training_df)
  landscape_test_withKN <- capture_zone_nutrient(testing_df)
  
  M4_train_withKN <- cbind(M2_train,as.data.frame(landscape_train_withKN))
  M4_test_withKN <- cbind(M2_test,as.data.frame(landscape_test_withKN))
  names(M4_test_withKN) <- names(M4_train_withKN)
  
  ## create the training and testing sets 
  ## build the model for map2
  #names(M4_train_withKN)[1:9]<-c("Soil", "Veg", "Landuse","Catchment", "GW_depth", "Distance", "DON","Longitude","Latitude")
  #names(M4_test_withKN)[1:9]<-c("Soil",  "Veg", "Landuse", "Catchment", "GW_depth", "Distance", "DON","Longitude","Latitude")
  
  #M4_train_withKN <- reclass3(M4_train_withKN,0.5,1.0)
  #M4_test_withKN <- reclass3(M4_test_withKN,0.5,1.0)
  
  M4_train_withKN<-M4_train_withKN[,-c(10)]
  M4_test_withKN<-M4_test_withKN[,-c(10)]
  
  M4_train_withKN$Distance<-log10(M4_train_withKN$Distance+0.01)
  
  M4_test_withKN$Distance<-log10(M4_test_withKN$Distance+0.01)
   
  M4_train_withKN$DOC_dep<-M4_train_withKN$GW_depth*M4_train_withKN$DOC_k
  M4_test_withKN$DOC_dep<-M4_test_withKN$GW_depth*M4_test_withKN$DOC_k

    for(i in c(5,6,8:11)){
    min_train<-min(M4_train_withKN[,i])
    max_train<-max(M4_train_withKN[,i])
    
    M4_train_withKN[,i]<-(M4_train_withKN[,i]-min_train)/(max_train-min_train)
    M4_test_withKN[,i]<-(M4_test_withKN[,i]-min_train)/(max_train-min_train)
    
    sd_train<-sd(M4_train_withKN[,i])
    mean_train<-mean(M4_train_withKN[,i])
    
    M4_train_withKN[,i]<-(M4_train_withKN[,i]-mean_train)/sd_train
    M4_test_withKN[,i]<-(M4_test_withKN[,i]-mean_train)/sd_train
    
      }
  
  set.seed(seeds)
  #M4_train_withKN$DON<-log10(M4_train_withKN$DON)
  #M4_test_withKN$DON<-log10(M4_test_withKN$DON)
  
  M4_train_withKN<-createDummyFeatures(M4_train_withKN,target = "DON")
  M4_test_withKN<-createDummyFeatures(M4_test_withKN,target = "DON")
  
  rf_DON_m4<-model_build(M4_train_withKN,"DON")
  
  ## map3 predict accuracy
  map4_predict<-predict(rf_DON_m4,newdata=M4_test_withKN)
  map4_train<-predict(rf_DON_m4,newdata=M4_train_withKN)
#  
  M4_rmse<-postResample(map4_predict$data$response,map4_predict$data$truth)[1]
  M4_r2<-postResample(map4_predict$data$response,map4_predict$data$truth)[2]
  M4_rmse_train<-postResample(map4_train$data$response,map4_train$data$truth)[1]
  M4_r2_train<-postResample(map4_train$data$response,map4_train$data$truth)[2]
  sing_acc<-data.frame(M1_r2,M2_r2,M4_r2,M1_r2_train,M2_r2_train,M4_r2_train)
  
  all_results<-rbind(all_results,sing_acc)
  
  print(all_results)
 
 }
 
 
