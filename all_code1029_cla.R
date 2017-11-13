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
ss <- raster("~/WP2/data/Ssubarea.ovr")
gs <- raster("~/WP2/data/Gsubarea.ovr")
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


reclass2 <- function(df) {
  df[, "DON"][df[, "DON"] == 1] <- "Low"
  df[, "DON"][df[, "DON"] == 2] <- "Medium"
  df[, "DON"][(df[, "DON"] != "Low") & (df[, "DON"] != "Medium")] <- "High"
  df[, "DON"] <- factor(df[, "DON"], levels = c("Low", "Medium", "High"))
  return(df)
}

reclass3 <- function(df, i, j) {
  df[, "DON_k"][df[, "DON_k"] <= i] <- "Low"
  df[, "DON_k"][df[, "DON_k"] < j] <- "Medium"
  df[, "DON_k"][(df[, "DON_k"] != "Low") & (df[, "DON_k"] != "Medium")] <- "High"
  df[, "DON_k"] <- factor(df[, "DON_k"], levels = c("Low", "Medium", "High"))
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
    SS=tail(names(sort(table(aa[,4]))),1)
    GS=tail(names(sort(table(aa[,5]))),1)
    Catchment=tail(names(sort(table(aa[,6]))),1)
    GW_depth=mean(aa[,7])
    Distance=mean(aa[,8])
    
    sing_land<-data.frame(Soil,Veg,Landuse,SS,GS,Catchment,GW_depth,Distance)
    landscape_all<-rbind(landscape_all,sing_land)
  }
  return(landscape_all)
}


## preprocess the landscape raster
Soil <- pre(Soil)
Veg <- pre(Veg)
Land_use <- pre(Land_use)
Cat <- pre(Cat)
ss <- pre(ss)
gs <- pre(gs)

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
landscapes<-stack(Soil,Veg,Land_use,ss,gs,Cat,depth_k,water_distance)
names(landscapes) <- c("Soil", "Veg", "Landuse","SS","GS", "Catchment", "GW_depth", "Distance")


## load the data 
set.seed(666)
seed.list<-sample(1:1000,50,replace =F)

all_points<-read.csv("~/WP2/data/all_data1210.csv",header = T)
#all_g<-read.csv("~/WP2_GIT/all_g.csv",header=T)
#all_g<-rbind(all_points[,c(1,2,3,4,5,7,8)],extra_all_N[,c(1,2,3,5,4,6,7)],GW_extra[,c(1,2,3,4,5,6,7)])
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
      class_rf = makeLearner("classif.RRF")
      #class_rf$par.vals<-list(importance=T)
      ctrl = makeTuneControlIrace(maxExperiments = 200L)
      rdesc = makeResampleDesc("CV", iters = 5)
      
      ## define the parameter spaces for RF      
      para_rf = makeParamSet(
        makeDiscreteParam("ntree", values=seq(300,500,20)),
        makeIntegerParam("nodesize", lower = 4, upper = 10),
        makeIntegerParam("mtry", lower = 4, upper =8),
        makeDiscreteParam("coefReg", values=seq(0.05,0.5,0.05))
      )

      model_build <- function(dataset, n_target) {
        set.seed(35)
          ## define the regression task for DON 
          WP3_target = makeClassifTask(id = "WP3_target", data = dataset, target = n_target)
          ## cross validation
          ## 10-fold cross-validation
          rin = makeResampleInstance(rdesc, task = WP3_target)
          res_rf = mlr::tuneParams(class_rf, WP3_target, resampling = rdesc, par.set = para_rf, control = ctrl,
                                   show.info = FALSE, measures = kappa)
          lrn_rf = setHyperPars(class_rf, par.vals = res_rf$x)
        
        ## train the final model 
        set.seed(35)
        rf <- mlr::train(lrn_rf, WP3_target)
        return(rf)
      }
    
    a1=1
    a2=2
    print(a1)
    print(a2)
 all_results<-data.frame()

for (tt in c(1:20)){
  
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
  #plot(var.smpl1)
  # Compute the variogram model by passing the nugget, sill and range value
  dat.fit1 <- fit.variogram(var.smpl1,vgm(c("Sph","Exp")))
  
  #plot(var.smpl1,dat.fit1)
  # Perform the krige interpolation (note the use of the variogram model
  kriging_DON_m1 <- krige(f.1, training_df, base_grid, dat.fit1) %>% raster(.) %>% raster::mask(., study_area)
  values(kriging_DON_m1) <- 10 ^ (values(kriging_DON_m1))
  dat.krg_DON<-kriging_DON_m1
  
  map1_predict <- data.frame(observed_DON=testing_df@data$DON,predicted_DON=raster::extract(kriging_DON_m1, testing_points))
  
      for (t in c(1,2)){
    map1_predict[, t][map1_predict[, t] <=a1] <- "Low"
    map1_predict[, t][map1_predict[, t] < a2] <- "Medium"
    map1_predict[, t][(map1_predict[, t] != "Low") & (map1_predict[, t] != "Medium")] <- "High"
    map1_predict[, t] <- factor(map1_predict[, t], levels = c("Low", "Medium", "High"))
    
       }
  
  M1_ACC<-postResample(map1_predict[,2],map1_predict[,1])[1]
  M1_kappa<-postResample(map1_predict[,2],map1_predict[,1])[2]
  
  map1_train <- data.frame(observed_DON=training_df@data$DON,predicted_DON=raster::extract(kriging_DON_m1, training_points))
  
      for (t in c(1,2)){
    map1_train[, t][map1_train[, t] <=a1] <- "Low"
    map1_train[, t][map1_train[, t] < a2] <- "Medium"
    map1_train[, t][(map1_train[, t] != "Low") & (map1_train[, t] != "Medium")] <- "High"
    map1_train[, t] <- factor(map1_train[, t], levels = c("Low", "Medium", "High"))
    
       }

    M1_ACC_train<-postResample(map1_train[,2],map1_train[,1])[1]

  ## M2, using RF to predict the DON
  landscape_train <- raster::extract(landscapes, training_points,buffer=800)
  landscape_test <- raster::extract(landscapes, testing_points,buffer=800)
  
  landscape_train<-get_landscape(landscape_train)
  landscape_test<-get_landscape(landscape_test)
  
  M2_train <- cbind(as.data.frame(landscape_train), training_df@data[c("DON","Longitude","Latitude")])
  M2_test <- cbind(as.data.frame(landscape_test), testing_df@data[c("DON","Longitude","Latitude")])
  
  names(M2_train) <- colnames(M2_test)
  
  transCat<-function(i){
  var_name<-names(M2_train)[i]
  a<-as.list(unique(M2_train[,var_name]))
  b<- data.frame()
  
  for (x in seq(1, length(a))) {
    value = a[[x]]
    data1 <- subset(M2_train, M2_train[,var_name] == value)
    c <- data.frame(value, mean(data1$DON), length(data1$DON))
    b <- rbind(b,c)
  }
  
  colnames(b) <- c(var_name, paste0("mean_DON_",var_name),"NO.")
  return(b)
 }


  b1<-transCat(1)
  b2<-transCat(2)
  b3<-transCat(3)
  b4<-transCat(4)
  b5<-transCat(5)
  b6<-transCat(6)
  
  soil_max = b1[which(b1$NO. == max(b1$NO.)),][, 2]
  veg_max = b2[which(b2$NO. == max(b2$NO.)),][, 2]
  landuse_max = b3[which(b3$NO. == max(b3$NO.)),][, 2]
  ss_max = b4[which(b4$NO. == max(b4$NO.)),][, 2]
  GS_max = b5[which(b5$NO. == max(b5$NO.)),][, 2]
  cat_max = b6[which(b6$NO. == max(b6$NO.)),][, 2]
  
  base1 <- (merge(b1, M2_train, by = 'Soil', all.y = T))
  base2 <- (merge(b2, base1, by = 'Veg', all.y = T))
  base2 <- base2[, - c(3, 6)]
  
  test1 <- (merge(b1, M2_test, by = 'Soil', all.y = T))
  test1[is.na(test1)] <- soil_max
  test1 <- test1[, - c(3)]
  
  test2 <- (merge(b2, test1, by = 'Veg', all.y = T))
  test2[is.na(test2)] <- veg_max
  test2 <- test2[, - c(3)]
  
  base3 <- (merge(b3, base2, by = 'Landuse', all.y = T))
  base4 <- (merge(b4, base3, by = 'SS', all.y = T))
  base4 <- base4[, - c(3, 6)]
  
  test3 <- (merge(b3, test2, by = 'Landuse', all.y = T))
  test3[is.na(test3)] <- landuse_max
  test3 <- test3[, - c(3)]
  
  test4 <- (merge(b4, test3, by = 'SS', all.y = T))
  test4[is.na(test4)] <- ss_max
  test4 <- test4[, - c(3)]
  
  base5 <- (merge(b5, base4, by = 'GS', all.y = T))
  base6 <- (merge(b6, base5, by = 'Catchment', all.y = T))
  base6 <- base6[, - c(3, 6)]
  
  test5 <- (merge(b5, test4, by = 'GS', all.y = T))
  test5[is.na(test5)] <- GS_max
  test5 <- test5[, - c(3)]
  
  test6 <- (merge(b6, test5, by = 'Catchment', all.y = T))
  test6[is.na(test6)] <- cat_max
  test6 <- test6[, - c(3)]
  
  ## create the training and testing sets 
  WP2Train <- base6[, c(12,10,8, 6, 4,2, 13:17)]
  WP2Test <- test6[, c(12,10,8, 6, 4,2, 13:17)]
  
  ## build the model for map2
  names(WP2Train)<-c("Soil", "Veg", "Landuse","SS","GS","Catchment", "GW_depth", "Distance", "DON","Longitude","Latitude")
  names(WP2Test)<-c("Soil",  "Veg", "Landuse","SS","GS", "Catchment", "GW_depth", "Distance", "DON","Longitude","Latitude")
  
  WP2Train<-reclass(WP2Train,a1,a2)
  WP2Test<-reclass(WP2Test,a1,a2)
  
  WP2Train<-WP2Train[,-c(4,5,10,11)]
  WP2Test<-WP2Test[,-c(4,5,10,11)]
  
  WP2Train$Distance<-log10(WP2Train$Distance+0.01)
  
  WP2Test$Distance<-log10(WP2Test$Distance+0.01)
  
  
  # WP2Train$dep_soil<-WP2Train$Soil*WP2Train$GW_depth
  # WP2Train$dep_veg<-WP2Train$Veg*WP2Train$GW_depth
  # WP2Train$dep_land<-WP2Train$Landuse*WP2Train$GW_depth
  # WP2Train$dep_cat<-WP2Train$Catchment*WP2Train$GW_depth
  
   #WP2Test$dep_soil<-WP2Test$Soil*WP2Test$GW_depth
   #WP2Test$dep_veg<-WP2Test$Veg*WP2Test$GW_depth
   #WP2Test$dep_land<-WP2Test$Landuse*WP2Test$GW_depth
   #WP2Test$dep_cat<-WP2Test$Catchment*WP2Test$GW_depth
  
  #WP2Train$log_lat<-WP2Train$Longitude/WP2Train$Latitude
  #WP2Test$log_lat<-WP2Test$Longitude/WP2Test$Latitude
  
      for(i in c(1:6)){

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
  rf_DON_m2 <- model_build(WP2Train,"DON")
  
  map2_predict <- predict(rf_DON_m2, newdata = WP2Test)
  map2_train <- predict(rf_DON_m2, newdata = WP2Train)

  M2_ACC<-postResample(map2_predict$data$response, map2_predict$data$truth)[1]
  M2_kappa<-postResample(map2_predict$data$response, map2_predict$data$truth)[2]
  
  M2_ACC_train<-postResample(map2_train$data$response, map2_train$data$truth)[1]

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
  landscape_train_withKN <- raster::extract(kriging_nutrietn, read_points(base6[,15:17]))
  landscape_test_withKN <- raster::extract(kriging_nutrietn, read_points(test6[,15:17]))
  
  M4_train_withKN <- cbind(base6[, c(12,10,8, 6, 4,2, 13:17)],as.data.frame(landscape_train_withKN))
  M4_test_withKN <- cbind(test6[, c(12,10,8, 6, 4,2, 13:17)],as.data.frame(landscape_test_withKN))  
  names(M4_test_withKN) <- names(M4_train_withKN)
  
  
  ## create the training and testing sets 
  ## build the model for map2
  names(M4_train_withKN)[1:11]<-c("Soil", "Veg", "Landuse","SS","GS","Catchment", "GW_depth", "Distance", "DON","Longitude","Latitude")
  names(M4_test_withKN)[1:11]<-c("Soil",  "Veg", "Landuse","SS","GS", "Catchment", "GW_depth", "Distance", "DON","Longitude","Latitude")
  
  M4_train_withKN<-reclass(M4_train_withKN,a1,a2)
  M4_test_withKN<-reclass(M4_test_withKN,a1,a2)
  
  #M4_train_withKN <- reclass3(M4_train_withKN,0.5,1.0)
  #M4_test_withKN <- reclass3(M4_test_withKN,0.5,1.0)
  
  M4_train_withKN<-M4_train_withKN[,-c(4,5,10,11,12)]
  M4_test_withKN<-M4_test_withKN[,-c(4,5,10,11,12)]
  
  M4_train_withKN$DOC_SOIL<-M4_train_withKN$DOC_k*M4_train_withKN$Soil
  M4_train_withKN$DOC_VEG<-M4_train_withKN$DOC_k*M4_train_withKN$Veg
  M4_train_withKN$DOC_LAND<-M4_train_withKN$DOC_k*M4_train_withKN$Landuse
  M4_train_withKN$DOC_CAT<-M4_train_withKN$Catchment*M4_train_withKN$DOC_k
  
  M4_test_withKN$DOC_SOIL<-M4_test_withKN$DOC_k*M4_test_withKN$Soil
  M4_test_withKN$DOC_VEG<-M4_test_withKN$DOC_k*M4_test_withKN$Veg
  M4_test_withKN$DOC_LAND<-M4_test_withKN$DOC_k*M4_test_withKN$Landuse
  M4_test_withKN$DOC_CAT<-M4_test_withKN$Catchment*M4_test_withKN$DOC_k
  
  #M4_train_withKN$dep_soil<-M4_train_withKN$Soil*M4_train_withKN$GW_depth
  #M4_train_withKN$dep_veg<-M4_train_withKN$Veg*M4_train_withKN$GW_depth
  #M4_train_withKN$dep_land<-M4_train_withKN$Landuse*M4_train_withKN$GW_depth
  #M4_train_withKN$dep_cat<-M4_train_withKN$Catchment*M4_train_withKN$GW_depth
  
  #M4_test_withKN$dep_soil<-M4_test_withKN$Soil*M4_test_withKN$GW_depth
  #M4_test_withKN$dep_veg<-M4_test_withKN$Veg*M4_test_withKN$GW_depth
  #M4_test_withKN$dep_land<-M4_test_withKN$Landuse*M4_test_withKN$GW_depth
  #M4_test_withKN$dep_cat<-M4_test_withKN$Catchment*M4_test_withKN$GW_depth
  #M4_train_withKN$log_lat<-M4_train_withKN$Longitude/M4_train_withKN$Latitude
  #M4_test_withKN$log_lat<-M4_test_withKN$Longitude/M4_test_withKN$Latitude
  
 # M4_train_withKN$DON<-log10(M4_train_withKN$DON)
  M4_train_withKN$Distance<-log10(M4_train_withKN$Distance+0.01)
  
  M4_test_withKN$Distance<-log10(M4_test_withKN$Distance+0.01)
   
  M4_train_withKN$DOC_dep<-M4_train_withKN$GW_depth*M4_train_withKN$DOC_k
  M4_test_withKN$DOC_dep<-M4_test_withKN$GW_depth*M4_test_withKN$DOC_k


      for(i in c(1:6,8:13)){
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
  rf_DON_m4<-model_build(M4_train_withKN,"DON")
  
  ## map3 predict accuracy
  map4_predict<-predict(rf_DON_m4,newdata=M4_test_withKN)
  map4_train<-predict(rf_DON_m4,newdata=M4_train_withKN)

  #  
  M4_kappa<-postResample(map4_predict$data$response,map4_predict$data$truth)[2]
  M4_ACC<-postResample(map4_predict$data$response,map4_predict$data$truth)[1]
  M4_ACC_train<-postResample(map4_train$data$response,map4_train$data$truth)[1]

  sing_acc<-data.frame(M1_ACC,M2_ACC,M4_ACC,M1_ACC_train,M2_ACC_train,M4_ACC_train)
  
  all_results<-rbind(all_results,sing_acc)
  
  print(all_results)
  print(rf_DON_m2$learner)
  print(rf_DON_m4$learner)
 
    }
  
