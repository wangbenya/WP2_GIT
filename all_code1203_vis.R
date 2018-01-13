  ## build the model for map2
  names(WP2Train)<-c("Soil", "Veg", "Landuse","SS","GS","Catchment", "GW_depth", "Distance", "DON","Longitude","Latitude")
  
  WP2Train<-reclass(WP2Train,a1,a2)
  
  WP2Train<-WP2Train[,-c(4,5,10,11)]
  
  ## scale the dataset using training data 
  name_list <- list("Soil", "Veg", "Landuse", "SS", "GS", "Catchment", "GW_depth", "Distance", "DON")
  value_list <- list(b1,b2,b3,b4,b5,b6)
  value_max_list <-list(soil_max,veg_max,landuse_max,ss_max,GS_max,cat_max)
  
  for (i in seq(1,6)) {
    sub_layer <- landscapes@layers[[i]]
    sub_data <- as.data.frame(sub_layer) 
    print(i)
    for (q in seq(1, 37044)) {
      print(q)
      if ((is.numeric(sub_data[q, ])) & (sub_data[q, ] %in% value_list[[i]][, 1])) {
        sub_data[q,] <- value_list[[i]][value_list[[i]][name_list[[i]]] == sub_data[q,], 2]
      } else if ((is.numeric(sub_data[q, ])) & (!(sub_data[q, ] %in% value_list[[i]][, 1]))) {
        sub_data[q,] <- value_max_list[[i]]
      }
    }
    values(landscapes@layers[[i]]) <- sub_data[, 1]
    names(landscapes@layers[[i]]) <- name_list[[i]]
  }
  
  plot(landscapes)

  WP2Train$Distance<-log10(WP2Train$Distance+0.01)

  WP2Train_raster<-stack(landscapes[[1]],landscapes[[2]],landscapes[[3]],landscapes[[6]],landscapes[[7]],landscapes[[8]])
  WP2Train_v<-as.data.frame(values(WP2Train_raster))
  WP2Train_v$Distance<-log10(WP2Train_v$Distance+0.01)
    
  for(i in c(1:6)){
    names(WP2Train_raster@layers[[i]]) <- names(WP2Train)[i]
    min_train<-min(WP2Train[,i])
    max_train<-max(WP2Train[,i])
    
    WP2Train[,i]<-(WP2Train[,i]-min_train)/(max_train-min_train)
    WP2Train_v[,i]<-(WP2Train_v[,i]-min_train)/(max_train-min_train)
    
    sd_train<-sd(WP2Train[,i])
    mean_train<-mean(WP2Train[,i])
    
    WP2Train[,i]<-(WP2Train[,i]-mean_train)/sd_train
    WP2Train_v[,i]<-(WP2Train_v[,i]-mean_train)/sd_train
    
  }
  
  set.seed(35)
  rf_DON_m2 <- model_build(WP2Train, "DON","cla")
  
  map2_predict <- predict(rf_DON_m2, newdata = WP2Train_v)
  
  map2<-dat.krg_DON
  values(map2)<-map2_predict$data$response
  #convert the raster to points for plotting
  map2_df<-raster::mask(map2,study_area_withW) %>% rasterToPoints(.) %>% data.frame(.)
  #Make appropriate column headings
  colnames(map2_df) <- c("Longitude", "Latitude", "DON")
  
  #Now make the map
  ggplot(data = map2_df, aes(y = Latitude, x = Longitude)) +
    geom_raster(aes(fill = as.factor(DON))) + theme_bw() +
    coord_equal() +
    theme(panel.grid = element_blank(), legend.position = "right", legend.key = element_blank())+
    scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9"),
                      name = "DON mg/L",
                      breaks = c("Low", "Medium", "High"),
                      labels = c("Low", "Medium", "High"))
  
  
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
  dat.krg_DOC <- krige(f.DOC, training_DOC, base_grid, dat.fit_DOC) %>% raster(.) 
  values(dat.krg_DOC) <- 10 ^ (values(dat.krg_DOC))
  
  ## create rasterstack with kriging data
  kriging_nutrietn<-stack(dat.krg_DOC)
  names(kriging_nutrietn) <- c("DOC_k")
  
  ## extract the data from landscapes_withN
  landscape_train_withKN <- raster::extract(kriging_nutrietn, read_points(base6[,15:17]))
  
  M4_train_withKN <- cbind(base6[, c(12,10,8, 6, 4,2, 13:17)],as.data.frame(landscape_train_withKN))
  
  ## create the training and testing sets 
  ## build the model for map2
  names(M4_train_withKN)[1:11]<-c("Soil", "Veg", "Landuse","SS","GS","Catchment", "GW_depth", "Distance", "DON","Longitude","Latitude")
  
  M4_train_withKN<-reclass(M4_train_withKN,a1,a2) %>% .[,-c(4,5,10,11)]

  M4_train_withKN$DOC_SOIL<-M4_train_withKN$DOC_k*M4_train_withKN$Soil
  M4_train_withKN$DOC_VEG<-M4_train_withKN$DOC_k*M4_train_withKN$Veg
  M4_train_withKN$DOC_LAND<-M4_train_withKN$DOC_k*M4_train_withKN$Landuse
  M4_train_withKN$DOC_CAT<-M4_train_withKN$Catchment*M4_train_withKN$DOC_k

  set.seed(35)
  rf_DON_m4<-model_build(M4_train_withKN,"DON","cla")
  
  ## map3 predict accuracy
  map4_predict<-predict(rf_DON_m4,newdata=M4_stack_v)
  #  map4_predict$data$response=map4_predict$data$response*sd_train_DON+mean_train_DON
  map4<-dat.krg_DON
  
  values(map4)<-map4_predict$data$response
  #convert the raster to points for plotting
  map4_df<-raster::mask(map4,study_area_withW) %>% rasterToPoints(.) %>% data.frame(.)
  #Make appropriate column headings
  colnames(map4_df) <- c("Longitude", "Latitude", "DON")
  
  #Now make the map
  ggplot(data = map4_df, aes(y = Latitude, x = Longitude)) +
    geom_raster(aes(fill = as.factor(DON))) + theme_bw() +
    coord_equal() +
    theme(panel.grid = element_blank(), legend.position = "right", legend.key = element_blank())+
    scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9"),
                      name = "DON mg/L",
                      breaks = c("Low", "Medium", "High"),
                      labels = c("Low", "Medium", "High")) 

  ### get high DON area 
  high_area<-map4_df
  high_area[high_area$DON==2,]<-1

  ggplot(data = map4_df, aes(y = map4_df$Latitude, x = map4_df$Longitude)) +
    geom_raster(aes(fill = as.factor(high_area$DON))) + theme_bw() +
    coord_equal() +
    theme(panel.grid = element_blank(), legend.position = "right", legend.key = element_blank())+
    scale_fill_manual(values = c("#999999","#56B4E9"),
                      name = "DON mg/L",
                      breaks = c("Low", "Medium", "High"),
                      labels = c("Low", "Medium", "High")) 
  
  ## get the uncertainty 
  m4_p<-map4_predict$data
  m4_p$main<-apply(m4_p[,1:3],1,max)
  m4_p$uncertainty<-1-m4_p$main
  
  map_uncertainty<-dat.krg_DON
  values(map_uncertainty)<-m4_p$uncertainty
  #convert the raster to points for plotting
  map_un_df<-raster::mask(map_uncertainty,study_area_withW) %>% rasterToPoints(.) %>% data.frame(.)
  
  ggplot(data = map4_df, aes(y = map4_df$Latitude, x = map4_df$Longitude)) +
    geom_raster(aes(fill = map_un_df$var1.pred)) + theme_bw() +
    coord_equal() +
    theme(panel.grid = element_blank(), legend.position = "right", legend.key = element_blank())+
    #geom_point(data=data.frame(training_df),aes(x=training_df$s1,y=training_df$s2),
     #          shape=5,size=1.5,col="red")+
    geom_point(data=data.frame(depth),aes(x=depth$s1,y=depth$s2),
               shape=3,size=2,col="red")
 
  # map4_predict$data$response=map4_predict$data$response*(max_train_DON-min_train_DON)+min_train_DON
  high_area2<-map4_df
  high_area2$uncertianty <-map_un_df$var1.pred
  high_area2$risk <-map_un_df$var1.pred
  high_area2$needSamp <-map_un_df$var1.pred

  high_area2[, "risk"][(high_area2[, "DON"] !=1)&(high_area2[,"uncertianty"]<=0.3)] <- "High_risk"
  high_area2[, "risk"][high_area2[, "risk"] !="High_risk"] <- "No_risk"
  
  high_area2[, "needSamp"][(high_area2[, "DON"] !=1)&(high_area2[,"uncertianty"]>=0.55)] <- "More_data"
  high_area2[, "needSamp"][high_area2[, "needSamp"] !="More_data"] <- "No_need"
  
  ## get high DON area with low uncertatiy 
  ggplot(data = map4_df, aes(y = map4_df$Latitude, x = map4_df$Longitude)) +
    geom_raster(aes(fill = as.factor(high_area2$risk))) + theme_bw() +
    coord_equal() +
    theme(panel.grid = element_blank(), legend.position = "right", legend.key = element_blank())+
    scale_fill_manual(values = c("#56B4E9","#D0D3D4"),
                      name = "DON mg/L",
                      breaks = c("Low", "Medium", "High"),
                      labels = c("Low", "Medium", "High")) +
  geom_point(data=data.frame(depth),aes(x=depth$s1,y=depth$s2),
             shape=3,size=3,col="red")
  
  ## get median and high DON areaa with high uncertainty 

  ggplot(data = map4_df, aes(y = map4_df$Latitude, x = map4_df$Longitude)) +
    geom_raster(aes(fill = as.factor(high_area2$needSamp))) + theme_bw() +
    coord_equal() +
    theme(panel.grid = element_blank(), legend.position = "right", legend.key = element_blank())+
    scale_fill_manual(values = c("#56B4E9","#D0D3D4"),
                      name = "DON mg/L",
                      breaks = c("Low", "Medium", "High"),
                      labels = c("Low", "Medium", "High"))+
    geom_point(data=data.frame(training_df),aes(x=training_df$s1,y=training_df$s2),
               shape=5,size=1)+
    geom_point(data=data.frame(depth),aes(x=depth$s1,y=depth$s2),
               shape=3,size=3,col="red")
  
  
  

    



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
parallelStartSocket(4)

WGS84 <- CRS("+proj=utm +zone=50 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

study_area <- shapefile("D://Program//work package 2//WP2_1128//WP2//data//study_area.shp")
water <- shapefile("D://Program//work package 2//WP2_1128//WP2//data//water.shp")

## load the veg, soil, land use, groundwater subarea,
## surface water subare, catchment
Soil <- raster("D://Program//work package 2//WP2_1128//WP2//data//soil1.ovr")
Veg <- raster("D://Program//work package 2//WP2_1128//WP2//data//vegetation.ovr")
Land_use <- raster("D://Program//work package 2//WP2_1128//WP2//data//landuse1.ovr")
Cat <- raster("D://Program//work package 2//WP2_1128//WP2//data//catch_name2.tif.ovr")
DEM<-raster("D://Program//work package 2//WP2_1128//WP2//data//topo_ProjectRaster3.tif.ovr")
## define the function 
## preprocess 
study_area <- spTransform(study_area, WGS84)
extent <- c(study_area@bbox[1, 1:2], study_area@bbox[2, 1:2])

water <- spTransform(water, WGS84)

# study_area_withW
study_area_withW <- symdif(study_area, water)

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



get_landscape2<-function(df){
  landscape_a<-data.frame()
  for (ii in seq(1,length(df))){
    aa<-as.data.frame(df[[ii]])
    aa<-subset(aa,aa$Soil!="NA")
    Soil=tail(names(sort(table(aa[,1]))),1)
    Veg=tail(names(sort(table(aa[,2]))),1)
    Landuse=tail(names(sort(table(aa[,3]))),1)
    GW_depth=mean(aa[,4])
    Distance_GWC=mean(aa[,5])
    slope=mean(aa[,6])
    s1=mean(aa[,7])
    s2=mean(aa[,8])
    sing_land<-data.frame(Soil,Veg,Landuse,GW_depth,Distance_GWC,slope,s1,s2)
    landscape_a<-rbind(landscape_a,sing_land)
  }
  return(landscape_a)
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

v_slope<-values(slope)
v_slope[v_slope %in% c(NA,NaN)]=0
values(slope)<-v_slope

v_aspect<-values(aspect)
v_aspect[v_aspect %in% c(NA,NaN)]=0
values(aspect)<-v_aspect


## M2, using RF to predict the DON
depth <- read.csv("D://Program//work package 2//WP2_1128//WP2//data//sampling_depth.csv",header=T) %>% read_pointDataframes(.)

# Define the 1st order polynomial equation
f_depth <- as.formula(sampling_d ~ 1)
# Add X and Y to training 
depth<-add_S1S2(depth)
# variogram on the de-trended data.
var.depth <- variogram(f_depth, depth)
plot(var.depth)
dat.fit_depth <- fit.variogram(var.depth,vgm(c("Sph","Exp")))
plot(var.depth,dat.fit_depth)

# created in the earlier step)
depth_k <- krige(f_depth, depth, base_grid, dat.fit_depth) %>% raster(.)
#plot(depth_k)
depth_k@data@names<-"GW_depth"

#Now make the map
### distance
water <- raster::rasterize(water, depth_k)
#water_distance <- raster::mask(distance(water),study_area)
water_distance <- raster::distance(water)

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

all_points<-read.csv("D://Program//work package 2//WP2_1128//WP2//data//all_data1127.csv",header = T)
extra_n<-read.csv("D://Program//work package 2//WP2_1128//WP2//data//extra_n.csv",header = T)
extra_n<-subset(extra_n,!(extra_n$WIN_Site_ID %in% all_points$WIN_Site_ID))
#DN_GW4<-read.csv("~/WP2_GIT/TN_GW4.csv",header = T)
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

all_points[all_points$DON==0.25,"DON"]=aaa


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
  makeDiscreteParam("ntree", values=seq(50,500,50)),
  makeIntegerParam("nodesize", lower = 60, upper = 65),
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

  training <- all_points 
  ## load the point data 
  training_df <- read_pointDataframes(training)
  
  training_points<- read_points(training)  
  ## map1, using kringing for DON interpolation
  f.1 <- as.formula(log10(DON) ~ 1)
  # Add X and Y to training 
  training_df<-add_S1S2(training_df)
  # Compute the sample variogram; note that the f.1 trend model is one of the
  var.smpl1 <- variogram(f.1, training_df)
  plot(var.smpl1)
  # Compute the variogram model by passing the nugget, sill and range value
 # dat.fit1 <- fit.variogram(var.smpl1,fit.method = FALSE,fit.sills = FALSE,fit.ranges = FALSE
  #                          ,vgm(model = "Sph",range = 5000,psill = 0.12,nugget = 0.05))
  dat.fit1 <- fit.variogram(var.smpl1,vgm(c("Sph","Exp")))
  
  plot(var.smpl1,dat.fit1)
  # Perform the krige interpolation (note the use of the variogram model
  kriging_DON_m1 <- krige(f.1, training_df, base_grid, dat.fit1) %>% raster(.) 
  values(kriging_DON_m1) <- 10 ^ (values(kriging_DON_m1))
  dat.krg_DON<-kriging_DON_m1

#convert the raster to points for plotting
  DON_map1_2<-raster::mask(kriging_DON_m1,study_area_withW)  
  map1_df <- rasterToPoints(DON_map1_2) %>% data.frame(.)
  #Make the points a dataframe for ggplot
  #Make appropriate column headings
  colnames(map1_df) <- c("Longitude", "Latitude", "DON")
  t=3
  map1_df[, t][map1_df[, t] <=a1] <- "Low"
  map1_df[, t][map1_df[, t] < a2] <- "Medium"
  map1_df[, t][(map1_df[, t] != "Low") & (map1_df[, t] != "Medium")] <- "High"
  map1_df[, t] <- factor(map1_df[, t], levels = c("Low", "Medium", "High"))
  
  
  #Now make the map
  p11<-ggplot(data = map1_df, aes(y = Latitude, x = Longitude)) +
    geom_raster(aes(fill = as.factor(DON))) + theme_bw() +
    coord_equal() +
    theme(panel.grid = element_blank(), legend.position = "right", legend.key = element_blank())+
    scale_fill_manual(values = c("#B8B8B8", "#FFD700", "#FF3030"),
                      name = "DON concentration",
                      breaks = c( "High","Medium", "Low"),
                      labels = c("High","Medium", "Low"))
  plot(p11)
   
  #convert the raster to points for plotting
  map1_df2 <- rasterToPoints(DON_map1_2) %>% data.frame(.)
  #Make the points a dataframe for ggplot
  #Make appropriate column headings
  colnames(map1_df2) <- c("Longitude", "Latitude", "DON")
  
  #Now make the map
  ggplot(data = map1_df2, aes(y = Latitude, x = Longitude)) +
    geom_raster(aes(fill = DON)) + theme_bw() +
    coord_equal() +
    theme(panel.grid = element_blank(), legend.position = "right", legend.key = element_blank())
    scale_fill_gradientn(colours = c("#B8B8B8", "#FFD700", "#FF3030"),breaks=c(0,1,2,3,4,5),limits=c(0,4))
  
  a=700
  b=1400
  
  capture_zone_land1<-function(df){
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
  
  capture_zone_land2<-function(df){
    num<-nrow(df)
    landscape_data<-data.frame()
    for (r in seq(1,num)){
      print(r)
      p1_long<-df@coords[r,1]
      p1_lat<-df@coords[r,2]
      pg<-spPolygons(rbind(c(p1_long,p1_lat),c(p1_long+a,p1_lat+b),c(p1_long+2*a,p1_lat+b),
                           c(p1_long+2*a,p1_lat-b),c(p1_long+a,p1_lat-b),c(p1_long,p1_lat)))  
      projection(pg)<- WGS84
      p1_landscape<-raster::extract(M2_map,pg)
      p1_landscape<-get_landscape2(p1_landscape)
      landscape_data<-rbind(landscape_data,p1_landscape)
    }
    return(landscape_data)
  }
  
  landscape_train <- capture_zone_land1(training_df)
  
  M2_train <- cbind(as.data.frame(landscape_train), training_df@data[c("DON","s1","s2")])
  
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
  
  b1=data.frame(table(M2_train$Soil))
  b2=data.frame(table(M2_train$Veg))
  b3=data.frame(table(M2_train$Landuse))
  b4=data.frame(table(M2_train$Catchment))
  value_list<-list(b1,b2,b3,b4)
  name_list <- list("Soil", "Veg", "Landuse","Catchment")

  landscapes_all<-landscapes
  for (ii in seq(1,4)) {
    sub_layer <- landscapes_all@layers[[ii]]
    sub_data <- as.data.frame(sub_layer) 
    print(ii)
    for (q in seq(1, 37044)) {
    if ((is.numeric(sub_data[q, ])) & (!(sub_data[q, ] %in% value_list[[ii]][, 1]))){
        sub_data[q,1] <- as.numeric(max_list[[ii]])
        print(q)
      }
    }
    values(landscapes_all@layers[[ii]]) <- sub_data[, 1]
    names(landscapes_all@layers[[ii]]) <- name_list[[ii]]
  }
  
  plot(landscapes_all)
  
  S1=landscapes@layers[[1]]
  values(S1)=rep(seq(357239,408058,length.out = 126),294)
  
  rep.row<-function(x,n){
    matrix(rep(x,each=n),nrow=n)
  }
  
  S2=landscapes@layers[[1]]
  values(S2)=as.vector(rep.row(seq(6394000,6530740,length.out = 294),126))
  
  M2_map<-stack(landscapes_all[[1]],landscapes_all[[2]],landscapes_all[[3]],landscapes_all[[5]],landscapes_all[[8]],landscapes_all[[9]],S1,S2)
  names(M2_map)[c(1,7,8)]<-c('Soil','s1','s2')
  
  #WP2Train_v<-as.data.frame(values(M2_map))
  
  #map2<-raster::mask(landscapes[[1]],study_area_withW)
  #convert the raster to points for plotting
  map2_all_points<- rasterToPoints(M2_map) %>% data.frame(.) 
  map2_all_points<-map2_all_points[,c(1,2,3)]
  #Make appropriate column headings
  colnames(map2_all_points) <- c("Longitude", "Latitude", "DON")
  SP <- SpatialPoints(map2_all_points[, 1:2], proj4string =  WGS84)
  SPD <- SpatialPointsDataFrame(SP, map2_all_points)
  SPD@bbox <- study_area@bbox
  
  landscape_all_points <- capture_zone_land2(SPD)
  
  landscape_all_points2<-cbind(landscape_all_points[,-c(7,8)],map2_all_points[,-3])
  #names(M2_map_data)[c(1,10,11)]<-c("Soil","s1","s2")
  M2_map_data<-landscape_all_points2
  names(M2_map_data)[c(7,8)]<-c("s1","s2")
  
  WP2Train<-M2_train[,-c(4,6,7,10)]
  
  WP2Train<-reclass(WP2Train,a1,a2)
  
  for(i in c("GW_depth","s1")) {
    
    min_train<-min(M2_train[,i])
    max_train<-max(M2_train[,i])
    
    M2_train[,i]<-(M2_train[,i]-min_train)/(max_train-min_train)
    M2_map_data[,i]<-(M2_map_data[,i]-min_train)/(max_train-min_train)
    
  }
  
  for(i in c("Distance_GWC","slope","s2")){
    
    min_train<-min(M2_train[,i])
    max_train<-max(M2_train[,i])
    
    M2_train[,i]<-(max_train-M2_train[,i])/(max_train-min_train)
    M2_map_data[,i]<-(max_train-M2_map_data[,i])/(max_train-min_train)
    
  }
  
  set.seed(666)
  rf_DON_m2 <- model_build2(WP2Train,"DON")
  
  map2_predict<-predict(rf_DON_m2,newdata=M2_map_data)
  
  map2<-landscapes_all[[1]]
  values(map2)<-map2_predict$data$response
  
  DON_map2_2<-raster::mask(map2,study_area_withW)
  
  map2_df <- rasterToPoints(DON_map2_2) %>% data.frame(.)
  #Make the points a dataframe for ggplot
  #Make appropriate column headings
  colnames(map2_df) <- c("Longitude", "Latitude", "DON")
  
  #Now make the map
  ggplot(data = map2_df, aes(y = Latitude, x = Longitude)) +
    geom_raster(aes(fill = DON)) + theme_bw() +
    coord_equal() +
    theme(panel.grid = element_blank(), legend.position = "right", legend.key = element_blank())+
    scale_fill_continuous(low = "lightgrey", high = "red2",limits=c(0,4))
  
  
  
  DON_map2_2<-raster::mask(map2,study_area_withW)
  
  map2_df <- rasterToPoints(DON_map2_2) %>% data.frame(.)
  #Make the points a dataframe for ggplot
  #Make appropriate column headings
  colnames(map2_df) <- c("Longitude", "Latitude", "DON")
  t=3
  map2_df[, t][map2_df[, t] <=a1] <- "Low"
  map2_df[, t][map2_df[, t] < a2] <- "Medium"
  map2_df[, t][(map2_df[, t] != "Low") & (map2_df[, t] != "Medium")] <- "High"
  map2_df[, t] <- factor(map2_df[, t], levels = c("Low", "Medium", "High"))
  
  #Now make the map
  ggplot(data = map2_df, aes(y = Latitude, x = Longitude)) +
    geom_raster(aes(fill = as.factor(DON))) + theme_bw() +
    coord_equal() +
    theme(panel.grid = element_blank(), legend.position = "right", legend.key = element_blank())+
    scale_fill_manual(values = c("#B8B8B8", "#FFD700", "#FF3030"),
                      name = "DON concentration",
                      breaks = c( "High","Medium", "Low"),
                      labels = c("High","Medium", "Low"))
  

  # kriging for DOC
  f.DOC <- as.formula(log10(DOC) ~ 1)
  
  training_DOC <- training[,c(1,2,3,7)] %>% rbind(.,extra_n[,c(1,2,3,4)]) %>%
    subset(.,.[,"DOC"]!="NA") %>% read_pointDataframes(.)
  
  training_DOC<-add_S1S2(training_DOC)
  var.smpl_DOC <- variogram(f.DOC, training_DOC)
  plot(var.smpl_DOC)
  
  dat.fit_DOC <- fit.variogram(var.smpl_DOC,vgm(c("Sph","Exp")))
  plot(var.smpl_DOC,dat.fit_DOC)
  # Perform the krige interpolation (note the use of the variogram model
  dat.krg_DOC <- krige(f.DOC, training_DOC, base_grid, dat.fit_DOC) %>% raster(.) 
  values(dat.krg_DOC) <- 10 ^ (values(dat.krg_DOC))
  
  ## create rasterstack with kriging data
  kriging_nutrietn_DOC<-stack(dat.krg_DOC,dat.krg_DON)
  names(kriging_nutrietn_DOC) <- c("DOC_k","DON_k")
  ## extract the data from landscapes_withN
  landscape_train_withKN <- raster::extract(kriging_nutrietn_DOC,training_df)
  landscape_all_points_DOC <-  raster::extract(kriging_nutrietn_DOC,SPD)

  M4_train_withKN <- cbind(WP2Train,as.data.frame(landscape_train_withKN))
  M4_train_all<-cbind(M2_map_data,as.data.frame(landscape_all_points_DOC))
  
  ## create the training and testing sets 
  #M4_test_withKN$DOC_dep<-M4_test_withKN$GW_depth*M4_test_withKN$DOC_k
  M4_train_withKN$DOC_k<-log10(M4_train_withKN$DOC_k)
  M4_train_all$DOC_k<-log10(M4_train_all$DOC_k)
  
  set.seed(666)
  rf_DON_m4<-model_build2(M4_train_withKN,"DON")
  
  #map4_predict$data$truth<-10^map4_predict$data$truth
  map4_predict<-predict(rf_DON_m4,newdata=M4_train_all)
  map4<-landscapes_all[[1]]
  values(map4)<-map4_predict$data$response
  
  DON_map4<-raster::mask(map4,study_area_withW)
  map4_df <- rasterToPoints(DON_map4) %>% data.frame(.)
  
  #Make the points a dataframe for ggplot
  #Make appropriate column headings
  colnames(map4_df) <- c("Longitude", "Latitude", "DON")
  
  #Now make the map
  p11<-ggplot(data = map4_df, aes(y = Latitude, x = Longitude)) +
    geom_raster(aes(fill = as.factor(DON))) + theme_bw() +
    coord_equal() +
    theme(panel.grid = element_blank(), legend.position = "right", legend.key = element_blank())+
    scale_fill_manual(values = c("#B8B8B8", "#FFD700", "#FF3030"),
                      name = "DON concentration",
                      breaks = c( "High","Medium", "Low"),
                      labels = c("High","Medium", "Low"))
  plot(p11)
  
  
  
  
  GW_depth_df <- raster::mask(landscapes[[5]],study_area_withW) %>% 
                  rasterToPoints(.) %>% data.frame(.)
  #Make the points a dataframe for ggplot
  #Make appropriate column headings
  colnames(GW_depth_df) <- c("Longitude", "Latitude", "GW_depth")
  
  #Now make the map
  ggplot(data = GW_depth_df, aes(y = Latitude, x = Longitude)) +
    geom_raster(aes(fill = GW_depth)) + theme_bw() +
    coord_equal() +
    theme(panel.grid = element_blank(), legend.position = "right", legend.key = element_blank())+
    scale_fill_continuous(low = "yellow", high = "red")
  
  
  Distance_df <- raster::mask(landscapes[[6]],study_area_withW) %>% 
    rasterToPoints(.) %>% data.frame(.)
  #Make the points a dataframe for ggplot
  #Make appropriate column headings
  colnames(Distance_df) <- c("Longitude", "Latitude", "distance")
  
  #Now make the map
  ggplot(data = Distance_df, aes(y = Latitude, x = Longitude)) +
    geom_raster(aes(fill = distance)) + theme_bw() +
    coord_equal() +
    theme(panel.grid = element_blank(), legend.position = "right", legend.key = element_blank())+
    scale_fill_gradientn(colours = terrain.colors(5),limits=c(0,60000))
  
  
  Distance_gwc_df <- raster::mask(landscapes[[7]],study_area_withW) %>% 
    rasterToPoints(.) %>% data.frame(.)
  #Make the points a dataframe for ggplot
  #Make appropriate column headings
  colnames(Distance_gwc_df) <- c("Longitude", "Latitude", "distance_gwc")
  
  #Now make the map
  ggplot(data = Distance_gwc_df, aes(y = Latitude, x = Longitude)) +
    geom_raster(aes(fill = distance_gwc)) + theme_bw() +
    coord_equal() +
    theme(panel.grid = element_blank(), legend.position = "right", legend.key = element_blank())+
    scale_fill_gradientn(colours = terrain.colors(5),limits=c(0,60000))
