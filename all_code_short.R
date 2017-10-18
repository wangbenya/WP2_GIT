
## load the data 
all_results<-data.frame()
set.seed(66)
seed.list<-sample(1:1000,50,replace =F)

all_points<-read.csv("~/WP2/data/all_data1210.csv",header = T)
extra_n<-read.csv("~/WP2/data/extra_n.csv",header = T)
extra_n<-subset(extra_n,!(extra_n$WIN_Site_ID %in% all_points$WIN_Site_ID))

for (tt in c(1,2,3,13)){
  
  print(tt)
  seeds<-seed.list[tt]
  set.seed(seeds)
  trainIndex <- createDataPartition(all_points$DON, p = .8, list = FALSE)
  
  training <- all_points[trainIndex,]
  testing <- all_points[-trainIndex,]
  
  ## load the point data 
  training_df <- training[,c(1,2,3,5)] %>% rbind(.,extra_n[,c(1,2,3,5)]) %>%
    subset(.,.[,"DON"]!="NA") %>% read_pointDataframes(.) 
  
  testing_df <- testing[,c(1,2,3,5)] %>% read_pointDataframes(.) 
  
  training_points<-training[,c(1,2,3,5)] %>% rbind(.,extra_n[,c(1,2,3,5)]) %>%
    subset(.,.[,"DON"]!="NA") %>% read_points(.) 
  testing_points <- read_points(testing)
  
  ## map1, using kringing for DON interpolation
  f.1 <- as.formula(log10(DON) ~ s1+s2)
  
  # Add X and Y to training 
  training_df<-add_S1S2(training_df)
  testing_df<-add_S1S2(testing_df)
  
  # Compute the sample variogram; note that the f.1 trend model is one of the
  var.smpl1 <- variogram(f.1, training_df)
  # Compute the variogram model by passing the nugget, sill and range value
  dat.fit1 <- fit.variogram(var.smpl1, vgm(c("Exp","Sph","Gau","Lin","Spl")))
  # Perform the krige interpolation (note the use of the variogram model
  map1 <- krige(f.1, training_df, base_grid, dat.fit1) %>% raster(.) %>% raster::mask(., study_area)
  values(map1) <- 10 ^ (values(map1))
  dat.krg_DON<-map1
  
  
  ## M2, using RF to predict the DON
  landscape_train <- raster::extract(landscapes, training_points)
  landscape_test <- raster::extract(landscapes, testing_points)
  
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
  
  
  #set.seed(seeds)
  #rf_DON_m2 <- model_build(WP2Train, "DON","reg")
  
  #map2_predict <- predict(rf_DON_m2, newdata = WP2Test)
  #print(postResample(map2_predict$data$response, map2_predict$data$truth))
  
  ## map4, kriging first and then rf
  # kriging for DOC
  f.DOC <- as.formula(log10(DOC) ~ s1+s2)
  
  training_DOC <- training[,c(1,2,3,4)] %>% rbind(.,extra_n[,c(1,2,3,4)]) %>%
    subset(.,.[,"DOC"]!="NA") %>% read_pointDataframes(.) 
  training_DOC<-add_S1S2(training_DOC)
  var.smpl_DOC <- variogram(f.DOC, training_DOC)
  #plot(var.smpl_DOC)
  
  dat.fit_DOC <- fit.variogram(var.smpl_DOC,vgm(c("Exp","Sph","Gau","Lin","Spl")))
  #plot(var.smpl_DOC,dat.fit_DOC)
  # Perform the krige interpolation (note the use of the variogram model
  dat.krg_DOC <- krige(f.DOC, training_DOC, base_grid, dat.fit_DOC) %>% raster(.) %>% raster::mask(., study_area)
  values(dat.krg_DOC) <- 10 ^ (values(dat.krg_DOC))
  
  # kriging for NH4
  f.NH4 <- as.formula(log10(NH4) ~ s1+s2)
  
  training_NH4 <- training[,c(1,2,3,8)] %>% rbind(.,extra_n[,c(1,2,3,7)]) %>%
    subset(.,.[,"NH4"]!="NA") %>% read_pointDataframes(.) 
  training_NH4<-add_S1S2(training_NH4)
  
  var.smpl_NH4 <- variogram(f.NH4, training_NH4)
  
  #plot(var.smpl_NH4)
  dat.fit_NH4 <- fit.variogram(var.smpl_NH4, vgm(c("Exp","Sph","Gau","Lin","Spl")))
  
  # Perform the krige interpolation (note the use of the variogram model
  dat.krg_NH4 <- krige(f.NH4, training_NH4, base_grid, dat.fit_NH4) %>% raster(.) %>% raster::mask(., study_area)
  values(dat.krg_NH4) <- 10 ^ (values(dat.krg_NH4))
  
  # kriging for NOx
  f.NOx <- as.formula(log10(NOx) ~ s1+s2)
  
  training_NOx <- training[,c(1,2,3,7)] %>% rbind(.,extra_n[,c(1,2,3,6)]) %>%
    subset(.,.[,"NOx"]!="NA") %>% read_pointDataframes(.) 
  training_NOx<-add_S1S2(training_NOx)
  
  var.smpl_NOx <- variogram(f.NOx, training_NOx)
  #plot(var.smpl_NOx)
  
  dat.fit_NOx <- fit.variogram(var.smpl_NOx,vgm(c("Sph","Exp","Gau","Lin","Spl")))
  # Perform the krige interpolation 
  dat.krg_NOx <- krige(f.NOx, training_NOx, base_grid, dat.fit_NOx) %>% raster(.) %>% raster::mask(., study_area)
  values(dat.krg_NOx) <- 10 ^ (values(dat.krg_NOx))
  
  ## create rasterstack with kriging data
  kriging_nutrietn<-stack(dat.krg_DON,dat.krg_DOC, dat.krg_NH4, dat.krg_NOx)
  names(kriging_nutrietn) <- c("DON_k", "DOC_k", "NH4_k", "NOx_k")
  
  ## extract the data from landscapes_withN
  landscape_train_withKN <- raster::extract(kriging_nutrietn, read_points(base6[,15:17]))
  landscape_test_withKN <- raster::extract(kriging_nutrietn, read_points(test6[,15:17]))
  
  M4_train_withKN <- cbind(base6[, c(12,10,8, 6, 4,2, 13:17)],as.data.frame(landscape_train_withKN))
  M4_test_withKN <- cbind(test6[, c(12,10,8, 6, 4,2, 13:17)],as.data.frame(landscape_test_withKN))  
  names(M4_test_withKN) <- names(M4_train_withKN)
  
  M4_train_withKN<-M4_train_withKN[,-c(4)]
  M4_test_withKN<-M4_test_withKN[,-c(4)]

  #in_withKN <- reclass(M4_train_withKN,0.6,1.2)
  #M4_test_withKN <- reclass(M4_test_withKN,0.6,1.2)
  
  #M4_train_withKN <- reclass3(M4_train_withKN,0.6,1.2)
  #M4_test_withKN <- reclass3(M4_test_withKN,0.6,1.2)
 
  ## 

  set.seed(seeds)
  DON_rf_m4<-model_build(M4_train_withKN,"DON","reg")
  
  ## map3 predict accuracy
  map4_predict<-predict(DON_rf_m4,newdata=M4_test_withKN)
  
  #map4_predict$data$response=map4_predict$data$response*sd_train_DON+mean_train_DON
  #map4_predict$data$truth=map4_predict$data$truth*sd_train_DON+mean_train_DON
  
  print(postResample(map4_predict$data$response,map4_predict$data$truth))
  
  predict_results<-data.frame(seeds,map4_predict$data)
  
  all_results<-rbind(all_results,predict_results)
  
}


dim(all_results)

seeds<-unique(all_results$seeds)
length(seeds)


all_acc<-data.frame()

for (qq in seeds){
  print(qq)
  sub_data<-subset(all_results,all_results$seeds==qq)   
  #p1<-reclass4(sub_data[,c(3,2)],0.5,1.0)
  #p2<-reclass4(sub_data[,c(5,4)],0.5,1.0)
  p4<-reclass4(sub_data[,c(3,2)],0.5,1.0)
  print(table(p4[,2]))
  # print(table(p2[,2]))
  #  print(table(p4[,2]))
  sing_acc<-data.frame(p4=postResample(p4[,1],p4[,2])[1])
  
  all_acc<-rbind(all_acc,sing_acc)
}

print(all_acc)
