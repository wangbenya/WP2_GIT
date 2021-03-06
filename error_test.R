library(reshape2)
library(xlsx)

all_data<-read.csv("/home/ubuntu/WP2/data/WQ-WA-SWRISHIST.csv",header=T)
dim(all_data)
colnames(all_data)
WP4<-all_data[,c(1,2,4,8,18,19,20,32,33)]
head(WP4)
levels((WP4$Determinand))[c(190,227)]
WP4<-subset(WP4,WP4$Determinand %in% c(levels((WP4$Determinand))[c(190,227:229,420:424,429:441,446,447,807,808)]))

head(WP4)
WP4<-WP4[,-3]

WP4_n<-dcast(WP4,WIN.Site.Id+AWRC.Ref+Collect.Date+Sample.Depth.M+Sample.Depths.M+Sample.Depth.To.M~Determinand)
dim(WP4_n)

for (i in seq(4,32)){
  sub_data<-data.frame(WP4_n[,i])
  print(i)
  print(names(WP4_n)[i])
  print(dim(subset(sub_data,sub_data[,1]!=0)))
  print("")
  print("")
}

WP4_n<-WP4_n[,c(1:3,5,7,9,10,11,26,27,30,32)]

dim(WP4_n)

WP4_n<-WP4_n[,-4]

head(WP4_n)

## 

all_DOC<-subset(WP4_n,WP4_n$`C (sol org) {DOC, DOC as NPOC} | mg/L`!=0)
all_NOx<-subset(WP4_n,WP4_n$`N (sum sol ox) {NOx-N, TON} | mg/L`!=0)
all_NH4<-subset(WP4_n,WP4_n$`NH3-N/NH4-N (sol) | mg/L`!=0)
all_DON<-subset(WP4_n,WP4_n$`N (sum sol org) {DON} | mg/L`!=0)
all_TN<-subset(WP4_n,WP4_n$`N (tot) {TN, pTN} | mg/L`!=0)





write.csv(all_DOC,"~/WP2/data/all_DOC.csv",row.names = F)
write.csv(all_NOx,"~/WP2/data/all_NOx.csv",row.names = F)
write.csv(all_NH4,"~/WP2/data/all_NH4.csv",row.names = F)
write.csv(all_DON,"~/WP2/data/all_DON.csv",row.names = F)
write.csv(all_TN,"~/WP2/data/all_TN.csv",row.names = F)



all_site<-read.xls("/home/ubuntu/WP2/data/full_site_listing.xlsx",sheet = 2,header=T)



for (i1 in seq(0,1,0.1)){
  for (i2 in seq(0,1,0.1)){
    for (i3 in seq(0,1,0.1)){
     threshold<-c(Low=i1,Medium=i2,High=i3)
     pred2<-setThreshold(map2_predict,threshold)
     new_acc2<-postResample(pred2$data$response,pred2$data$truth)[1]
     
     pred4<-setThreshold(map4_predict,threshold)
     new_acc4<-postResample(pred4$data$response,pred4$data$truth)[1]
     if (new_acc4-M1_ACC>0){
       print(c(i1,i2,i3))
       print(c(M1_ACC,new_acc2,new_acc4))
     }
      
    }
  }
}



all_results<-read.csv("~/WP2/data/all_result_12.csv",header=T)

dim(all_results)

for (i in seq(1:6)){
  print(names(all_results)[i])
  print(mean(all_results[,i]))
}




plot(study_area)
points(training_df)
points(testing_df,col="red")

ggplot(data=as.data.frame(training_df),aes(x=training_df$DON))+
  geom_density()+geom_density(data=as.data.frame(testing_df),aes(x=testing_df$DON),col="red")+
  theme_bw()
summary(testing_df$DON)

cor(WP2Train[,-c(1:4)])
cor(WP2Test[,-c(1:4)])

ggplot(data=as.data.frame(training_df),aes(x=training_df$Longitude,y=training_df$Latitude))+
  geom_point(size=training_df$DON)+
  geom_point(data=as.data.frame(testing_df),aes(x=testing_df$Longitude,testing_df$Latitude),col="red",size=testing_df$DON)+
  theme_bw()

training_df$error<-(map2_train$data$truth-map2_train$data$response)^2
testing_df$error<-(map2_predict$data$truth-map2_predict$data$response)^2

as.data.frame(testing_df[,c(2,3,8,14)]) %>% round(.,3)

ggplot(data=all_points,aes(x="DON",y=all_points$DON))+geom_boxplot()


plot(map2_train$data$truth,map2_train$data$response)
plot(map2_predict$data$truth,map2_predict$data$response)


depth_k2<-raster::mask(depth_k,water,inverse=T)
plot(depth_k2)
points(402000,6495000)
points(402000,6475000)
points(402000,6460000)
points(402000,6448000)
points(402000,6403000)
distance_LP <-raster::mask(distance(left_up),study_area)
distance_LP@data@names<-"Distance_LP"


ggplot(data=all_points,aes(x=all_points$Longitude,y=all_points$Latitude))+
  geom_text(aes(label=all_points$Collect_Month,size=all_points$DON,col=all_points$Collect_Year))+
  theme_bw()

g_train<-M4_train_withKN
g_test<-M4_test_withKN
g_predict_test<-map4_predict
g_predict_train<-map4_train
g_training<-as.data.frame(training_df)
g_testing<-as.data.frame(testing_df)
g_training$error<-(g_predict_train$data$truth-g_predict_train$data$response)^2
g_testing$error<-(g_predict_test$data$truth-g_predict_test$data$response)^2



b_train<-M4_train_withKN
b_test<-M4_test_withKN
b_predict_test<-map4_predict
b_predict_train<-map4_train
b_training<-as.data.frame(training_df)
b_testing<-as.data.frame(testing_df)
b_training$error<-(b_predict_train$data$truth-b_predict_train$data$response)^2
b_testing$error<-(b_predict_test$data$truth-b_predict_test$data$response)^2


ggplot(data=g_testing,aes(x=g_testing$Longitude,y=g_testing$Latitude))+
  geom_point(col="red",aes(size=g_testing$error))+
geom_point(data=g_training,col="black",aes(x=g_training$Longitude,y=g_training$Latitude,size=g_training$error))

ggplot(data=b_testing,aes(x=b_testing$Longitude,y=b_testing$Latitude))+
  geom_point(col="red",aes(size=b_testing$error))+
  geom_point(data=b_training,col="black",aes(x=b_training$Longitude,y=b_training$Latitude,size=b_training$error))



ggplot(data=as.data.frame(g_predict_test$data),aes(x=g_predict_test$data$truth,y=g_predict_test$data$response))+
 geom_point()+lims(x=c(0,3),y=c(0,3))+geom_abline(slope = 1,intercept = 0,size=1)



ggplot(data=as.data.frame(b_predict_test$data),aes(x=b_predict_test$data$truth,y=b_predict_test$data$response))+
  geom_point()+lims(x=c(0,3),y=c(0,3))+geom_abline(slope = 1,intercept = 0,size=1)


head(M4_train_withKN)

tune_land<-read.csv("D://Program//work package 2//WP2_1128//WP2//data//tune_land_results.csv",header=T)

error_surface<-tune_land[,c(1:3)] 
train_surface<-tune_land[,c(1,2,4)] 

ggplot(data=error_surface,aes(x=error_surface$a,y=error_surface$b))+
  geom_tile(aes(fill=error_surface$M2_r2))+scale_fill_continuous(low = "grey",high = "black")+
  scale_x_continuous(breaks = seq(0,1500,100))+  scale_y_continuous(breaks = seq(0,1500,100))


ggplot(data=error_surface,aes(x=train_surface$a,y=train_surface$b))+
  geom_tile(aes(fill=train_surface$M2_r2))+scale_fill_continuous(low = "grey",high = "black")+
  scale_x_continuous(breaks = seq(0,1500,100))+  scale_y_continuous(breaks = seq(0,1500,100))



Res_df <- raster::mask(kriging_DON_res,study_area_withW) %>% 
  rasterToPoints(.) %>% data.frame(.)
#Make the points a dataframe for ggplot
#Make appropriate column headings
colnames(Res_df) <- c("Longitude", "Latitude", "residual")

#Now make the map
ggplot(data = Res_df, aes(y = Latitude, x = Longitude)) +
  geom_raster(aes(fill = residual)) + theme_bw() +
  coord_equal() +
  theme(panel.grid = element_blank(), legend.position = "right", legend.key = element_blank())+
  scale_fill_gradientn(colours = terrain.colors(10))+
  geom_point(data=data.frame(training_df),aes(x=training_df$s1,y=training_df$s2),size=2,shape=1)


all_results<-read.csv("D:/Program/work package 2/WP2_1128/WP2/data/all_results.csv",header=T)

all_results<-melt(all_results)

all_results_r2<-subset(all_results,all_results$variable %in% c('M1_r2','M2_r2','M4_r2'))
all_results_rmse<-subset(all_results,all_results$variable %in% c('M1_rmse','M2_rmse','M4_rmse'))

ggplot(data=all_results_r2,aes(all_results_r2$variable,all_results_r2$value))+geom_boxplot()+theme_bw()

ggplot(data=all_results_rmse,aes(all_results_rmse$variable,all_results_rmse$value))+geom_boxplot()+theme_bw()

