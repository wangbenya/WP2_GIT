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
  geom_point(size=training_df$error)+
  geom_point(data=as.data.frame(testing_df),aes(x=testing_df$Longitude,testing_df$Latitude),col="red",size=testing_df$error)+
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

GW_center<-data.frame(Latitude=c(6495000,6475000,6460000,6448000,6403000),LongitudeLatitu=rep(402000,5),values=1)

GW_center <- SpatialPoints(GW_center[, c(2:1)], proj4string = WGS84)
GW_center@bbox <- study_area@bbox
  
distance_LP@data@names<-"Distance_LP"
  

d1<-distanceFromPoints(left_up,GW_center)


plot(d1)


left_up<-water 
values(left_up)<-1
values(left_up)[[1]]<-1

