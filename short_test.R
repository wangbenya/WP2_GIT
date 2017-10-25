
#dim(all_results)

#seeds<-unique(all_results$seeds)
#length(seeds)

#all_acc<-data.frame()

for (a1 in c(0.5,1.0,1.5)){
   for (a2 in c(1.0,1.5,2.0,2.5,3.0,3.5)){
  if (a2-a1>=0.5){
  p1<-reclass4(all_results[,c(3,2)],a1,a2)
  p2<-reclass4(all_results[,c(5,4)],a1,a2)
  p4<-reclass4(all_results[,c(7,6)],a1,a2)
  print(table(p1[,2]))
  sing_acc<-data.frame(p1=postResample(p1[,1],p1[,2])[1],p2=postResample(p2[,1],p2[,2])[1],p4=postResample(p4[,1],p4[,2])[1])
  print(a1)
  print(a2)
  print(sing_acc)
  }
 }}
  
