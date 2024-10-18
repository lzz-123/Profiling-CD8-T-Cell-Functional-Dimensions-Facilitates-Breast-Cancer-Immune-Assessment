#FuncAggre score development
#TCGA is the dataset containing 
#subtype: HR+HER2-, HER2+, TNBC
#stage: I, II, III
library(randomForestSRC)
#LOOCV
for (i in 1:nrow(TCGA)){
  train<-TCGA[-1,]
  test<-TCGA[-train,]
  fit<-rfsrc(Surv(PFI.time,PFI)~reactivity+cytotoxicity+IFNG+proliferation+apoptosis+as.factor(subtype)+as.factor(stage),data=train,ntree=500,seed=123,importance=TRUE)
  fit2<-predict(fit,test)
  error[i,1]<-fit2$err.rate[8]
}
#the model with lowest error rate was used and the relative importance of each functional dimension is derived for weighted aggregation
fit<-rfsrc(Surv(PFI.time,PFI)~reactivity+cytotoxicity+IFNG+proliferation+apoptosis+as.factor(subtype)+as.factor(stage),data=TCGA[-train,],ntree=500,seed=123,importance=TRUE)
fit[["importance"]]
#2 sets of weights are developed for non-metastatic and metastatic breast cancer, respectively

#FuncAggre score calculation
#prognosis
#for non-metastatic BC
FuncAggre<-reactivity*1+cytotoxicity*0.2+IFNG*0.2+proliferation*0.1-apoptosis
#for metastatic BC
FuncAggre<-reactivity*0.3+cytotoxicity*0.8+IFNG*0.7+proliferation*0.1-apoptosis*0.9

#immunotherapy response
#for non-metastatic BC
FuncAggre<-reactivity*1+cytotoxicity*0.2+IFNG*0.2+proliferation*0.1+apoptosis
#for metastatic BC
FuncAggre<-reactivity*0.3+cytotoxicity*0.8+IFNG*0.7+proliferation*0.1+apoptosis*0.9