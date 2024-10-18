#scRNA-seq based functional evaluation
#functional dimension evaluation at single cell level
#AUCell score calculation 
#matrix is the expression matrix of CD8+ T cells derived from Seurat.Object

#identification of CD8+ T cell exhibiting strong activity in tumor reactivity
#reactivity is the signature for tumor reactivity derived from 10.1126/science.abl5447)
library(AUCell)
set.seed(123)
AUC1<-AUCell_run(matrix,reactivity,BPPARAM=BiocParallel::MulticoreParam(5))
assignment1<-as.data.frame(AUC1@assays@data@listData)
assignment1<-as.data.frame(t(assignment1))
colnames(assignment1)[1]<-"AUC"

set.seed(123)
cells_assignment1 <- AUCell_exploreThresholds(AUC1,plotHist=TRUE, assign=TRUE)
assignment1$reactivity<-ifelse(assignment1$AUC>cells_assignment1[["V1"]][["aucThr"]][["thresholds"]][1],1,0)
assignment1<-aggregate(assignment1$reactivity,by=list(assignment1$orig.ident),FUN=sum,na.rm=TRUE)

#identification of CD8+ T cell exhibiting strong activity in cytotoxicity
#cytotoxicity is the signature for cytotoxicity derived from T cell mediated cytotoxicity (GO: 0001913)
set.seed(123)
AUC2<-AUCell_run(matrix,cytotoxicity,BPPARAM=BiocParallel::MulticoreParam(5))
assignment2<-as.data.frame(AUC2@assays@data@listData)
assignment2<-as.data.frame(t(assignment2))
colnames(assignment2)[1]<-"AUC"

set.seed(123)
cells_assignment2 <- AUCell_exploreThresholds(AUC2,plotHist=TRUE, assign=TRUE)
assignment2$cytotoxicity<-ifelse(assignment2$AUC>cells_assignment2[["V1"]][["aucThr"]][["thresholds"]][1],1,0)
assignment2<-aggregate(assignment2$cytotoxicity,by=list(assignment2$orig.ident),FUN=sum,na.rm=TRUE)

#identification of CD8+ T cell exhibiting strong activity in IFNG signaling
#IFNG is the signature for IFNG signaling derived from IFN-γ signaling (GO: 0060333)
set.seed(123)
AUC3<-AUCell_run(matrix,IFNG,BPPARAM=BiocParallel::MulticoreParam(5))
assignment3<-as.data.frame(AUC3@assays@data@listData)
assignment3<-as.data.frame(t(assignment3))
colnames(assignment3)[1]<-"AUC"

set.seed(123)
cells_assignment3 <- AUCell_exploreThresholds(AUC3,plotHist=TRUE, assign=TRUE)
assignment3$IFNG<-ifelse(assignment3$AUC>cells_assignment3[["V1"]][["aucThr"]][["thresholds"]][1],1,0)
assignment3<-aggregate(assignment3$IFNG,by=list(assignment3$orig.ident),FUN=sum,na.rm=TRUE)

#identification of CD8+ T cell exhibiting strong activity in proliferation
#proliferation is the signature for IFNG signaling derived from IFN-γ signaling (GO: 0060333)
set.seed(123)
AUC4<-AUCell_run(matrix,proliferation,BPPARAM=BiocParallel::MulticoreParam(5))
assignment4<-as.data.frame(AUC4@assays@data@listData)
assignment4<-as.data.frame(t(assignment4))
colnames(assignment4)[1]<-"AUC"

set.seed(123)
cells_assignment4 <- AUCell_exploreThresholds(AUC4,plotHist=TRUE, assign=TRUE)
assignment4$proliferation<-ifelse(assignment4$AUC>cells_assignment4[["V1"]][["aucThr"]][["thresholds"]][1],1,0)
assignment4<-aggregate(assignment4$proliferation,by=list(assignment4$orig.ident),FUN=sum,na.rm=TRUE)

#identification of CD8+ T cell exhibiting strong activity in apoptosis
#apoptosis is the signature for apoptosis derived from apoptotic process (GO: 0006915)
set.seed(123)
AUC5<-AUCell_run(matrix,apoptosis,BPPARAM=BiocParallel::MulticoreParam(5))
assignment5<-as.data.frame(AUC5@assays@data@listData)
assignment5<-as.data.frame(t(assignment5))
colnames(assignment5)[1]<-"AUC"

set.seed(123)
cells_assignment5 <- AUCell_exploreThresholds(AUC5,plotHist=TRUE, assign=TRUE)
assignment5$apoptosis<-ifelse(assignment5$AUC>cells_assignment5[["V1"]][["aucThr"]][["thresholds"]][1],1,0)
assignment5<-aggregate(assignment5$apoptosis,by=list(assignment5$orig.ident),FUN=sum,na.rm=TRUE)

#total is the total cell counts for each participant
total<-left_join(total,assignment1,by="orig.ident")
total<-left_join(total,assignment2,by="orig.ident")
total<-left_join(total,assignment3,by="orig.ident")
total<-left_join(total,assignment4,by="orig.ident")
total<-left_join(total,assignment5,by="orig.ident")

total$reactivity<-total$reactivity/total$total
total$cytotoxicity<-total$cytotoxicity/total$total
total$IFNG<-total$IFNG/total$total
total$proliferation<-total$proliferation/total$total
total$apoptosis<-total$apoptosis/total$total