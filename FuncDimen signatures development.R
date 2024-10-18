#FuncDimen signatures development
#For each functional dimension, the same gene selection strategy is applied
#step 1 spearman correlation
#matrix is the combination of expression profile and abundance 
for (i in 1:m){
  a<-cor.test(matrix[,i+6],matrix[,2],method="spearman")
  cor[i,1]<-a$estimate
  p[i,1]<-a$p.value
  a<-cor.test(matrix[,i+6],matrix[,3],method="spearman")
  cor[i,2]<-a$estimate
  p[i,2]<-a$p.value
  a<-cor.test(matrix[,i+6],matrix[,4],method="spearman")
  cor[i,3]<-a$estimate
  p[i,3]<-a$p.value
  a<-cor.test(matrix[,i+6],matrix[,5],method="spearman")
  cor[i,4]<-a$estimate
  p[i,4]<-a$p.value
  a<-cor.test(matrix[,i+6],matrix[,6],method="spearman")
  cor[i,5]<-a$estimate
  p[i,5]<-a$p.value
}
#FDR
for (i in 1:5){
  p[,i]<-p.adjust(p[,i],method="fdr")
}

#step 2 FuncDImen signature generation
#matrixA is the expression matrix of the genes correlated with tumor reactivity
#Y1 is the abundance of CD8+ T cells exhibiting strong activity in tumor reactivity
X1<-matrixA
Y1<-total$reactivity
#lasso
set.seed(123)
lamodel <- glmnet(X1, Y1,family="gaussian") 
#cross validation for lasso
set.seed(123)
lasso_model <- cv.glmnet(X1,Y1,nfolds =10)
best_lambda <- lasso_model$lambda.min
set.seed(123)
best_model1 <- glmnet(X1, Y1, alpha = 1, lambda = best_lambda1, family="gaussian")
selection1a1<-as.matrix(coef(best_model1))
selection1a1<-as.data.frame(selection1a1)
selection1a1$Gene<-rownames(selection1a1)
selection1a1<-selection1a1[-which(selection1a1$s0==0),]

#randomforest
#matrixA is the gene expression matrix and is integrated with the abundance of CD8+ T cells exhibiting strong activity in tumor reactivity
library(randomForest)
set.seed(123)
rf<- randomForest(x=X1,y=Y1,ntree=500,important=TRUE,proximity=TRUE)
set.seed(123)
cv1<-rfcv(matrixA[,-1], matrixA$reactivity, cv.fold = 10)
with(cv1, plot(n.var, error.cv, log="x", type="o", lwd=2))
set.seed(123)
cv <- replicate(10, rfcv(matrixA[,-1], matrixA$reactivity, cv.fold = 10), simplify = FALSE)
cvtest <- data.frame(sapply(cv, '[[', 'error.cv'))
cvtest$otus <- rownames(cvtest)
cvtest <- reshape2::melt(cvtest, id = 'otus')
cvtest$otus <- as.numeric(as.character(cvtest$otus))
cvmean <- aggregate(cvtest$value, by = list(cvtest$otus), FUN = mean)
selection1b1<-as.data.frame(rf$importance)
selection1b1$Gene<-rownames(selection1b1)
selection1b1<-arrange(selection1b1,IncNodePurity)
#the optimal number of critical genes selected by randomforest
selection1b1<-as.data.frame(selection1b1[c(173:309),])

#xgboost
library("xgboost")
rownames(matrixA)<-matrixA$patient
#remove the column of patients
matrixA<-matrixA[,-2]
model_martix_train1 <- model.matrix(reactivity ~ . - 1, matrixA)
data_train1 <- xgb.DMatrix(model_martix_train1, label = matrixA$reactivity,nthread=2)
param <-list(max_depth = 6, eta = 0.3,objective = "reg:linear")
set.seed(123)
cv<-xgb.cv(param,data_train1,nrounds=500,nfold=10)
set.seed(123)
xgb_model1<- xgb.train(param, data_train1, nrounds = x)#the iteration based on the cross-validation results
xgb1 <- xgb.importance(xgb_model1[["feature_names"]], model = xgb_model1)
selection1c1<-xgb.importance(xgb_model1[["feature_names"]], model = xgb_model1)

#selected by any of 2 machine learning algorithms
inter1<-selection1a1[selection1a1$Gene%in%selection1b1$Gene,]
inter2<-selection1a1[selection1a1$Gene%in%selection1c1$Gene,]
inter3<-selection1b1[selection1b1$Gene%in%selection1c1$Gene,]
inter<-rbind(inter1,inter2)
inter<-rbind(inter,inter3)
inter<-distinct(inter,Gene,.keep_all=TRUE)

#step3 FuncDimen signatures refinement
#size is the gene number 
library(caret)
data_train1a <- matrixA[,colnames(matrixA) %in% inter$Gene]
data_train1a$reactivity<-matrixA$reactivity
set.seed(123)
sig1<-rfe(x = data_train1a[,-1],
          y = data_train1a$reactivity,
          sizes = c(1:24),
          rfeControl = rfeControl(functions = rfFuncs,
                                  number=10,method = 'LOOCV'))
sig1<-as.data.frame(sig1$optVariables)

#FuncDimen score calculation
#tpm is the gene expression matrix. column is ID and row is gene name
#formula: average expression of the genes that are positively correlated with functional dimension minus average expression of the genes that are negatively correlated with functional dimension)
#score1a is average expression of the genes that are positively correlated with functional dimension
#score1b is average expression of the genes that are negatively correlated with functional dimension

#tumor reactivity
up1<-c("MARCO","R3HDM1","UBD","CXCR6","KHDRBS1","HLA-DQB1","B4GALT5","CD8B","OPTN","C12orf75","RNF19A","RUNX3","ZFP69B","CCDC65","GBP1")
down1<-c("COA3","TMEM54","SCNN1A","ADIRF","RHOD","USP27X-AS1")
score1a<-tpm[tpm$geneSymbol%in%$up1,]
score1a<-as.data.frame(t(score1a))
score1a<-apply(score1a,1,mean)
score1b<-tpm[tpm$geneSymbol%in%$down1,]
score1b<-as.data.frame(t(score1b))
score1b<-apply(score1b,1,mean)
score1<-left_join(score1a,score1b,by="ID")
score1$reactivity<-score1$score1a-score1$score1b

#cytotoxicity
up2<-c("HLA-DQB1","GRAP2","GBP2","FASLG","RNF19A","ZFF683","FNBP1","CST1","C12orf75","GZMA","GRSF1","IKZF4","CCDC59")
down2<-c("NT5C3B","MPC2","PRSS22","JPH1","USP27X-AS1")
score2a<-tpm[tpm$geneSymbol%in%$up2,]
score2a<-as.data.frame(t(score2a))
score2a<-apply(score2a,1,mean)
score2b<-tpm[tpm$geneSymbol%in%$down2,]
score2b<-as.data.frame(t(score2b))
score2b<-apply(score2b,1,mean)
score2<-left_join(score2a,score2b,by="ID")
score2$cytotoxicity<-score2$score2a-score2$score2b

#IFNG
up3<-c("HLA-DQA1","CD8A","B4GALT5","GPA33","MIIP","HLA-DQB1","CCR5","CST1","OPTN","FASLG","FPR3","UXS1","LAG3","ZFF683","IL15")
down3<-c("CELSR2","USP27X-AS1")
score3a<-tpm[tpm$geneSymbol%in%$up3,]
score3a<-as.data.frame(t(score3a))
score3a<-apply(score3a,1,mean)
score3b<-tpm[tpm$geneSymbol%in%$down3,]
score3b<-as.data.frame(t(score3b))
score3b<-apply(score3b,1,mean)
score3<-left_join(score3a,score3b,by="ID")
score3$IFN<-score3$score3a-score3$score3b

#proliferation
up4<-c("HLA-DQA1","CD8A","B4GALT5","R3HDM1","UBD","OPTN","RNF19A","YME1L1")
down4<-c("NT5C3B","JPH1","ADIRF")
score4a<-tpm[tpm$geneSymbol%in%$up4,]
score4a<-as.data.frame(t(score4a))
score4a<-apply(score4a,1,mean)
score4b<-tpm[tpm$geneSymbol%in%$down4,]
score4b<-as.data.frame(t(score4b))
score4b<-apply(score4b,1,mean)
score4<-left_join(score4a,score4b,by="ID")
score4$proliferation<-score4$score4a-score4$score4b

#apoptosis
up5<-c("GPA33","CD8A","B4GALT5","HLA-DQB1","MARCO","OPTN","KLHDC4","CD70","GZMB","PLCG1","HSPA8","AKIRIN1","JSRP1")
down5<-c("AGR3","TRIM45","KLC4","MPC2","MED20")
score5a<-tpm[tpm$geneSymbol%in%$up5,]
score5a<-as.data.frame(t(score5a))
score5a<-apply(score5a,1,mean)
score5b<-tpm[tpm$geneSymbol%in%$down5,]
score5b<-as.data.frame(t(score5b))
score5b<-apply(score5b,1,mean)
score5<-left_join(score5a,score5b,by="ID")
score5$apoptosis<-score5$score5a-score5$score5b