#multiplex immunofluorescence  based functional dimension evaluation (counted by pathologists)

#imaging mass cytometry based functional dimension evaluation(protein)
#determination of the cut off values for CD8 positive cells
#sce1 is the processed single cell matrix (protein)
#the cutoff value used for determination of CD8 positive cell is 0.35
library(CATALYST)
sce1 <- estCutoffs(sce1)
metadata(sce1)$sep_cutoffs
sce1 <- AddMetaData(sce1, metadata = FetchData(sce1, vars = "CD8a"), col.name = "CD8A")
sce1a<-subset(sce1,subset =CD8A >0.35)

#determination of the cut off values for GZMB, MKI67 and CASP3 positive cells
#the cutoff values used for determination of GZMB, MKI67 and CASP3 positive cell are 0.85, 0.15, and 0.05, respectively
sce1a <- estCutoffs(sce1a)
metadata(sce1a)$sep_cutoffs
sce1a <- AddMetaData(sce1a, metadata = FetchData(sce1a, vars = "Granzyme_B"), col.name = "GB")
sce1a <- AddMetaData(sce1a, metadata = FetchData(sce1a, vars = "Ki-67"), col.name = "Ki67")
sce1a <- AddMetaData(sce1a, metadata = FetchData(sce1a, vars = "Cleaved_CP"), col.name = "CAPS3")
meta1<-sce1a@meta.data
meta1$GB<-ifelse(meta1$GB>0.85,1,0)
meta1$Ki67<-ifelse(meta1$Ki67>0.15,1,0)
meta1$CAPS3<-ifelse(meta1$CAPS3>0.05,1,0)

#abundance of CD8+ T cell expressing specific marker
#total is the total cell counts for each participant
assignment1<-aggregate(meta1$GB,by=list(meta1$orig.ident),FUN=sum,na.rm=TRUE)
assignment2<-aggregate(meta1$Ki67,by=list(meta1$orig.ident),FUN=sum,na.rm=TRUE)
assignment3<-aggregate(meta1$CAPS3,by=list(meta1$orig.ident),FUN=sum,na.rm=TRUE)
total<-left_join(total,assignment1,by="orig.ident")
total<-left_join(total,assignment2,by="orig.ident")
total<-left_join(total,assignment3,by="orig.ident")
total$cytotoxicity<-total$GB/total$total
total$proliferation<-total$Ki67/total$total
total$apoptosis<-total$CAPS3/total$total

#imaging mass cytometry based functional dimension evaluation(RNAscope)
#determination of the cut off values for CD8 positive cells
#sce2 is the processed single cell matrix (RNAscope)
#the cutoff value used for determination of CD8 positive cell is 0.36
sce2 <- estCutoffs(sce2)
metadata(sce2)$sep_cutoffs
sce2 <- AddMetaData(sce2, metadata = FetchData(sce1, vars = "CD8a"), col.name = "CD8A")
sce2a<-subset(sce2,subset =CD8A >0.6)

#determination of the cut off values for CXCL13 and CXCL9 positive cells
#the cutoff values used for determination of CXCL13 and CXCL9 positive cell are 1.2 and 0.25, respectively
sce2a <- estCutoffs(sce2a)
metadata(sce2a)$sep_cutoffs
sce2a <- AddMetaData(sce2a, metadata = FetchData(sce2a, vars = "CXCL13"), col.name = "CXCL13")
sce2a <- AddMetaData(sce2a, metadata = FetchData(sce2a, vars = "CXCL9"), col.name = "CXCL9")
meta2<-sce2a@meta.data
meta2$CXCL13<-ifelse(meta2$CXCL13>1.2,1,0)
meta2$CXCL9<-ifelse(meta2$CXCL9>0.25,1,0)

#abundance of CD8+ T cell expressing specific marker
#total is the total cell counts for each participant
assignment1<-aggregate(meta2$CXCL13,by=list(meta2$orig.ident),FUN=sum,na.rm=TRUE)
assignment2<-aggregate(meta2$CXCL9,by=list(meta2$orig.ident),FUN=sum,na.rm=TRUE)
total<-left_join(total,assignment1,by="orig.ident")
total<-left_join(total,assignment2,by="orig.ident")
total$reactivity<-total$CXCL13/total$total
total$IFNG<-total$CXCL9/total$total
