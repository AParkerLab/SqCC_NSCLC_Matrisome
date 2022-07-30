#Pancancer analysis for both Lung ECM

#files downloaded from https://gdc.cancer.gov/about-data/publications/pancanatlas

# Upload data -------------------------------------------------------------
pancan_rnaseq<-read.table("./EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv", header = TRUE) #samples identified by full barcode
#RNAseq looks like FPKM or similar; not count data
pancan_RPPA<-read.table("./TCGA-RPPA-pancan-clean.txt", header = TRUE, sep = "\t")

library("readxl")
pancan_clin <- read_excel("./TCGA-CDR-SupplementalTableS1.xlsx", sheet = "TCGA-CDR") #imports as a tibble; patients identified by minimal barcode
pancan_clinEE <- read_excel("./TCGA-CDR-SupplementalTableS1.xlsx", sheet = "ExtraEndpoints") #imports as tibble

pancan_additionalclin<-read.table("./clinical_PANCAN_patient_with_followup.tsv", sep = "\t", header=TRUE, fill= TRUE)
#contains smoking behaviour etc

#extract minimal details from the additional clinical info
pancan_additionalclin_min<-pancan_additionalclin[,c(1,2,4,20,24,29,32:35, 37, 44:51, 56:57, 78:80,92,118, 119,179:196)]

# log2 transform expression data ----------------------------------------------------------
rownames(pancan_rnaseq)<-pancan_rnaseq[,1]

#reformat rownames so that they are HUGO format only
#note that there are two rows with same hugo id but dofferent entrez due to different isoforms
rownames(pancan_rnaseq)[16302]<-"SLC35E2B" #row16302 is SLC35E2|728661 corresponding to SLC35E2B
rownames(pancan_rnaseq)[16303]<-"SLC35E2A" #row16303 is SLC35E2|9906 corresponding to SLC35E2A
pancan_rnaseq<-pancan_rnaseq[-(1:29),]

#reformat row names
rownames(pancan_rnaseq)<-gsub("\\|.*", "", rownames(pancan_rnaseq))

#now remove the gene id column
pancan_rnaseq<-pancan_rnaseq[,-1]

#log2 transform the data with offset
pancan_rnaseqoff<-pancan_rnaseq+0.0001
pancan_rnaseql<-log2(pancan_rnaseqoff) 
# extract tumour only data for expression matrix ---------------------------
pancan_rnaseqlT<-pancan_rnaseql[,grep("01A", colnames(pancan_rnaseql))] #indicates 9572 Tumour cases
pancan_rnaseqlN<-pancan_rnaseql[,grep("11A", colnames(pancan_rnaseql))] #indicates 720 Tumour cases

#reformat column names so it's just the acc
colnames(pancan_rnaseqlT)<-gsub("\\.", "-", colnames(pancan_rnaseqlT))
colnames(pancan_rnaseqlT)<-substr(colnames(pancan_rnaseqlT), 1,12)

colnames(pancan_rnaseqlN)<-gsub("\\.", "-", colnames(pancan_rnaseqlN))
colnames(pancan_rnaseqlN)<-substr(colnames(pancan_rnaseqlN), 1,12)

#derive corresponding clinical data for these matrices
#Tonly expression matrix
pancan_clinTex<-pancan_clin[which(pancan_clin$bcr_patient_barcode %in% colnames(pancan_rnaseqlT)),]
pancan_clinTex<-pancan_clinTex[order(match(pancan_clinTex$bcr_patient_barcode, colnames(pancan_rnaseqlT))),]
#corresponding clinical data
pancan_rnaseqlTmin<-pancan_rnaseqlT[,which(colnames(pancan_rnaseqlT) %in% pancan_clinTex$bcr_patient_barcode)]
pancan_rnaseqlTmin<-pancan_rnaseqlTmin[,order(match(colnames(pancan_rnaseqlTmin),pancan_clinTex$bcr_patient_barcode))]
# zscale and define all samples (T+NT combined) ---------------------------

#define samples into their tumour types
allaccs<-substr(colnames(pancan_rnaseql),1,12)
allaccs<-gsub("\\.", "-", allaccs)

accdf<-data.frame(Accs = allaccs,
                  sample_barcode = colnames(pancan_rnaseql))
#match the accs to their cancer types
pancan_clinall<-merge(accdf, pancan_clin, by.x = "Accs", by.y = "bcr_patient_barcode", all = TRUE)

tumourtypeall<-unique(pancan_clinall$type)

#remove rows with no tumour type match
pancan_clinall2<-pancan_clinall[-which(is.na(pancan_clinall$type)==TRUE),]

#now redefine tumour types
tumourtypeall<-unique(pancan_clinall2$type)

alltumouraccs<-list()
for (i in 1:length(tumourtypeall)){
  alltumouraccs[[i]]<-pancan_clinall$sample_barcode[grep(tumourtypeall[i],pancan_clinall$type)]
}
names(alltumouraccs)<-tumourtypeall

#now iterate through the tumour types to z scale the tumour+non-tumour together
#zscale each tumour type individually across all Tumour samples
pancan_rnaseqlz<-t(pancan_rnaseql)
for (i in 1:length(alltumouraccs)){
  pancan_rnaseqlz[which(colnames(pancan_rnaseql) %in% alltumouraccs[[i]]),]<-apply(pancan_rnaseql[,which(colnames(pancan_rnaseql) %in% alltumouraccs[[i]])],1,scale) #because scaling by rows transposes the matrix
}

#apply the LUSCrisk score

#LUSC models given by 
#least stringent model
LUSC_lasso1<-read.csv("./LUSC_TvNT_DEGcoremat_lassomodelalphahalflambda1se_univariatecoef.csv",header = TRUE, row.names = 1)
#more stringent model
LUSC_lasso2<-read.csv("./LUSC_TvNT_DEGcoremat_lassomodelalpha1lambda1se_univariatecoef.csv", header= TRUE, row.names = 1)

#first create a submatrix 
pancan_rnaseqlz_lasso1<-pancan_rnaseqlz[,which(colnames(pancan_rnaseqlz) %in% rownames(LUSC_lasso1))]
pancan_rnaseqlz_lasso1<-pancan_rnaseqlz_lasso1[,order(match(colnames(pancan_rnaseqlz_lasso1), rownames(LUSC_lasso1)))]


riskscore<-matrix(NA, nrow = nrow(pancan_rnaseqlz_lasso1), ncol =1)
for (i in 1:length(tumourtypeall)){
  riskscore[which(rownames(pancan_rnaseqlz_lasso1) %in% alltumouraccs[[i]])]<-LRscore(LUSC_lasso1[-1,], pancan_rnaseqlz_lasso1[which(rownames(pancan_rnaseqlz_lasso1) %in% alltumouraccs[[i]]),])[[2]] #extract the actual score not the matrix itself
}

#append the score to the original matrix
pancan_LUSCriskscore<-data.frame(Accs = rownames(pancan_rnaseqlz_lasso1), 
                                 LUSCRiskScore = riskscore)

#create a TvNT column to indicate tumour status
TvNT<-rep(NA, nrow(pancan_LUSCriskscore))
TvNT[grep(".11[[:upper:]].", substring(pancan_LUSCriskscore$Acc, 1,17))]<-0
TvNT[grep(".0[[:digit:]][[:upper:]]", substring(pancan_LUSCriskscore$Acc, 1,17))] <-1
#create list for subtypes of solid tumours
solidtumourtype<-substring(pancan_LUSCriskscore$Acc, 14,15) 
#note that 01 = solid tumour, 02 = recurrent solid tumour, 03 = primary blood derived cancer - peripheral blood, 
#05 = additional new primary, 06 = metastatic tumour, 07 = additional metastatic, 11 = solid tumour normal

#add TvNT and tumour type data to the dataframe
pancan_LUSCriskscore$TvNT = as.factor(TvNT)
pancan_LUSCriskscore$solidtissuetype = solidtumourtype

#ROC analysis for each tumour type
library(pROC)
# Compute roc
pancan_LUSCriskscore.model.roc<-vector(mode = "list", length = length(tumourtypeall))
pancan_LUSCriskscore.model.roc.data<-vector(mode = "list", length = length(tumourtypeall))
pancan_LUSCriskscore.model.roc.aucdata<-vector(mode = "list", length = length(tumourtypeall))

names(pancan_LUSCriskscore.model.roc)<-tumourtypeall
names(pancan_LUSCriskscore.model.roc.data)<-tumourtypeall
names(pancan_LUSCriskscore.model.roc.aucdata)<-tumourtypeall

for (i in 1: length(tumourtypeall)){
  
  if (length(levels(droplevels(pancan_LUSCriskscore$TvNT[which(pancan_LUSCriskscore$Accs %in% alltumouraccs[[i]])])))>1){
  pancan_LUSCriskscore.model.roc[[i]]<- roc(pancan_LUSCriskscore$TvNT[which(pancan_LUSCriskscore$Accs %in% alltumouraccs[[i]])], pancan_LUSCriskscore[which(pancan_LUSCriskscore$Accs %in% alltumouraccs[[i]]),2])#second column is the LUSCRiskScore
  
  pancan_LUSCriskscore.model.roc.data[[i]] <- data.frame(thresholds = pancan_LUSCriskscore.model.roc[[i]]$thresholds,
                                                         sensitivity = pancan_LUSCriskscore.model.roc[[i]]$sensitivities,
                                                         specificity = pancan_LUSCriskscore.model.roc[[i]]$specificities)
  pancan_LUSCriskscore.model.roc.aucdata[[i]]<-auc(pancan_LUSCriskscore.model.roc[[i]])
  }
}

pancan_LUSCriskscore.model.rocmin<-pancan_LUSCriskscore.model.roc[!sapply(pancan_LUSCriskscore.model.roc,is.null)]
pancan_LUSCriskscore.model.rocmin<-pancan_LUSCriskscore.model.rocmin[-2]
length(pancan_LUSCriskscore.model.rocmin)

for (i in 1:length(pancan_LUSCriskscore.model.rocmin)){
plot.roc(pancan_LUSCriskscore.model.rocmin[[i]], print.auc = TRUE, add = FALSE, reuse.auc = FALSE)
title(main = names(pancan_LUSCriskscore.model.rocmin)[i])
}
dev.off() 

#plot the AUC as lollipop plot
#create dataframe just of tumour name and AUC
LUSC_LUSCriskscore_plottingdf<-data.frame(Tumour = names(unlist(pancan_LUSCriskscore.model.roc.aucdata)),
                                          AUC = unlist(pancan_LUSCriskscore.model.roc.aucdata))

LUSC_LUSCriskscore_plottingdf<-LUSC_LUSCriskscore_plottingdf[order(LUSC_LUSCriskscore_plottingdf$AUC, decreasing = TRUE),]

#remove SKCM from the plotting df as there is only one NT sample in that cohort
LUSC_LUSCriskscore_plottingdf<-LUSC_LUSCriskscore_plottingdf[-which(rownames(LUSC_LUSCriskscore_plottingdf)=="SKCM"),]

pancan_LUSC_riskscore_AUC_lollipop<-ggplot(LUSC_LUSCriskscore_plottingdf, aes(x=reorder(Tumour, -AUC), y=AUC)) +
  geom_point(color=ifelse(LUSC_LUSCriskscore_plottingdf$Tumour %in% "LUSC", "red2", "black"), 
             size=ifelse(LUSC_LUSCriskscore_plottingdf$Tumour %in% "LUSC", 3, 2)) + 
  xlab("Tumour Type") +
  geom_segment( aes(x=Tumour, xend=Tumour, y=0, yend=AUC),
                color=ifelse(LUSC_LUSCriskscore_plottingdf$Tumour %in% "LUSC", "red2", "black"), 
                size=ifelse(LUSC_LUSCriskscore_plottingdf$Tumour %in% "LUSC", 1.3, 0.7))+
  theme_bw()+
  theme(axis.text.x=element_text(angle = 90, size = 10, face="bold", vjust = 0.5 ),
        axis.text.y=element_text(size = 10, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 12, face="bold", hjust= 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey80"),
        panel.grid.minor= element_blank(),
        legend.position = "none")

#to create the violin plot we need the data in long format
pancan_LUSCriskscore2<-merge(pancan_LUSCriskscore, pancan_clinall2,by.x = "Accs", by.y ="sample_barcode")
#keep only tumour types that have tumour and non-tumour
pancan_LUSCriskscore2TvNT<-pancan_LUSCriskscore2[which(pancan_LUSCriskscore2$type %in% names(unlist(pancan_LUSCriskscore.model.roc.aucdata))),]
#remove SKCM from list since it only has one NT sample
pancan_LUSCriskscore2TvNT<-pancan_LUSCriskscore2TvNT[-which(pancan_LUSCriskscore2TvNT$type=="SKCM"),]

#test if significantly different for each tumour type
pancan_LUSCriskscore2TvNTsplit<-split(pancan_LUSCriskscore2TvNT, pancan_LUSCriskscore2TvNT$type)

test<-wilcox.test(LUSCRiskScore~TvNT,pancan_LUSCriskscore2TvNTsplit[[1]])$p.value
pancan_LUSCriskscore2TvNTsplit_wilcoxpval<-lapply(pancan_LUSCriskscore2TvNTsplit, function(x) wilcox.test(LUSCRiskScore~TvNT,x)$p.value)
pancan_LUSCriskscore2TvNTsplit_wilcoxpvaladj<-p.adjust(unlist(pancan_LUSCriskscore2TvNTsplit_wilcoxpval), method = "BH")

for (i in 1: length(tumourtypeall)){
  pancan_LUSCriskscore2TvNT
  if (length(levels(droplevels(pancan_LUSCriskscore$TvNT[which(pancan_LUSCriskscore$Accs %in% alltumouraccs[[i]])])))>1){
    pancan_LUSCriskscore.model.roc[[i]]<- roc(pancan_LUSCriskscore$TvNT[which(pancan_LUSCriskscore$Accs %in% alltumouraccs[[i]])], pancan_LUSCriskscore[which(pancan_LUSCriskscore$Accs %in% alltumouraccs[[i]]),2])#second column is the LUSCRiskScore
    
    pancan_LUSCriskscore.model.roc.data[[i]] <- data.frame(thresholds = pancan_LUSCriskscore.model.roc[[i]]$thresholds,
                                                           sensitivity = pancan_LUSCriskscore.model.roc[[i]]$sensitivities,
                                                           specificity = pancan_LUSCriskscore.model.roc[[i]]$specificities)
    pancan_LUSCriskscore.model.roc.aucdata[[i]]<-auc(pancan_LUSCriskscore.model.roc[[i]])
  }
}

#reorder so it matches the order of the ROC curve
pancan_LUSCriskscore2TvNT$type<-factor(pancan_LUSCriskscore2TvNT$type, levels = LUSC_LUSCriskscore_plottingdf$Tumour)

pancan_LUSCriskscore_vio<-ggplot(pancan_LUSCriskscore2TvNT, aes(x=type , y=LUSCRiskScore, fill = as.factor(TvNT)))+
  theme_bw()+
  geom_violin(trim=FALSE, position = position_dodge(width = 0.5), width = 3, size = 0.1) +
  scale_y_continuous(name = "LUSC ECM Risk Score") +
  scale_x_discrete(name = "Tumour Type")+
  scale_fill_manual(values = c("#F8766D", "#00BFC4"),name = "", limits = c("0","1"), labels = c("NT", "T"))+ #for default colours
  geom_boxplot(lwd = 0.1,width = 0.1, position = position_dodge(width = 0.5), outlier.shape=NA)+
  theme(axis.text.x=element_text(angle = 90, size = 10, face="bold", vjust = 0.5 ),
        axis.text.y=element_text(size = 10, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 12, face="bold", hjust= 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor= element_blank(),
        legend.position = "top")

ggplot_build(pancan_LUSCriskscore_vio)$data #to get the colour palette used

#combine the violin plot and ROC plot together
library(ggpubr)
Pancan_LUSCriskscore_plots<-ggarrange(pancan_LUSCriskscore_vio,
                            pancan_LUSC_riskscore_AUC_lollipop, 
                            ncol = 1, nrow = 2)

# Define which samples belong to which tumour types -----------------------

tumourtypes<-unique(pancan_clinTex$type)

tumourrows<-list()
tumouraccs<-list()
for (i in 1: length(tumourtypes)){
  tumourrows[[i]]<-grep(tumourtypes[i],pancan_clinTex$type)
  tumouraccs[[i]]<-pancan_clinTex$bcr_patient_barcode[tumourrows[[i]]]
  
  names(tumourrows)[i]<-tumourtypes[i]
  names(tumouraccs)[i]<-tumourtypes[i]
}

# Assign Molecular Subtypes -----------------------------------------------
pancan_rnaseqlTminz_clust<-pancan_rnaseqlTminz[,which(colnames(pancan_rnaseqlTminz) %in% rownames(LUSC_Km_centroid))]
pancan_rnaseqlTminz_clust<-pancan_rnaseqlTminz_clust[, order(match(colnames(pancan_rnaseqlTminz_clust), rownames(LUSC_Km_centroid)))]

#have checked that all genes in the centroids are represented in the expression matrix
#for each sample calculate the distances to each centroid
#apply the matrix function as euclidean distance
pancan_rnaseqlTminz_LUSCclust_euc<-calc_mat2mat_dist(pancan_rnaseqlTminz_clust, t(LUSC_Km_centroid)) #note one cluster per row for centroid matrix

#relabel the cluster assignments according to the column names
pancan_rnaseqlTminz_LUSCclust_euc<-as.data.frame(pancan_rnaseqlTminz_LUSCclust_euc)
pancan_rnaseqlTminz_LUSCclust_euc[,ncol(pancan_rnaseqlTminz_LUSCclust_euc)]<-as.factor(pancan_rnaseqlTminz_LUSCclust_euc[,ncol(pancan_rnaseqlTminz_LUSCclust_euc)])
levels(pancan_rnaseqlTminz_LUSCclust_euc[,ncol(pancan_rnaseqlTminz_LUSCclust_euc)])<-colnames(pancan_rnaseqlTminz_LUSCclust_euc)[1:3]

#make the column names informative for the distance method only
colnames(pancan_rnaseqlTminz_LUSCclust_euc)<-c(paste0(colnames(pancan_rnaseqlTminz_LUSCclust_euc)[1:3], "_Dist"), "Distmethod_cluster")

#combine the assignments together
pancan_rnaseqlTminz_LUSCclust<-data.frame(pancan_rnaseqlTminz_LUSCclust_euc)

#merge with the clinical dataframe
pancan_rnaseqlTminz_survsig_scoredf3<-merge(pancan_rnaseqlTminz_survsig_scoredf3, pancan_rnaseqlTminz_LUSCclust, by.x = "Accs", by.y = "row.names")
