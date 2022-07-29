# Program for assessing TCGA lung data for matrisome influence

#clear all existing data
rm(list=ls())

#clear all plots
dev.off()

#packages

#for random forest/machine learning
library(tidyverse)
library(caret) #for machine learning workflows
library(randomForest)
library(gower)
library(ModelMetrics)
library(survminer)
library(survival)
library(edgeR)
library(maftools)
library(ggbiplot)
library(GSVA)
library(cluster) #silhouette function
library(fpc) #cluser.stats function

#load data from desktop project

#base packages
library(ggplot2)

# GDAC Clinical Data -read and reformat -----------------------------------

#read in clinical data
LUSC_clin <- read.table('./LUSC.clin.merged.txt',
                        header=F,row.names=1,sep='\t', stringsAsFactors = TRUE, fill = TRUE)

LUSC_clint<-transpose(as.data.frame(LUSC_clin))
colnames(LUSC_clint)<-rownames(LUSC_clin)
LUSCbarcode<-LUSC_clint$patient.bcr_patient_barcode
LUSCbarcode<-gsub("-",".", LUSCbarcode, fixed = TRUE)
LUSCbarcode<- gsub("tcga","TCGA",LUSCbarcode)
LUSC_clint<-data.frame(LUSCbarcode, LUSC_clint)
colnames(LUSC_clint)<-c("Acc", rownames(LUSC_clin))

#extract minimal clinical information
##LUSC
LUSC_OS<-data.frame(Acc = LUSC_clint$Acc,
                       Stage = LUSC_clint$patient.stage_event.pathologic_stage,
                       Age = as.numeric(LUSC_clint$patient.age_at_initial_pathologic_diagnosis),
                       Gender = as.factor(LUSC_clint$patient.gender),
                       Packyrs = as.numeric(LUSC_clint$patient.number_pack_years_smoked),
                       Smokingstatus = as.factor(LUSC_clint$patient.tobacco_smoking_history),
                       Race = as.factor(LUSC_clint$patient.race_list.race),
                       Fev1_fvc = as.numeric(LUSC_clint$patient.pre_bronchodilator_fev1_fvc_percent),
                       Fev1 = as.numeric(LUSC_clint$patient.pre_bronchodilator_fev1_percent),
                       dlco = as.numeric(LUSC_clint$patient.dlco_predictive_percent),
                       AnatomicalLocation = LUSC_clint$patient.location_in_lung_parenchyma,
                       Recurrence = LUSC_clint$patient.new_tumor_events.new_tumor_event.new_neoplasm_event_types.new_neoplasm_event_type,
                       Recurrence2 = LUSC_clint$'patient.new_tumor_events.new_tumor_event.new_neoplasm_event_types.new_neoplasm_event_type-2'
)

#reformat the Accs and stage
LUSC_OS$Stage<-substring(LUSC_OS$Stage, 7, 100000L)
LUSC_OS[,1:2]<-apply(LUSC_OS[,1:2], 2, toupper)
LUSC_OS[,1]<-paste0(LUSC_OS[,1], "_", "1")
rownames(LUSC_OS)<-LUSC_OS[,1]

#Survival Information
LUSC_dead<-which(is.na(LUSC_clint$patient.days_to_last_followup))
LUSC_timetoevent<-as.numeric(LUSC_clint$patient.days_to_last_followup)
LUSC_timetoevent[LUSC_dead]<-as.numeric(LUSC_clint$patient.days_to_death[LUSC_dead])

#limit the days to followup or death to five years
rep<-which(LUSC_timetoevent>1826.25)
LUSC_timetoevent5yr<-LUSC_timetoevent
LUSC_timetoevent5yr[rep]<-1826.25

#generate a column specifying alive(0) and dead(1)
LUSC_alivestatus<-which(LUSC_clint$patient.follow_ups.follow_up.vital_status=="alive")
LUSC_deadstatus<-which(LUSC_clint$patient.follow_ups.follow_up.vital_status=="dead")
LUSC_status<-matrix(NA, ncol = 1, nrow = nrow(LUSC_clint))
LUSC_status[LUSC_alivestatus]<-0
LUSC_status[LUSC_deadstatus]<-1

#Specify survival for specific time periods
LUSC_status5yrs<-LUSC_status
LUSC_status5yrs[rep]<-0

LUSC_OS<-data.frame(LUSC_OS,
                    Timetoevent_allyrs = as.numeric(LUSC_timetoevent),
                    Timetoevent_5yrs = as.numeric(LUSC_timetoevent5yr),
                    vital_status_numeric = as.numeric(LUSC_status),
                    vital_status5yrs_numeric = as.numeric(LUSC_status5yrs))
LUSC_vital_status1yr<-LUSC_OS$vital_status_numeric
LUSC_Timetoevent_1yr<-LUSC_OS$Timetoevent_allyrs
LUSC_vital_status2yr<-LUSC_OS$vital_status_numeric
LUSC_Timetoevent_2yr<-LUSC_OS$Timetoevent_allyrs

LUSC_post1yrevents<-which(LUSC_OS$Timetoevent_allyrs>365)
LUSC_post2yrevents<-which(LUSC_OS$Timetoevent_allyrs>730)

LUSC_vital_status1yr[LUSC_post1yrevents]<-0
LUSC_Timetoevent_1yr[LUSC_post1yrevents]<-365
LUSC_vital_status2yr[LUSC_post2yrevents]<-0
LUSC_Timetoevent_2yr[LUSC_post2yrevents]<-730

LUSC_OS<-data.frame(LUSC_OS,
                    Timetoevent_1yr = as.numeric(LUSC_Timetoevent_1yr),
                    Timetoevent_2yrs = as.numeric(LUSC_Timetoevent_2yr),
                    vital_status1yr_numeric = as.numeric(LUSC_vital_status1yr),
                    vital_status2yrs_numeric = as.numeric(LUSC_vital_status2yr))

# GDAC raw gene expression counts - read and reformat ---------------------

LUSC_rnaseq <- read.table('/home/z3403160/LungTCGA/CodeTest/LUSC.rnaseq__illuminahiseq_rnaseq__unc_edu__Level_3__gene_expression__data.data.txt',#nrows=20533,
                          header=T,row.names=1,sep='\t', stringsAsFactors = FALSE)

#extract raw count data
counts_cols<-which(LUSC_rnaseq[1,]=='raw_counts')
LUSC_rnaseq_raw<-LUSC_rnaseq[2:nrow(LUSC_rnaseq),counts_cols]
LUSC_rnaseq_rawnum<-apply(LUSC_rnaseq_raw,2, as.numeric)
rownames(LUSC_rnaseq_rawnum)<-rownames(LUSC_rnaseq_raw)

# Process gene expression counts ------------------------------------------
#Normalize and filter the data using EdgeR
library("edgeR")
LUSC_rnaseq_DGElist<-DGEList(counts = LUSC_rnaseq_rawnum)

#Filter out genes with less than 10 counts i.e. keep genes with more than 10 counts in at least two samples
#library sizes are approximately 60 million so cpm=0.6 corresponds to 10 counts per sample
LUSC_keep<-rowSums(cpm(LUSC_rnaseq_DGElist)>0.6) >=2
LUSC_rnaseq_DGElistf<-LUSC_rnaseq_DGElist[LUSC_keep, , keep.lib.sizes = FALSE]
#TMM Normalisation and log cpm transformation
LUSC_rnaseq_DGElistfTMM<-calcNormFactors(LUSC_rnaseq_DGElistf)
LUSC_rnaseq_logcpm<-cpm(LUSC_rnaseq_DGElistfTMM, log = TRUE, prior.count = 3)
rownames(LUSC_rnaseq_logcpm)<-gsub('\\|.*', '', rownames(LUSC_rnaseq_logcpm))

# subset Tumor vs Normal samples ------------------------------------------

#Subset the tumor and non-tumor values
LUSC_n_index <- which(substr(colnames(LUSC_rnaseq_logcpm),14,14) == '1')
LUSC_t_index <- which(substr(colnames(LUSC_rnaseq_logcpm),14,14) == '0')

LUSC_TNT<-matrix(0, nrow=1, ncol=ncol(LUSC_rnaseq_logcpm))
LUSC_TNT[LUSC_t_index]<-1
LUSC_rnaseq_logcpmTNT<-rbind(LUSC_TNT, LUSC_rnaseq_logcpm)
rownames(LUSC_rnaseq_logcpmTNT)[1]="TvNT"

rownames(LUSC_rnaseq_DGElistfTMM$counts)<-gsub('\\|.*', '', rownames(LUSC_rnaseq_DGElistfTMM$counts))

#reformat Accs for cibersort processing
LUSC_temp<-LUSC_rnaseq_logcpmTNT
colnames(LUSC_temp)<-paste0(colnames(LUSC_temp), "_", LUSC_temp[1,])
LUSC_temp<-data.frame(Gene = rownames(LUSC_temp),
                      LUSC_temp)

#set the group info of the DGEList object to Tumour vs Non-tumour
LUSC_TvNTgrp<-rep(0, length = nrow(LUSC_rnaseq_DGElistfTMM$samples))
LUSC_TvNTgrp[which(substr(rownames(LUSC_rnaseq_DGElistfTMM$samples),14,14) == '0')]<-1
LUSC_rnaseq_DGElistfTMM$samples$group<-as.factor(LUSC_TvNTgrp)

LUSC_rnaseq_logcpmTt<-t(LUSC_rnaseq_logcpmTNT[,which(LUSC_rnaseq_logcpmTNT[1,]==1)])

# Matrisome protein Databases -------------------------------------------------
#read in Naba matrisome gene list from MSigDb
matrisome_genelist<- read.table('./NabaMatrisome.txt',
           header=T,row.names=NULL,sep='\t', stringsAsFactors = FALSE)
matrisome_genelist<-matrisome_genelist[-1,]

##Read in the Matrisome master list from the matrisome website maintained by Richard O'Hynes and Alex naba
matrisome_master<-read.csv(file = './matrisome_hs_masterlist.csv', header = TRUE)
matrisome_master_min<-matrisome_master[1:4]

#subset matrisome master list into core matrisome
core_genenames<-as.character(matrisome_master_min[which(matrisome_master_min$Division=="Core matrisome"), 3])


# Merge with Clinical Data --------------------------------------
#transpose the datasets so genes are in columns
LUSC_rnaseq_logcpmTNTt<-t(as.matrix(LUSC_rnaseq_logcpmTNT))
LUSC_rnaseq_logcpmTNTtz<-apply(LUSC_rnaseq_logcpmTNTt[,3:ncol(LUSC_rnaseq_logcpmTNTt)], 2, scale)
rownames(LUSC_rnaseq_logcpmTNTtz)<-rownames(LUSC_rnaseq_logcpmTNTt)

#append a column for the Accs barcodes
LUSC_rnaseq_logcpmTNTt<-data.frame(substr(colnames(LUSC_rnaseq_logcpmTNT),1,12),as.matrix(LUSC_rnaseq_logcpmTNTt))
colnames(LUSC_rnaseq_logcpmTNTt)<-c("Acc", rownames(LUSC_rnaseq_logcpmTNT))
rownames(LUSC_rnaseq_logcpmTNTt)<-colnames(LUSC_rnaseq_logcpmTNT)

LUSC_logcpmt_T<-LUSC_rnaseq_logcpmTNTt[LUSC_rnaseq_logcpmTNTt$TvNT==1,] #Tumor
LUSC_logcpmt_N<-LUSC_rnaseq_logcpmTNTt[LUSC_rnaseq_logcpmTNTt$TvNT==0,] #NonTumor

#merge the clinical data with the tumor only data
LUSC_clin_logcpmT<-merge(LUSC_logcpmt_T, LUSC_clint, by = "Acc")

# TvNT DGE -----------------------------------------
# #subset the matrix by the the all matrisome genes
targetgenes_rows_matrisome<-which(rownames(LUSC_rnaseq_logcpmTNT) %in% matrisome_genelist)
LUSC_rnaseq_logcpmTNT_matrisome<-LUSC_rnaseq_logcpmTNT[targetgenes_rows_matrisome,]
#define the design matrix as tumour vs non-tumour
TNT<-LUSC_rnaseq_logcpmTNT[1,]
tumor_values<-which(TNT==1)
nontumor_values<-which(TNT=="0")
TNT[tumor_values]<-"T"
TNT[nontumor_values]<-"NT"
design_LUSC_TNTv2<-model.matrix(~0+TNT)
colnames(design_LUSC_TNTv2)<-c("NT", "T")

#DGE for tumor vs nontumor for matrisomal genes
fit_LUSCmatrisome_TNT <- lmFit(LUSC_rnaseq_logcpmTNT_matrisome, design_LUSC_TNTv2)
LUSCmatrisome_TNT_contrast<-makeContrasts(coef1 =T- NT, levels = design_LUSC_TNTv2)
LUSCmatrisome_TNT_contrastfit<-contrasts.fit(fit_LUSCmatrisome_TNT, LUSCmatrisome_TNT_contrast)
LUSCmatrisome_TNT_contrastfit<-eBayes(LUSCmatrisome_TNT_contrastfit)
LUSCmatrisome_TNT_contrastfitTable<-topTable(LUSCmatrisome_TNT_contrastfit, n= 'Inf')

#use contrasts for tumor v nontumor for all genes
fit_LUSCallgenes_TNT <- lmFit(LUSC_rnaseq_logcpmTNT, design_LUSC_TNTv2)
LUSCallgenes_TNT_contrast<-makeContrasts(coef1 =T- NT, levels = design_LUSC_TNTv2)
LUSCallgenes_TNT_contrastfit<-contrasts.fit(fit_LUSCallgenes_TNT, LUSCallgenes_TNT_contrast)
LUSCallgenes_TNT_contrastfit<-eBayes(LUSCallgenes_TNT_contrastfit)
LUSCallgenes_TNT_contrastfitTable<-topTable(LUSCallgenes_TNT_contrastfit, n= 'Inf')

#define which genes are significantly different between tumor and nontumor
LUSC_TvNTDGEsigIDs<- LUSCallgenes_TNT_contrastfitTable[which(LUSCallgenes_TNT_contrastfitTable$adj.P.Val<0.05),1]
LUSCallgenes_TNT_contrastfitTablesig<-LUSCallgenes_TNT_contrastfitTable[which(LUSCallgenes_TNT_contrastfitTable$adj.P.Val<0.05),]

#subset by core matrisome only
LUSCallgenes_TNT_contrastfitTable_coremat<-LUSCallgenes_TNT_contrastfitTable[which(LUSCallgenes_TNT_contrastfitTable$ID %in% core_genenames),]
colnames(LUSCallgenes_TNT_contrastfitTable_coremat)[1]= "Gene.Symbol"
LUSCallgenes_TNT_contrastfitTable_coremat<-merge(LUSCallgenes_TNT_contrastfitTable_coremat, matrisome_master_min, by = "Gene.Symbol")
levels(LUSCallgenes_TNT_contrastfitTable_coremat$Category)[2:4] = c("Glycoproteins", "Regulators", "Affiliated")

LUSCallgenes_TNT_contrastfitTablesig_mat<-LUSCallgenes_TNT_contrastfitTablesig[which(LUSCallgenes_TNT_contrastfitTablesig$ID %in% matrisome_master_min$Gene.Symbol),]
LUSCallgenes_TNT_contrastfitTablesig_coremat<-LUSCallgenes_TNT_contrastfitTablesig[which(LUSCallgenes_TNT_contrastfitTablesig$ID %in% matrisome_master_min$Gene.Symbol[which(matrisome_master_min$Division=="Core matrisome")]),]
LUSCallgenes_TNT_contrastfitTablesig_matassoc<-LUSCallgenes_TNT_contrastfitTablesig[which(LUSCallgenes_TNT_contrastfitTablesig$ID %in% matrisome_master_min$Gene.Symbol[which(matrisome_master_min$Division=="Matrisome-associated")]),]

LUSCallgenes_TNT_contrastfitTablesig_corematcollagen<-LUSCallgenes_TNT_contrastfitTablesig_coremat[which(LUSCallgenes_TNT_contrastfitTablesig_coremat$ID %in% matrisome_master_min$Gene.Symbol[which(matrisome_master_min$Category=="Collagens")]),]
LUSCallgenes_TNT_contrastfitTablesig_corematglycoproteins<-LUSCallgenes_TNT_contrastfitTablesig_coremat[which(LUSCallgenes_TNT_contrastfitTablesig_coremat$ID %in% matrisome_master_min$Gene.Symbol[which(matrisome_master_min$Category=="ECM Glycoproteins")]),]
LUSCallgenes_TNT_contrastfitTablesig_corematproteoglycans<-LUSCallgenes_TNT_contrastfitTablesig_coremat[which(LUSCallgenes_TNT_contrastfitTablesig_coremat$ID %in% matrisome_master_min$Gene.Symbol[which(matrisome_master_min$Category=="Proteoglycans")]),]

LUSCallgenes_TNT_contrastfitTablesig_matregulators<-LUSCallgenes_TNT_contrastfitTablesig_mat[which(LUSCallgenes_TNT_contrastfitTablesig_mat$ID %in% matrisome_master_min$Gene.Symbol[which(matrisome_master_min$Category=="ECM Regulators")]),]
LUSCallgenes_TNT_contrastfitTablesig_mataffiliated<-LUSCallgenes_TNT_contrastfitTablesig_mat[which(LUSCallgenes_TNT_contrastfitTablesig_mat$ID %in% matrisome_master_min$Gene.Symbol[which(matrisome_master_min$Category=="ECM-affiliated Proteins")]),]
LUSCallgenes_TNT_contrastfitTablesig_matsecretedfactors<-LUSCallgenes_TNT_contrastfitTablesig_mat[which(LUSCallgenes_TNT_contrastfitTablesig_mat$ID %in% matrisome_master_min$Gene.Symbol[which(matrisome_master_min$Category=="Secreted Factors")]),]

LUSCallgenes_TNT_contrastfitTablesig_matmerge<-merge(LUSCallgenes_TNT_contrastfitTablesig, matrisome_master_min, by.x = "ID", by.y = "Gene.Symbol")

# Figure 1B PCA plots of core matrisomal genes -------------------------------------------
#core matrisome
LUSC_rnaseq_logcpmTNTt_coremat<-LUSC_rnaseq_logcpmTNTt[,c(2, which(colnames(LUSC_rnaseq_logcpmTNTt) %in% core_genenames))]

LUSC_rnaseq_logcpmTNTt_coremat.PCA<-prcomp(LUSC_rnaseq_logcpmTNTt_coremat[,2:ncol(LUSC_rnaseq_logcpmTNTt_coremat)],
                                           center = TRUE,
                                           scale = TRUE)
summary(LUSC_rnaseq_logcpmTNTt_coremat.PCA)

library(ggbiplot)
LUSC_rnaseq_logcpmTNTt_coremat[,1]<-as.factor(LUSC_rnaseq_logcpmTNTt_coremat[,1])
levels(LUSC_rnaseq_logcpmTNTt_coremat[,1])<-c("NT", "T")

LUSC_allcoremat_PCA<-ggbiplot(LUSC_rnaseq_logcpmTNTt_coremat.PCA,
                              groups = as.factor(LUSC_rnaseq_logcpmTNTt_coremat[,1]),
                              ellipse = TRUE,
                              var.axes= FALSE) +
                      theme_bw()

# Figure 1C Matrisomal Composition of DEGs---------------------------------------------------------------
#create a piechart of the all DEGs compared with matrisomal DEGs
allcategories<-c("All DEGs" = length(LUSC_TvNTDGEsigIDs),
                 "All Matrisomal" = nrow(LUSCallgenes_TNT_contrastfitTablesig_mat))

matcategories<-c("Matrisome-Associated" = nrow(LUSCallgenes_TNT_contrastfitTablesig_matassoc),
                 "Core Matrisome" = nrow(LUSCallgenes_TNT_contrastfitTablesig_coremat))
corematcategories<-c("Collagens" = nrow(LUSCallgenes_TNT_contrastfitTablesig_corematcollagen),
                 "Glycoproteins" = nrow(LUSCallgenes_TNT_contrastfitTablesig_corematglycoproteins),
                 "Proteoglycans" = nrow(LUSCallgenes_TNT_contrastfitTablesig_corematproteoglycans))

matassoccategories<-c("Regulators" = nrow(LUSCallgenes_TNT_contrastfitTablesig_matregulators),
                     "Affiliated" = nrow(LUSCallgenes_TNT_contrastfitTablesig_mataffiliated),
                     "Secreted Factors" = nrow(LUSCallgenes_TNT_contrastfitTablesig_matsecretedfactors))

#par(mfrow=c(2,2))
dev.off()
pdf(file = ".LUSC_DEGs_Piechart_AllCategories.pdf", width = 5, height = 5)
pie(allcategories,labels = names(allcategories), col = c("black", "red2"), border = "white")
dev.off()

dev.off()
pdf(file = "./LUSC_DEGs_Piechart_AllMatCategories.pdf", width = 5, height = 5)
pie(matcategories,labels = names(matcategories), col = c("purple", "red2"), border = "white")
dev.off()

dev.off()
pdf(file = "./LUSC_DEGs_Piechart_CoreMatCategories.pdf", width = 5, height = 5)
pie(corematcategories,labels = names(corematcategories), col = c("gold", "orange", "red2"), border = "white")
dev.off()

dev.off()
pdf(file = "./LUSC_DEGs_Piechart_MatAssocCategories.pdf", width = 5, height = 5)
pie(matassoccategories,labels = names(matassoccategories), col = c("navy", "royalblue", "purple"), border = "white")
dev.off()

# Fig 1D: correlation matrix for all genes ----------------------------------------
#use only the tumour tissue matrisomal gene expression
LUSC_logcpmt_matrisome<-LUSC_logcpmt_T[,which(colnames(LUSC_logcpmt_T) %in% matrisome_genelist)]

#generate the correlation data using the Hmisc package
library(Hmisc)

#only plot the significantly different genes
LUSCmatrisome_sig<-rownames(LUSCmatrisome_TNT_contrastfitTable)[LUSCmatrisome_TNT_contrastfitTable$adj.P.Val<0.05]
LUSCcorematrisome_sig<-LUSCmatrisome_sig[LUSCmatrisome_sig %in% as.character(core_genenames)]

#filter the expression matrix by columns of significantly different genes
LUSC_logcpmt_corematrisomesig<-LUSC_logcpmt_matrisome[,colnames(LUSC_logcpmt_matrisome) %in% LUSCcorematrisome_sig]

LUSC_T_corematsig_corr <- Hmisc::rcorr(as.matrix(LUSC_logcpmt_corematrisomesig), type = "spearman")

#set the values that are not significant as 0 so you can set them to white in the colour mapping
LUSC_T_corematsig_corr_notsig<-LUSC_T_corematsig_corr
LUSC_T_corematsig_corr_notsig$r[which(LUSC_T_corematsig_corr$P>0.05)]<-0

#Figure 1D: plot the heatmap
library(ComplexHeatmap)
col_fun <- circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red")) #all nonsig set to zero so will be white
LUSC_matsig_corr_complexheatmap<-Heatmap(LUSC_T_corematsig_corr_notsig$r,
                                         col = col_fun,
                                         name = "r",
                                         border = TRUE,
                                         column_names_gp = gpar(fontsize = 4),
                                         row_names_gp = gpar(fontsize = 4),
                                         show_row_names = TRUE,
                                         row_split = 4,
                                         column_split = 4)

#extract the gene clusters identified in the heatmap for transcription factor analysis
LUSC_corrorder<-rep(list(NA),4)
LUSC_corrorder[[1]]<-rownames(LUSC_T_corematsig_corr_notsig$r)[row_order(LUSC_matsig_corr_complexheatmap)[[1]]]
LUSC_corrorder[[2]]<-rownames(LUSC_T_corematsig_corr_notsig$r)[row_order(LUSC_matsig_corr_complexheatmap)[[2]]]
LUSC_corrorder[[3]]<-rownames(LUSC_T_corematsig_corr_notsig$r)[row_order(LUSC_matsig_corr_complexheatmap)[[3]]]
LUSC_corrorder[[4]]<-rownames(LUSC_T_corematsig_corr_notsig$r)[row_order(LUSC_matsig_corr_complexheatmap)[[4]]]

# #assign the genes to the clusters they are in
LUSC_matcorrelationcluster<-rep(list(NA),length(LUSC_corrorder) )
LUSC_matcorrelationcluster[[1]]<-rep(1, length = length(LUSC_corrorder[[1]]))
LUSC_matcorrelationcluster[[2]]<-rep(2, length = length(LUSC_corrorder[[2]]))
LUSC_matcorrelationcluster[[3]]<-rep(3, length = length(LUSC_corrorder[[3]]))
LUSC_matcorrelationcluster[[4]]<-rep(4, length = length(LUSC_corrorder[[4]]))

LUSC_corrorder<-unlist(LUSC_corrorder)
LUSC_matcorrelationcluster<-unlist(LUSC_matcorrelationcluster)

LUSC_corrorder_names<-data.frame(genes = LUSC_corrorder,
                                 correlationcluster =LUSC_matcorrelationcluster)

# Supp Fig 1B Overlapping Core matrisome DEGs between LUAD and LUSC -------------------
LUADLUSC_coremat_table<-merge(LUADmatrisome_TNT_contrastfitTable_core, LUSCallgenes_TNT_contrastfitTable_coremat, by = "Gene.Symbol")
rownames(LUADLUSC_coremat_table) <- LUADLUSC_coremat_table$Gene.Symbol

#split into lists for both, or only significant in one subtype
LUADLUSC_coremat_commonsig<-LUADLUSC_coremat_table[which(LUADLUSC_coremat_table$adj.P.Val.x<0.05 & LUADLUSC_coremat_table$adj.P.Val.y<0.05),]
LUADLUSC_coremat_LUADsigonly<-LUADLUSC_coremat_table[which(LUADLUSC_coremat_table$adj.P.Val.x<0.05 & LUADLUSC_coremat_table$adj.P.Val.y>=0.05),]
LUADLUSC_coremat_LUSCsigonly<-LUADLUSC_coremat_table[which(LUADLUSC_coremat_table$adj.P.Val.x>=0.05 & LUADLUSC_coremat_table$adj.P.Val.y<0.05),]

#append a column indicating which histological subtype it's significant in
LUADLUSC_coremat_commonsig$whichsig = rep("Both", length = nrow(LUADLUSC_coremat_commonsig))
LUADLUSC_coremat_LUADsigonly$whichsig = rep("LUAD", length = nrow(LUADLUSC_coremat_LUADsigonly))
LUADLUSC_coremat_LUSCsigonly$whichsig = rep("LUSC", length = nrow(LUADLUSC_coremat_LUSCsigonly))

#append a column to the big dataframe
LUADLUSC_coremat_table$whichsig<-rep(NA, length = nrow(LUADLUSC_coremat_table))
LUADLUSC_coremat_table$pointcolour<-rep("black", length = nrow(LUADLUSC_coremat_table))

LUADLUSC_coremat_table$whichsig[which(LUADLUSC_coremat_table$adj.P.Val.x<0.05 & LUADLUSC_coremat_table$adj.P.Val.y<0.05)]<-"Both"
LUADLUSC_coremat_table$whichsig[which(LUADLUSC_coremat_table$adj.P.Val.x<0.05 & LUADLUSC_coremat_table$adj.P.Val.y>=0.05)]<-"LUAD"
LUADLUSC_coremat_table$whichsig[which(LUADLUSC_coremat_table$adj.P.Val.x>=0.05 & LUADLUSC_coremat_table$adj.P.Val.y<0.05)]<-"LUSC"

LUADLUSC_coremat_table$pointcolour[LUADLUSC_coremat_table$whichsig=="LUSC"]<-"green"
LUADLUSC_coremat_table$pointcolour[LUADLUSC_coremat_table$whichsig=="LUAD"]<-"royalblue"
LUADLUSC_coremat_table$pointcolour[LUADLUSC_coremat_table$whichsig=="Both"]<-"red2"

#for only the significant list
LUADLUSC_coremat_commondiffsig<-rbind(LUADLUSC_coremat_commonsig,
                                      LUADLUSC_coremat_LUADsigonly,
                                      LUADLUSC_coremat_LUSCsigonly)
LUADLUSC_coremat_commondiffsig$colour<-rep("black", length = nrow(LUADLUSC_coremat_commondiffsig))
LUADLUSC_coremat_commondiffsig$colour[LUADLUSC_coremat_commondiffsig$whichsig=="LUAD"]<-"royalblue"
LUADLUSC_coremat_commondiffsig$colour[LUADLUSC_coremat_commondiffsig$whichsig=="LUSC"]<-"forestgreen"
LUADLUSC_coremat_commondiffsig$colour[LUADLUSC_coremat_commondiffsig$whichsig=="Both"]<-"red2"

col <-c( "LUAD" = "royalblue" , "LUSC" = "green", "Both" = "red2") 
LUADLUSC_coremat_commondiffsig_scatter<-ggplot(LUADLUSC_coremat_table, aes(x=logFC.y, y=logFC.x, color = LUADLUSC_coremat_table$whichsig))+
  geom_point(color=LUADLUSC_coremat_table$pointcolour)+
  theme_bw()+
  scale_y_continuous(name = "LUAD logFC") +
  scale_x_continuous(name = "LUSC logFC")+
  ggtitle("Core Matrisomal T vs NT")+
  geom_vline(xintercept = 0, colour = "gray50",linetype='dotted') + 
  geom_hline(yintercept = 0, colour = "gray50",linetype='dotted') +
  geom_abline(intercept =0 , slope = 1, colour = "gray50")+
  theme(axis.text=element_text(size = 8, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 12, face="bold", hjust= 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor= element_blank(),
        legend.position = "right")

# Supp Fig 1D Correlation with LUAD/LUSC Markers --------------------------------------
#TP63 as squamous marker
tp63_col<-which(colnames(LUSC_rnaseq_logcpmTt)=="TP63")

LUSC_tp63corr_val<-matrix(NA, ncol = ncol (LUSC_rnaseq_logcpmTt), nrow = 2)
for (i in 3: ncol(LUSC_rnaseq_logcpmTt)){
  LUSC_tp63corr<-cor.test(LUSC_rnaseq_logcpmTt[,i], LUSC_rnaseq_logcpmTt[,tp63_col], method = "spearman")
  LUSC_tp63corr_val[1,i]<-LUSC_tp63corr$estimate
  LUSC_tp63corr_val[2,i]<-LUSC_tp63corr$p.value
}
LUSC_tp63corr_val<-data.frame(t(LUSC_tp63corr_val),
                              pvaladj = p.adjust(LUSC_tp63corr_val[2,], method = "BH"))
LUSC_tp63corr_valt<-t(LUSC_tp63corr_val)
colnames(LUSC_tp63corr_valt) <-colnames(LUSC_rnaseq_logcpmTt) 
rownames(LUSC_tp63corr_valt)[1:2] <-c("r_spearman","pval")

#filter for matrisomal list
LUSC_tp63corr_valt_coremat<-LUSC_tp63corr_valt[,which(colnames(LUSC_tp63corr_valt) %in% core_genenames)]

#Supp Fig 1D- plot LUSC correlation between TP63 and COL4A6
col4a6_col<-which(colnames(LUSC_rnaseq_logcpmTt)=="COL4A6")

LUSC_TP63_COL4A6_scatter<-ggplot(as.data.frame(LUSC_rnaseq_logcpmTt[, c(tp63_col, col4a6_col)]), aes(x=TP63, y=COL4A6))+
  geom_point()+
  theme_bw()+
  scale_y_continuous(name = "Relative COL4A6 Expression") +
  scale_x_continuous(name = "Relative TP63 Expression")+
  ggtitle("LUSC TP63 Correlation")+
  theme(axis.text=element_text(size = 8, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 12, face="bold", hjust= 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor= element_blank())
LUSC_TP63_COL4A6_scatter

# Supp Fig 1K Transcription factor enrichment analysis --------------------------------
#read in data output from chEA3 web server
LUSC_chea_1<-read.table(file = "./TCGA_LUSC_CoreMatCluster1_chEA3_Integrated_meanRank.tsv", sep = "\t", header = TRUE)
LUSC_chea_2<-read.table(file = "./TCGA_LUSC_CoreMatCluster2_chEA3_Integrated_meanRank.tsv", sep = "\t", header = TRUE)
LUSC_chea_3<-read.table(file = "./TCGA_LUSC_CoreMatCluster3_chEA3_Integrated_meanRank.tsv", sep = "\t", header = TRUE)
LUSC_chea_4<-read.table(file = "./TCGA_LUSC_CoreMatCluster4_chEA3_Integrated_meanRank.tsv", sep = "\t", header = TRUE)

#create a column indicating the number of overlapping genes
library(stringr)
LUSC_chea_1$numberofoverlappinggenes<-str_count(LUSC_chea_1$Overlapping_Genes, ",") +1
LUSC_chea_2$numberofoverlappinggenes<-str_count(LUSC_chea_2$Overlapping_Genes, ",") +1
LUSC_chea_3$numberofoverlappinggenes<-str_count(LUSC_chea_3$Overlapping_Genes, ",") +1
LUSC_chea_4$numberofoverlappinggenes<-str_count(LUSC_chea_4$Overlapping_Genes, ",") +1

#convert the number of overlapping genes to a percentage of genes in teh list
cluster1length<-nrow(LUSC_corrorder_names[which(LUSC_corrorder_names$correlationcluster=="1"),])
cluster2length<-nrow(LUSC_corrorder_names[which(LUSC_corrorder_names$correlationcluster=="2"),])
cluster3length<-nrow(LUSC_corrorder_names[which(LUSC_corrorder_names$correlationcluster=="3"),])
cluster4length<-nrow(LUSC_corrorder_names[which(LUSC_corrorder_names$correlationcluster=="4"),])

LUSC_chea_1$ProportionOverlappingGenes<-100*LUSC_chea_1$numberofoverlappinggenes/cluster1length
LUSC_chea_2$ProportionOverlappingGenes<-100*LUSC_chea_2$numberofoverlappinggenes/cluster2length
LUSC_chea_3$ProportionOverlappingGenes<-100*LUSC_chea_3$numberofoverlappinggenes/cluster3length
LUSC_chea_4$ProportionOverlappingGenes<-100*LUSC_chea_4$numberofoverlappinggenes/cluster4length

#create dot plot of top TFs
LUSC_chea_1$CorrelationCluster<-rep(1, length = nrow(LUSC_chea_1))
LUSC_chea_2$CorrelationCluster<-rep(2, length = nrow(LUSC_chea_1))
LUSC_chea_3$CorrelationCluster<-rep(3, length = nrow(LUSC_chea_1))
LUSC_chea_4$CorrelationCluster<-rep(4, length = nrow(LUSC_chea_1))

#top 10 ranked TFs
LUSC_chea<-rbind(LUSC_chea_1[1:10,],LUSC_chea_2 [1:10,],LUSC_chea_3[1:10,],LUSC_chea_4[1:10,])

LUSC_chea$TForder<-seq(from = 1, to = nrow(LUSC_chea), by = 1)

LUSC_chea_bubbleplot<-ggplot(LUSC_chea,
                             aes(x = reorder(TF, TForder), 
                                 y = reorder(CorrelationCluster, -CorrelationCluster),
                                 colour = ProportionOverlappingGenes,
                                 size = 1/Score)) +
  geom_point() +
  labs(x = "Transcription Factor", y = "Correlation Cluster") +
  theme_bw()+
  theme(axis.text.x=element_text(size = 8, face="bold", angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y=element_text(size = 10, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 12, face="bold", hjust= 0.5),
        panel.grid.major = element_line(colour = "grey95", size = 0.5 ),
        panel.grid.minor= element_blank())


# Matrix Risk Signature Fig2ABD SuppFig2A ----------------------------------------------------
library(tidyverse)
library(caret)
library(glmnet)
library(pROC)
library(scales)

#LUSC
#reformat this data so it only contains DEG matrisomal genes
alpha<-0.05
LUSCmatrisome_TNT_contrastfitTablesig<-LUSCmatrisome_TNT_contrastfitTable[LUSCmatrisome_TNT_contrastfitTable$adj.P.Val<alpha,]

#Core matrisomal genes
LUSCmatrisome_TNT_contrastfitTablesigcore<-LUSCmatrisome_TNT_contrastfitTablesig[which(rownames(LUSCmatrisome_TNT_contrastfitTablesig) %in% as.character(core_genenames)),]

#Significantly DE core matrisomal genes
LUSC_rnaseq_logcpmTNTsigcore<-LUSC_rnaseq_logcpmTNT[rownames(LUSC_rnaseq_logcpmTNT) %in% rownames(LUSCmatrisome_TNT_contrastfitTablesigcore),]

#z-scale data
LUSC_rnaseq_logcpmTNTsigcorez<-apply(LUSC_rnaseq_logcpmTNTsigcore,1, scale)
rownames(LUSC_rnaseq_logcpmTNTsigcorez)<-colnames(LUSC_rnaseq_logcpmTNTsigcore)

#re-order the columns according to most signfiicantly DE genes
LUSC_rnaseq_logcpmTNTsigcorez_ord<-as.data.frame(LUSC_rnaseq_logcpmTNTsigcorez[,match(c("TvNT", rownames(LUSCmatrisome_TNT_contrastfitTablesigcore)),colnames(LUSC_rnaseq_logcpmTNTsigcorez))])
LUSC_rnaseq_logcpmTNTsigcorez_ord[,1]<-ifelse(substr(rownames(LUSC_rnaseq_logcpmTNTsigcorez_ord), 14,14)==1,0,1)
colnames(LUSC_rnaseq_logcpmTNTsigcorez_ord)[1]<-"TvNT"

#Bootstrap data into test and training data
set.seed(123)
LUSC_rnaseq_logcpmTNTsigcorez_ord$TvNT<-as.numeric(LUSC_rnaseq_logcpmTNTsigcorez_ord$TvNT)
LUSC_training.samples <- LUSC_rnaseq_logcpmTNTsigcorez_ord$TvNT %>% 
  createDataPartition(p = 0.8, list = FALSE)
LUSC_rnaseq_logcpmTNTsigcorez_ord_train.data  <- LUSC_rnaseq_logcpmTNTsigcorez_ord[LUSC_training.samples, ]
LUSC_rnaseq_logcpmTNTsigcorez_ord_test.data <- LUSC_rnaseq_logcpmTNTsigcorez_ord[-LUSC_training.samples, ]

#lasso penalised regression model

x <- model.matrix(TvNT~., LUSC_rnaseq_logcpmTNTsigcorez_ord_train.data )[,-1]
y<-LUSC_rnaseq_logcpmTNTsigcorez_ord_train.data$TvNT

#find the optimal value of lambda that minimises crossvalidation error
set.seed(123)
cv.lasso <- cv.glmnet(x, y, alpha = 1, family = "binomial")
plot(cv.lasso)
cv.lasso$lambda.min
cv.lasso$lambda.1se #value of lambda to give accurate model with smallest number of variables

#coefficients using minimal lambda
coef(cv.lasso, cv.lasso$lambda.min)
#coefficients using lambda within 1xstandard error of minimal lambda
coef(cv.lasso, cv.lasso$lambda.1se)

#compute the final model using minimal lambda
# Final model with lambda.min
LUSC_lasso.model <- glmnet(x, y, alpha = 1, family = "binomial",
                           lambda = cv.lasso$lambda.min)
coef(LUSC_lasso.model)
# Make prediction on test data
x.test <- model.matrix(TvNT ~., LUSC_rnaseq_logcpmTNTsigcorez_ord_test.data)[,-1]
probabilities <- LUSC_lasso.model %>% predict(newx = x.test, type = "response")
predicted.classes <- ifelse(probabilities > 0.5, 1, 0)
# Model accuracy
observed.classes <- LUSC_rnaseq_logcpmTNTsigcorez_ord_test.data$TvNT
mean(predicted.classes == observed.classes) #gives 100% accuracy

# Compute roc
LUSC_lasso.model.roc <- roc(LUSC_rnaseq_logcpmTNTsigcorez_ord_test.data$TvNT, as.numeric(probabilities))
LUSC_lasso.model.rocplot<-plot.roc(LUSC_lasso.model.roc, print.auc = TRUE)
LUSC_lasso.model.roc.data <- data.frame(
  thresholds = LUSC_lasso.model.roc$thresholds,
  sensitivity = LUSC_lasso.model.roc$sensitivities,
  specificity = LUSC_lasso.model.roc$specificities
)
LUSC_lasso.model.roc.auc<-auc(LUSC_lasso.model.roc)

##if alpha=0.5 is used
# Final model with lambda.min
set.seed(123)
cv.lasso <- cv.glmnet(x, y, alpha = 0.5, family = "binomial")
plot(cv.lasso)
cv.lasso$lambda.min
cv.lasso$lambda.1se #value of lambda to give accurate model with smallest number of variables

coef(cv.lasso, cv.lasso$lambda.1se)

LUSC_lasso.model3 <- glmnet(x, y, alpha = 0.5, family = "binomial",
                            lambda = cv.lasso$lambda.1se)
coef(LUSC_lasso.model3)
rownames(coef(LUSC_lasso.model3))[which(coef(LUSC_lasso.model3)!=0)] #21

# Make prediction on test data
x.test <- model.matrix(TvNT ~., LUSC_rnaseq_logcpmTNTsigcorez_ord_test.data)[,-1]
probabilities <- LUSC_lasso.model3 %>% predict(newx = x.test, type = "response")
predicted.classes <- ifelse(probabilities > 0.5, 1, 0)
# Model accuracy
observed.classes <- LUSC_rnaseq_logcpmTNTsigcorez_ord_test.data$TvNT
mean(predicted.classes == observed.classes) #also gives 100% accuracy

# Compute roc
LUSC_lasso.model3.roc <- roc(LUSC_rnaseq_logcpmTNTsigcorez_ord_test.data$TvNT, as.numeric(probabilities))
LUSC_lasso.model3.rocplot<-plot.roc(LUSC_lasso.model3.roc, print.auc = TRUE)
LUSC_lasso.model3.roc.data <- data.frame(
  thresholds = LUSC_lasso.model3.roc$thresholds,
  sensitivity = LUSC_lasso.model3.roc$sensitivities,
  specificity = LUSC_lasso.model3.roc$specificities
)
LUSC_lasso.model3.roc.auc<-auc(LUSC_lasso.model3.roc)

LUSC_lasso3genes<-rownames(coef(LUSC_lasso.model3))[which(coef(LUSC_lasso.model3)!=0)]
LUSC_lasso3genes<-LUSC_lasso3genes[-1]

#genes with non-zero coefficients
LUSC_lasso3genes<-c("MFAP4", "RSPO1", "SVEP1", "FBLN5","NTN4",
                    "WISP2", "VWF", "TNXB", "SLIT2", "MMRN2", 
                    "LAMB2", "CRIM1", "GLDN", "ABI3BP", "MMRN1",
                    "EMILIN2", "SPARCL1", "LGI4", "RSPO2", "PRG4",
                    "SPOCK2", "TNR", "CILP2", "COL10A1", "CTHRC1",
                    "SPP1", "COL11A1", "COL7A1")
#perform univariate ROC on each gene individually in the 28 gene signature defined by lasso3
LUSC_rnaseq_logcpmTNTsigcorez_ord_test.data_lasso3<-LUSC_rnaseq_logcpmTNTsigcorez_ord_test.data[,which(colnames(LUSC_rnaseq_logcpmTNTsigcorez_ord_test.data) %in% LUSC_lasso3genes)]

LUSC_lasso.model3.univar_roc.auc<-c()
for (i in 1: ncol(LUSC_rnaseq_logcpmTNTsigcorez_ord_test.data_lasso3)){
  LUSC_lasso.model3.univar_roc <- roc(LUSC_rnaseq_logcpmTNTsigcorez_ord_test.data$TvNT, LUSC_rnaseq_logcpmTNTsigcorez_ord_test.data_lasso3[,i])
  LUSC_lasso.model3.univar_roc.auc[i]<-LUSC_lasso.model3.univar_roc$auc
}

LUSC_lasso.model3.univar_roc_df<-data.frame(Gene = colnames(LUSC_rnaseq_logcpmTNTsigcorez_ord_test.data_lasso3),
                                            ROC.AUC = LUSC_lasso.model3.univar_roc.auc)

#subset the original matrix to include matris risk signature genes
LUSC_rnaseq_logcpmTNTsigcorez_lasso3<-LUSC_rnaseq_logcpmTNTsigcorez_ord[,c(1,which(colnames(LUSC_rnaseq_logcpmTNTsigcorez_ord) %in% LUSC_lasso3genes))]

#Calculate Odds Ratios for matrix risk signature genes
library(logistf)
LUSC_uniLR_ORCI_lasso3<-matrix(NA, nrow = ncol(LUSC_rnaseq_logcpmTNTsigcorez_lasso3), ncol = 3)
for (i in 2:ncol(LUSC_rnaseq_logcpmTNTsigcorez_lasso3)){
  LUSC_TvNT_lasso3glmfit <- logistf(formula= as.numeric(LUSC_rnaseq_logcpmTNTsigcorez_lasso3[,1]) ~LUSC_rnaseq_logcpmTNTsigcorez_lasso3[,i] ,
                                    data = LUSC_rnaseq_logcpmTNTsigcorez_lasso3,
                                    family = binomial(link = "logit"))
   LUSC_uniLR_ORCI_lasso3[i,]<-cbind(Coef = coef(LUSC_TvNT_lasso3glmfit)[2],
                                    LowerConfint = confint(LUSC_TvNT_lasso3glmfit)[2,1],
                                    UpperConfint = confint(LUSC_TvNT_lasso3glmfit)[2,2])
}
#convert coeff and ci to linear space
LUSC_uniLR_ORCI_lasso3<-data.frame(LUSC_uniLR_ORCI_lasso3,
                                   exp(LUSC_uniLR_ORCI_lasso3[,1:3]))
colnames(LUSC_uniLR_ORCI_lasso3)<-c("Coef", "LowerConfint", "UpperConfint",
                                    "ExpCoef", "ExpLowerConfint", "ExpUpperConfint")
rownames(LUSC_uniLR_ORCI_lasso3)<-colnames(LUSC_rnaseq_logcpmTNTsigcorez_lasso3)

LUSC_uniLR_ORCI_lasso3<-LUSC_uniLR_ORCI_lasso3[order(LUSC_uniLR_ORCI_lasso3$ExpCoef, decreasing = TRUE),]
LUSC_uniLR_ORCI_lasso3<-data.frame(LUSC_uniLR_ORCI_lasso3,
                                   Pos = seq(1,nrow(LUSC_uniLR_ORCI_lasso3), by = 1),
                                   Gene = rownames(LUSC_uniLR_ORCI_lasso3))
LUSC_uniLR_ORCI_lasso3$Gene<-factor(LUSC_uniLR_ORCI_lasso3$Gene, levels = LUSC_uniLR_ORCI_lasso3$Gene[order(LUSC_uniLR_ORCI_lasso3$Pos)])

LUSC_lasso3_OR <- ggplot(LUSC_uniLR_ORCI_lasso3[-nrow(LUSC_uniLR_ORCI_lasso3),], aes(x = ExpCoef, y = Gene)) + 
  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") + 
  geom_errorbarh(aes(xmax = ExpUpperConfint, xmin = ExpLowerConfint), size = .5, height = 
                   .2, color = "gray50") +
  geom_point(size = 3.5, color = "red2") +
  theme_bw()+
  theme(panel.grid.minor = element_blank()) +
  ylab("") +
  xlab("Odds Ratio") +
  ggtitle("LUSC Tumour Classification")
LUSC_lasso3_OR


LUSC_lasso3_ORlog <- ggplot(LUSC_uniLR_ORCI_lasso3[-nrow(LUSC_uniLR_ORCI_lasso3),], aes(x = ExpCoef, y = Gene)) + 
  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") + 
  geom_errorbarh(aes(xmax = ExpUpperConfint, xmin = ExpLowerConfint), size = .5, height = 
                   .2, color = "gray50") +
  geom_point(size = 3.5, color = "red2") +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "b")+
  theme_bw()+
  theme(panel.grid.minor = element_blank()) +
  ylab("") +
  xlab("Odds Ratio") +
  ggtitle("LUSC Tumour Classification")
LUSC_lasso3_ORlog


#Supplementary Figure 2A
#plot these as a heatmap to check that they distinguish T from NT
library(ComplexHeatmap)
colours <- list("Tissue"=c("NT"="darkgrey","T"="black"))
Tissue<-as.factor(LUSC_rnaseq_logcpmTNTsigcorez_lasso3[,1])
levels(Tissue)<-c("NT", "T")
LUSC_lasso_colAnn <- HeatmapAnnotation(Tissue = Tissue, which="col",
                                       col=colours,
                                       annotation_width=unit(c(1, 4), "cm"), gap=unit(1, "mm"), na_col = "white")

library(circlize)
col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")) 

LUSC_rnaseq_logcpmTNTsigcorez_lasso3v2<-LUSC_rnaseq_logcpmTNTsigcorez_lasso3[,-1]

#reorder columns
LUSC_rnaseq_logcpmTNTsigcorez_lasso3v2<-LUSC_rnaseq_logcpmTNTsigcorez[,which(colnames(LUSC_rnaseq_logcpmTNTsigcorez) %in% LUSC_lasso3genes)]
LUSC_rnaseq_logcpmTNTsigcorez_lasso3v2<-LUSC_rnaseq_logcpmTNTsigcorez_lasso3v2[,order(match(colnames(LUSC_rnaseq_logcpmTNTsigcorez_lasso3v2), LUSC_lasso3genes))]

Tissue<-as.factor(LUSC_rnaseq_logcpmTNT[1,])
levels(Tissue)<-c("NT", "T")
LUSC_lasso_colAnn <- HeatmapAnnotation(Tissue = Tissue, which="col",
                                       col=colours,
                                       annotation_width=unit(c(1, 4), "cm"), gap=unit(1, "mm"), na_col = "white")

LUSC_TvNT_lasso3_ccmap<-Heatmap(t(LUSC_rnaseq_logcpmTNTsigcorez_lasso3v2),
                                col = col_fun,
                                name = "z",
                                border = TRUE,#black border
                                cluster_rows = FALSE,
                                cluster_columns = TRUE, #this does not work
                                column_names_gp = gpar(fontsize = 0),
                                row_names_gp = gpar(fontsize = 8),
                                show_row_names = TRUE,
                                top_annotation=LUSC_lasso_colAnn)

#apply the matrix risk score
LUSC_uniLR_ORCI_lasso3_reord<-LUSC_uniLR_ORCI_lasso3[order(LUSC_uniLR_ORCI_lasso3$Coef, decreasing = FALSE),]

LRscore <- function(coefficients, zscoreexp){
  #genes in columns
  coef<-coefficients[,1]
  #zscoreexp<-LUSC_rnaseq_logcpmTNTsigcorez_lasso3_ord
  matrix<-matrix(NA, ncol = ncol(zscoreexp), nrow = nrow(zscoreexp))
  for (i in 1:ncol(zscoreexp)){
    for (j in 1:nrow(zscoreexp)){
      matrix[j,i]<-coef[i]*zscoreexp[j,i]
    }
  }
  score<-apply(matrix, 1, sum)
  return(list(matrix, score))
}

b<-LRscore(LUSC_uniLR_ORCI_lasso3_reord[-nrow(LUSC_uniLR_ORCI_lasso3_reord),], LUSC_rnaseq_logcpmTNTsigcorez_lasso3v2)
LUSC_rnaseq_logcpmTNTsigcorez_lasso3_scorematrixv2<-b[[1]]
LUSC_rnaseq_logcpmTNTsigcorez_lasso3_scorev2<-b[[2]]

matrixscore<-data.frame(Acc= substr(rownames(LUSC_rnaseq_logcpmTNTsigcorez_lasso3v2),1,12),
                        TvNT = ifelse(substr(rownames(LUSC_rnaseq_logcpmTNTsigcorez_lasso3v2),14,14)==1,0,1),
                        Score = LUSC_rnaseq_logcpmTNTsigcorez_lasso3_scorev2)
rownames(matrixscore)<-rownames(LUSC_rnaseq_logcpmTNTsigcorez_lasso3v2)

#Assess Matrix score association with stage
matrixscoreT<-matrixscore[grep("1", matrixscore$TvNT),]
matrixscoreT$Acc<-paste0(matrixscoreT$Acc,"_1", "")
matrixscore_clin<-merge(matrixscoreT, LUSC_OS, by.x = "Acc", by.y = "row.names")

matrixscore_clinv2<-matrixscore_clin[-which(matrixscore_clin$Stage=="II" | is.na(matrixscore_clin$Stage==TRUE)),]
#comparison by stage
LUSC_lasso3_matrixscore_stage_kruskal<-kruskal.test(matrixscore_clinv2$Score~as.factor(matrixscore_clinv2$Stage)) 
LUSC_lasso3_matrixscore_stage<-ggplot(matrixscore_clinv2, aes(x=as.factor(Stage), y=Score, fill = as.factor(Stage)))+
  theme_bw()+
  geom_violin(trim=FALSE) +
  scale_y_continuous(name = "Matrix Risk Score") +
  ggtitle("LUSC Matrix Score (a=0.5)")+
  geom_boxplot(width = 0.1)+
  xlab("Stage")+
  labs(fill = "Tumour and Non-Tumour")+
  theme(axis.text=element_text(size = 10, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 10, face="bold", hjust= 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor= element_blank(),
        legend.position = "none")
LUSC_lasso3_matrixscore_stage

#Fig 2B
LUSC_lasso3_matrixscore_wilcox<-wilcox.test(matrixscore$Score~as.factor(matrixscore$TvNT))
LUSC_lasso3_matrixscore<-ggplot(matrixscore, aes(x=as.factor(TvNT), y=Score, fill = as.factor(TvNT)))+
  theme_bw()+
  geom_violin(trim=FALSE) +
  scale_y_continuous(name = "Matrix Risk Score") +
  scale_x_discrete(name = "", limits = c("0", "1"), labels = c("Non-Tumour", "Tumour"))+
  scale_fill_discrete(name = "", limits = c("0", "1"), labels = c("Non-Tumour", "Tumour"))+
  ggtitle("LUSC Matrix Score (a=0.5)")+
  geom_boxplot(width = 0.1)+
  labs(fill = "Tumour and Non-Tumour")+
  theme(axis.text=element_text(size = 10, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 10, face="bold", hjust= 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor= element_blank(),
        legend.position = "none")
LUSC_lasso3_matrixscore

# Supp Fig 2C Risk score vs age -------------------------------------------------------
LUSC_riskcorrage<-cor.test(matrixscore_clin$Score, matrixscore_clin$Age, method = "spearman")
#rho = -0.2287163; p=0.0007275

#scatter plot of these values
LUSC_riskscore_age_scatter<-ggplot(matrixscore_clin, aes(x=Age, y=Score))+
  geom_point()+
  theme_bw()+
  scale_y_continuous(name = "Matrix Risk Score") +
  scale_x_continuous(name = "Age at Diagnosis")+
  ggtitle("LUSC Risk Score Correlation with Age at Diagnosis; p=0.0007275, r=-0.2287")+
   geom_smooth(method = lm)+ 
  theme(axis.text=element_text(size = 10, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 12, face="bold", hjust= 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor= element_blank())
LUSC_riskscore_age_scatter

# Fig3A Consensus Cluster LUSC -------------------------
alpha<-0.05
LUSCmatrisome_TNT_contrastfitTablesig<-LUSCmatrisome_TNT_contrastfitTable[LUSCmatrisome_TNT_contrastfitTable$adj.P.Val<alpha,]

#significantly DE core matrisome genes only
LUSCmatrisome_TNT_contrastfitTablesigcore<-LUSCmatrisome_TNT_contrastfitTablesig[which(rownames(LUSCmatrisome_TNT_contrastfitTablesig) %in% as.character(core_genenames)),]

LUSC_rnaseq_logcpmTNTsigcore<-LUSC_rnaseq_logcpmTNT[rownames(LUSC_rnaseq_logcpmTNT) %in% rownames(LUSCmatrisome_TNT_contrastfitTablesigcore),]

LUSC_rnaseq_logcpmTNTsigcorez<-apply(LUSC_rnaseq_logcpmTNTsigcore,1, scale) #note that this has already transformed the matrix

LUSC_accs<-paste0(substr(colnames(LUSC_rnaseq_logcpmTNT),1,12),"_", LUSC_rnaseq_logcpmTNT[1,])

rownames(LUSC_rnaseq_logcpmTNTsigcorez)<-LUSC_accs

#Consensus clustering on the z scaled data
library(M3C)
library(ComplexHeatmap)

#distance between a vector and each row of a matrix
calc_vec2mat_dist = function(x, ref_mat) {
  # compute row-wise vec2vec distance 
  apply(ref_mat, 1, function(r) sum((r - x)^2))
}

#calculate the euclidean distance for each row of the input matrix
calc_mat2mat_dist = function(input_mat, ref_mat) {
  
  dist_mat = apply(input_mat, 1, function(r) calc_vec2mat_dist(r, ref_mat))
  
  # transpose to have each row for each input datapoint
  # each column for each centroids
  cbind(t(dist_mat), max.col(-t(dist_mat)))
}

#pearson correlation between a vector and row of a matrix
calc_vec2mat_pearson = function(x, ref_mat) {
  # compute row-wise vec2vec distance 
  pearson<-apply(ref_mat, 1, function(r) cor.test(r,x, method = "pearson")$estimate)
  #pearson$estimate
}
#calculate the pearson correlation for each row of the input matrix
calc_mat2mat_pearson = function(input_mat, ref_mat) {
  
  dist_mat = apply(input_mat, 1, function(r) calc_vec2mat_pearson(r, ref_mat))
  
  # transpose to have each row for each input datapoint
  # each column for each centroids
  cbind(t(dist_mat), max.col(t(dist_mat)))
}

LUSC_TvNT<-data.frame(TvNT = as.factor(gsub(".*_", "", rownames(LUSC_rnaseq_logcpmTNTsigcorez))),
                      ID = rownames(LUSC_rnaseq_logcpmTNTsigcorez))
#kmeans clustering algorithm
LUSC_TvNT_sigmatrisome_cckm <- M3C(t(LUSC_rnaseq_logcpmTNTsigcorez), des = LUSC_TvNT, doanalysis=TRUE,analysistype='chi',
                                   printres = TRUE, variable = "TvNT", clusteralg = "km")
LUSC_TvNT_sigmatrisome_cckm$scores 
LUSC_TvNT_sigmatrisome_cckm$clinicalres 
#Indicates 3 clusters gives best cluster separation with least over-fitting
LUSC_TvNT_sigmatrisome_ccdata <- LUSC_TvNT_sigmatrisome_cckm$realdataresults[[3]]$ordered_data # this is the data
LUSC_TvNT_sigmatrisome_ccannot <- LUSC_TvNT_sigmatrisome_cckm$realdataresults[[3]]$ordered_annotation # this is the annotation
LUSC_TvNT_sigmatrisome_ccmatrix <- LUSC_TvNT_sigmatrisome_cckm$realdataresults[[3]]$consensus_matrix # this is the consensus matrix
LUSC_TvNT_sigmatrisome_cck3pca <- pca(LUSC_TvNT_sigmatrisome_cckm, K=3)
LUSC_TvNT_sigmatrisome_cck3tsne <- tsne(LUSC_TvNT_sigmatrisome_cckm, K=3, perplex=15)

LUSC_TvNT_sigmatrisome_ccannotcompiled<-LUSC_TvNT_sigmatrisome_cckm$realdataresults[[3]]$ordered_annotation
LUSC_TvNT_sigmatrisome_ccannotcompiled<-LUSC_TvNT_sigmatrisome_ccannotcompiled[match(rownames(LUSC_TvNT_sigmatrisome_cc$realdataresults[[4]]$ordered_annotation), rownames(LUSC_TvNT_sigmatrisome_ccannotcompiled)),]

#Add clinical predictive factors into this matrix along with stage
LUSC_TvNT_sigmatrisome_ccannot2<-merge(LUSC_TvNT_sigmatrisome_ccannotcompiled, LUSC_OS, by = "row.names", all.x = TRUE)
rownames(LUSC_TvNT_sigmatrisome_ccannot2)<-LUSC_TvNT_sigmatrisome_ccannot2$Row.names
colnames(LUSC_TvNT_sigmatrisome_ccannot2)[1]<-"Km_K3"
LUSC_TvNT_sigmatrisome_ccannot2<-LUSC_TvNT_sigmatrisome_ccannot2[match(rownames(LUSC_TvNT_sigmatrisome_ccannot),rownames(LUSC_TvNT_sigmatrisome_ccannot2)),] 
LUSC_TvNT_sigmatrisome_ccannot2<-LUSC_TvNT_sigmatrisome_ccannot2[,c(-1)]
LUSC_TvNT_sigmatrisome_ccannot3<-LUSC_TvNT_sigmatrisome_ccannot2[,c(2,1,4)] #for tissue, cluster and Stage 

colnames(LUSC_TvNT_sigmatrisome_ccannot3)<-c("Tissue", "Cluster", "Stage")
levels(LUSC_TvNT_sigmatrisome_ccannot3$Tissue)<-c("NT", "T")
levels(LUSC_TvNT_sigmatrisome_ccannot3$Cluster)<-c("K1", "K2", "K3")

colours <- list("Tissue"=c("NT"="darkgrey","T"="black"),
                "Cluster"=c("K1"="gold","K2"="royalblue", "K3"="green", "K4"="red", "K5" = "purple"),
                "Stage" = c("I" = "lightskyblue", "IA" ="lightblue", "IB" = "navyblue","II" = "orange", "IIA" = "orangered", "IIB" = "red2","III" ="plum1", "IIIA" = "purple", "IIIB" = "purple4", "IV" = "black"))

LUSC_colAnn <- HeatmapAnnotation(df=LUSC_TvNT_sigmatrisome_ccannot3, which="col",
                                 col=colours,
                                 annotation_width=unit(c(1, 4), "cm"), gap=unit(1, "mm"), na_col = "white")

LUSC_colAnnclust <- HeatmapAnnotation(df=LUSC_TvNT_sigmatrisome_ccannot2[,c(2,6,3)], which="col",
                                      col=colours2,
                                      annotation_width=unit(c(1, 4), "cm"), gap=unit(1, "mm"), na_col = "white")

library(ComplexHeatmap)
col_fun <- circlize::colorRamp2(c(-4, 0, 4), c("blue", "white", "red")) #all nonsig set to zero so will be white

#Supp Figure 4C consensus matrix for Kmeans with LUSC
consensus_col_fun <- circlize::colorRamp2(c(0,0.5, 1), c("blue", "white", "red")) 
LUSC_TvNT_sigmatrisome_kmconsensus<-Heatmap(LUSC_TvNT_sigmatrisome_kmccmatrix,
                                            col = consensus_col_fun,
                                            name = "Consensus",
                                            border = TRUE,#black border
                                            cluster_rows = FALSE,
                                            cluster_columns = FALSE, #this does not work
                                            column_names_gp = gpar(fontsize = 4),
                                            row_names_gp = gpar(fontsize = 4),
                                            show_row_names = FALSE,
                                            show_column_names = FALSE)

#Figure 3A plot the heatmap
#filter out the low variance genes from this list
LUSC_rnaseq_logcpmTNTsigcore_var<-apply(LUSC_rnaseq_logcpmTNTsigcore, 1, var)
LUSC_rnaseq_logcpmTNTsigcore_mean<-apply(LUSC_rnaseq_logcpmTNTsigcore, 1, mean)
hist(LUSC_rnaseq_logcpmTNTsigcore_var, breaks = 50)
hist(LUSC_rnaseq_logcpmTNTsigcore_mean, breaks = 50)
#filter out genes with low variance <1.8 and mean expression<-0.2 for appearances
var_out<-which(LUSC_rnaseq_logcpmTNTsigcore_var<1.8)
mean_out<-which(LUSC_rnaseq_logcpmTNTsigcore_mean<(0.2))
filt_out<-unname(c(var_out, mean_out))

LUSC_rnaseq_logcpmTNTsigcore_filt<-LUSC_rnaseq_logcpmTNTsigcore[-filt_out,]
colnames(LUSC_rnaseq_logcpmTNTsigcore_filt)<-LUSC_accs

LUSC_TvNT_sigmatrisome_ccdata_filtkmorder<-LUSC_rnaseq_logcpmTNTsigcore_filt[,match(rownames(LUSC_TvNT_sigmatrisome_cckm$realdataresults[[3]]$ordered_annotation), colnames(LUSC_rnaseq_logcpmTNTsigcore_filt))]
LUSC_TvNT_sigmatrisome_ccdata_filtkmorderz<-t(apply(LUSC_TvNT_sigmatrisome_ccdata_filtkmorder, 1, scale))
colnames(LUSC_TvNT_sigmatrisome_ccdata_filtkmorderz)<-colnames(LUSC_TvNT_sigmatrisome_ccdata_filtkmorder)

LUSC_TvNT_sigmatrisome_ccmapKmfilt<-Heatmap(LUSC_TvNT_sigmatrisome_ccdata_filtkmorderz,
                                            col = col_fun,
                                            name = "z",
                                            border = TRUE,#black border
                                            cluster_rows = TRUE,
                                            cluster_columns = FALSE, #this does not work
                                            column_names_gp = gpar(fontsize = 0),
                                            row_names_gp = gpar(fontsize = 4),
                                            show_row_names = FALSE,
                                            top_annotation=LUSC_colAnn)

# GSEA Output  - Formatting for cluster analysis -----------------------------------------------------

LUSC_rnaseq_logcpmTNT_gsea<-LUSC_rnaseq_logcpmTNT[19:nrow(LUSC_rnaseq_logcpmTNT),]

LUSC_rnaseq_logcpmTNT_gsea<-data.frame(NAME = rownames(LUSC_rnaseq_logcpmTNT_gsea),
                                       DESCRIPTION = rep(NA, length = nrow(LUSC_rnaseq_logcpmTNT_gsea)),
                                       LUSC_rnaseq_logcpmTNT_gsea)
LUSC_rnaseq_logcpmTNT_gseaTonly<-LUSC_rnaseq_logcpmTNT_gsea[, c(1,2, grep("_1",colnames(LUSC_rnaseq_logcpmTNT_gsea)))]

#write out design matrix
LUSC_TvNT_sigmatrisome_ccannot2_gseaord<-LUSC_TvNT_sigmatrisome_ccannot2[order(match(rownames(LUSC_TvNT_sigmatrisome_ccannot2), colnames(LUSC_rnaseq_logcpmTNT))),]
#filter out nontumour so you can just compare K1 and K3
LUSC_TvNT_sigmatrisome_ccannot2_gseaordTonly<-LUSC_TvNT_sigmatrisome_ccannot2_gseaord[grep("_1",rownames(LUSC_TvNT_sigmatrisome_ccannot2_gseaord)),]

# Association of clusters with clinical features - Stage ------------------
LUSC_conting_clusterstage<- table(LUSC_TvNT_sigmatrisome_ccannot3$Cluster,LUSC_TvNT_sigmatrisome_ccannot3$Stage)

#removeK1 from the contingency table as is it only for NT tissue
LUSC_conting_clusterstage<-LUSC_conting_clusterstage[-1,]
LUSC_conting_clusterstagef<-fisher.test(LUSC_conting_clusterstage)
LUSC_conting_clusterstagechi<-chisq.test(LUSC_conting_clusterstage)
LUSC_conting_clusterstagechi$expected
LUSC_conting_clusterstagechi$observed

# Fig 3D Supp Fig 4EF Survival for core matrisome clusters ------------------------------------
#merge the survival data with the cluster information 
LUSC_TvNT_sigmatrisome_ccannotcompiledT<-LUSC_TvNT_sigmatrisome_ccannot2[grep("_1",rownames(LUSC_TvNT_sigmatrisome_ccannot2)),]
rownames(LUSC_TvNT_sigmatrisome_ccannotcompiledT)<-substring(rownames(LUSC_TvNT_sigmatrisome_ccannotcompiledT),1,12)

LUSC_OS_cc<-LUSC_TvNT_sigmatrisome_ccannotcompiledT
colnames(LUSC_OS_cc)[1]<-"Km_K3"
LUSC_OS_cc<-droplevels(LUSC_OS_cc)

library("survival")
library("survminer")

#test associations with clinical factors
LUSC_Stage_fit<-survfit(Surv(Timetoevent_5yrs, vital_status5yrs_numeric)~as.factor(LUSC_OS_cc$Stage),
                        data = LUSC_OS_cc)
LUSC_CoreMatriClust_KmK3_fit<-survfit(Surv(Timetoevent_5yrs, vital_status5yrs_numeric)~as.factor(LUSC_OS_cc$Km_K3),
                                      data = LUSC_OS_cc)
LUSC_CoreMatriClust_KmK3_coxph<-coxph(Surv(Timetoevent_5yrs, vital_status5yrs_numeric)~as.factor(LUSC_OS_cc$Km_K3),
                                      data = LUSC_OS_cc)
LUSC_CoreMatriClust_KmK3_coxph$coefficients
summary(LUSC_CoreMatriClust_KmK3_coxph)$conf.int
summary(pancan_molsubtype_coxph[[i]])$conf.int

#Early Stage: Combine all stage I and II samples 
LUSC_OS_cc$Stage_EarlyLate<-as.factor(LUSC_OS_cc$Stage)
levels(LUSC_OS_cc$Stage_EarlyLate)<-c("Early", #IA
                                      "Early", #IB
                                      "Early", #II
                                      "Early", #IIA
                                      "Early", #IIB
                                      "Late", #IIIA
                                      "Late", #IIIB
                                      "Late") #IV
summary(LUSC_OS_cc$Stage_EarlyLate)

#Perform cox regression on subsetted survival matrix
LUSC_CoreMatriClust_Early_KmK3_coxph<-coxph(Surv(Timetoevent_5yrs, vital_status5yrs_numeric)~Km_K3,
                                            data = subset(LUSC_OS_cc, Stage_EarlyLate=="Early"))
LUSC_CoreMatriClust_Late_KmK3_coxph<-coxph(Surv(Timetoevent_5yrs, vital_status5yrs_numeric)~Km_K3,
                                           data = subset(LUSC_OS_cc, Stage_EarlyLate=="Late"))

#incrementally by each stage
#IA
LUSC_CoreMatriClust_StageIA_KmK3_coxph<-coxph(Surv(Timetoevent_5yrs, vital_status5yrs_numeric)~Km_K3,
                                              data = subset(LUSC_OS_cc, Stage=="IA"))
#IB
LUSC_CoreMatriClust_StageIB_KmK3_coxph<-coxph(Surv(Timetoevent_5yrs, vital_status5yrs_numeric)~Km_K3,
                                              data = subset(LUSC_OS_cc, Stage=="IB"))
#IIA
LUSC_CoreMatriClust_StageIIA_KmK3_coxph<-coxph(Surv(Timetoevent_5yrs, vital_status5yrs_numeric)~Km_K3,
                                               data = subset(LUSC_OS_cc, Stage=="IIA"))
#IIB
LUSC_CoreMatriClust_StageIIB_KmK3_coxph<-coxph(Surv(Timetoevent_5yrs, vital_status5yrs_numeric)~Km_K3,
                                               data = subset(LUSC_OS_cc, Stage=="IIB"))
#IIIA
LUSC_CoreMatriClust_StageIIIA_KmK3_coxph<-coxph(Surv(Timetoevent_5yrs, vital_status5yrs_numeric)~Km_K3,
                                                data = subset(LUSC_OS_cc, Stage=="IIIA"))
#IIIB #only one event so not enough to run the analysis
LUSC_CoreMatriClust_StageIIIB_KmK3_coxph<-coxph(Surv(Timetoevent_5yrs, vital_status5yrs_numeric)~Km_K3,
                                                data = subset(LUSC_OS_cc, Stage=="IIIB"))

LUSC_matcluster_bystage<-list(LUSC_CoreMatriClust_Early_KmK3_coxph,
                              LUSC_CoreMatriClust_Late_KmK3_coxph,
                              LUSC_CoreMatriClust_StageIA_KmK3_coxph,
                              LUSC_CoreMatriClust_StageIB_KmK3_coxph,
                              LUSC_CoreMatriClust_StageIIA_KmK3_coxph,
                              LUSC_CoreMatriClust_StageIIB_KmK3_coxph,
                              LUSC_CoreMatriClust_StageIIIA_KmK3_coxph)

LUSC_matcluster_bystage_df<-matrix(NA, ncol = 5, nrow = length(LUSC_matcluster_bystage))
for (i in 1: length(LUSC_matcluster_bystage)){
  LUSC_matcluster_bystage_df[i,1] <-summary(LUSC_matcluster_bystage[[i]])$wald["pvalue"] #pvalue
  LUSC_matcluster_bystage_df[i,2] <-summary(LUSC_matcluster_bystage[[i]])$coef[1,1]#coefficient beta
  LUSC_matcluster_bystage_df[i,3] <-summary(LUSC_matcluster_bystage[[i]])$coef[1,2]# HR
  LUSC_matcluster_bystage_df[i,4] <- summary(LUSC_matcluster_bystage[[i]])$conf.int[,"lower .95"][2] #lower confidence interval
  LUSC_matcluster_bystage_df[i,5] <- summary(LUSC_matcluster_bystage[[i]])$conf.int[,"upper .95"][2] #upper confidence interval
}

#create a df of these parameters
rownames(LUSC_matcluster_bystage_df)<- c("Early", "Late", "IA", "IB", "IIA", "IIB", "IIIA")
colnames(LUSC_matcluster_bystage_df)<- c("P.value", "Beta", "HR", "LowerCI", "UpperCI")     

#Supp Fig 4E plot the survival analysis for the early stage subset of patients
LUSC_CoreMatriClust_Early_KmK3_fit<-survfit(Surv(Timetoevent_5yrs, vital_status5yrs_numeric)~Km_K3,
                                            data = subset(LUSC_OS_cc, Stage_EarlyLate=="Early"))
LUSC_CoreMatClusters_EarlyStage_KM<-ggsurvplot(LUSC_CoreMatriClust_Early_KmK3_fit,
                                               pval = TRUE,
                                               risk.table = FALSE,
                                               title = "Early Stages (I + II)",
                                               xlab = "Days Since Diagnosis",
                                               ylab = "Survival Probability",
                                               palette = c("royalblue", "forestgreen"),
                                               risk.table.col = "strata",
                                               legend = c(0.9,0.9),
                                               legend.title = "Matreotype",
                                               legend.labs = c("ECM-Low", "ECM-High"),
                                               font.main = 20,
                                               font.legend = 20,
                                               font.x = 20,
                                               font.y = 20,
                                               font.p = 20,
                                               font.tickslab = c(20))

#plot the survival analysis for the late stage subset of patients
LUSC_CoreMatriClust_Late_KmK3_fit<-survfit(Surv(Timetoevent_5yrs, vital_status5yrs_numeric)~Km_K3,
                                           data = subset(LUSC_OS_cc, Stage_EarlyLate=="Late"))
LUSC_CoreMatClusters_LateStage_KM<-ggsurvplot(LUSC_CoreMatriClust_Late_KmK3_fit,
                                              pval = TRUE,
                                              risk.table = FALSE,
                                              title = "Late Stages (III + IV)",
                                              xlab = "Days Since Diagnosis",
                                              ylab = "Survival Probability",
                                              palette = c("royalblue", "forestgreen"),
                                              risk.table.col = "strata",
                                              legend = c(0.9,0.9),
                                              legend.title = "Matreotype",
                                              legend.labs = c("ECM-Low", "ECM-High"),
                                              font.main = 20,
                                              font.legend = 20,
                                              font.x = 20,
                                              font.y = 20,
                                              font.p = 20,
                                              font.tickslab = c(20))

# Table 2 Overlap with Canonical Molecular Subtypes -------------------------------
#read in canonical molecular subtype centroids from Wilkerson et al.
LUSC_molsubtypes<-read.csv("./LUSC.predictor.centroids.csv", header = TRUE, row.names = 1)

#LUSC
#Assign LUSC samples to canonical molecular subtypes by pearson correlation
LUSC_molsubtypes_overlap<-LUSC_molsubtypes[which(rownames(LUSC_molsubtypes) %in% colnames(LUSC_rnaseq_logcpmTNTtz)),]

LUSC_rnaseq_logcpmTNTtz_molsub<-LUSC_rnaseq_logcpmTNTtz[, which(colnames(LUSC_rnaseq_logcpmTNTtz) %in% rownames(LUSC_molsubtypes))]
LUSC_rnaseq_logcpmTNTtz_molsub_ord<-LUSC_rnaseq_logcpmTNTtz_molsub[, order(match(colnames(LUSC_rnaseq_logcpmTNTtz_molsub), rownames(LUSC_molsubtypes_overlap)))]

#calculate pearson correlation coefficient between sample and molecular subtype
LUSC_canonicalmolsubtypes<-calc_mat2mat_pearson(LUSC_rnaseq_logcpmTNTtz_molsub_ord, t(LUSC_molsubtypes_overlap))
LUSC_canonicalmolsubtypes<-as.data.frame(LUSC_canonicalmolsubtypes)
LUSC_canonicalmolsubtypes[,ncol(LUSC_canonicalmolsubtypes)]<-as.factor(LUSC_canonicalmolsubtypes[,ncol(LUSC_canonicalmolsubtypes)])
levels(LUSC_canonicalmolsubtypes[,ncol(LUSC_canonicalmolsubtypes)])<-colnames(LUSC_canonicalmolsubtypes)[1:4]
colnames(LUSC_canonicalmolsubtypes)[5]<-"LUSC_MolSubtype"
rownames(LUSC_canonicalmolsubtypes)<-paste0(LUSC_rnaseq_logcpmTNTt$Acc, "_", LUSC_rnaseq_logcpmTNTt$TvNT)

LUSC_TvNT_sigmatrisome_ccannot3<-merge(LUSC_TvNT_sigmatrisome_ccannot2, LUSC_canonicalmolsubtypes, by = "row.names")

LUSC_TvNT_sigmatrisome_ccannot3$primvsnot<-ifelse(LUSC_TvNT_sigmatrisome_ccannot3$LUSC_MolSubtype=="primitive",1,0)
LUSC_TvNT_sigmatrisome_ccannot3$classvsnot<-ifelse(LUSC_TvNT_sigmatrisome_ccannot3$LUSC_MolSubtype=="classical",1,0)
LUSC_TvNT_sigmatrisome_ccannot3$secvsnot<-ifelse(LUSC_TvNT_sigmatrisome_ccannot3$LUSC_MolSubtype=="secretory",1,0)
LUSC_TvNT_sigmatrisome_ccannot3$basalvsnot<-ifelse(LUSC_TvNT_sigmatrisome_ccannot3$LUSC_MolSubtype=="basal",1,0)

#convert cluster assignment to binary for each column
LUSC_TvNT_sigmatrisome_ccannot3$K1vsnot<-ifelse(LUSC_TvNT_sigmatrisome_ccannot3$Km_K3=="1",1,0)
LUSC_TvNT_sigmatrisome_ccannot3$K2vsnot<-ifelse(LUSC_TvNT_sigmatrisome_ccannot3$Km_K3=="2",1,0)
LUSC_TvNT_sigmatrisome_ccannot3$K3vsnot<-ifelse(LUSC_TvNT_sigmatrisome_ccannot3$Km_K3=="3",1,0)

LUSC_TvNT_sigmatrisome_ccannot3T<-droplevels(LUSC_TvNT_sigmatrisome_ccannot3[grep("_1", LUSC_TvNT_sigmatrisome_ccannot3$Row.names),])

#Fishers exact test across whole cohort
LUSC_TvNT_sigmatrisome_ccannot3T_molsubtypeK<-table(LUSC_TvNT_sigmatrisome_ccannot3T$Km_K3, LUSC_TvNT_sigmatrisome_ccannot3T$LUSC_MolSubtype)
LUSC_TvNT_sigmatrisome_ccannot3T_molsubtypeKprop<-prop.table(LUSC_TvNT_sigmatrisome_ccannot3T_molsubtypeK,2)
LUSC_TvNT_sigmatrisome_ccannot3T_molsubtypeKf<-fisher.test(LUSC_TvNT_sigmatrisome_ccannot3T_molsubtypeK, workspace = 2e8)

#subset by molecular subtypes
#primitive
LUSC_TvNT_sigmatrisome_ccannot3T_primK1<-table(LUSC_TvNT_sigmatrisome_ccannot3T$K1vsnot, LUSC_TvNT_sigmatrisome_ccannot3T$primvsnot)
LUSC_TvNT_sigmatrisome_ccannot3T_primK3<-table(LUSC_TvNT_sigmatrisome_ccannot3T$K3vsnot, LUSC_TvNT_sigmatrisome_ccannot3T$primvsnot)

LUSC_TvNT_sigmatrisome_ccannot3T_primK1f<-fisher.test(LUSC_TvNT_sigmatrisome_ccannot3T_primK1) 
LUSC_TvNT_sigmatrisome_ccannot3T_primK3f<-fisher.test(LUSC_TvNT_sigmatrisome_ccannot3T_primK3) 

#classical
LUSC_TvNT_sigmatrisome_ccannot3T_classK1<-table(LUSC_TvNT_sigmatrisome_ccannot3T$K1vsnot, LUSC_TvNT_sigmatrisome_ccannot3T$classvsnot)
LUSC_TvNT_sigmatrisome_ccannot3T_classK3<-table(LUSC_TvNT_sigmatrisome_ccannot3T$K3vsnot, LUSC_TvNT_sigmatrisome_ccannot3T$classvsnot)

LUSC_TvNT_sigmatrisome_ccannot3T_classK1f<-fisher.test(LUSC_TvNT_sigmatrisome_ccannot3T_classK1) 
LUSC_TvNT_sigmatrisome_ccannot3T_classK3f<-fisher.test(LUSC_TvNT_sigmatrisome_ccannot3T_classK3) 

#secretory
LUSC_TvNT_sigmatrisome_ccannot3T_secK1<-table(LUSC_TvNT_sigmatrisome_ccannot3T$K1vsnot, LUSC_TvNT_sigmatrisome_ccannot3T$secvsnot)
LUSC_TvNT_sigmatrisome_ccannot3T_secK3<-table(LUSC_TvNT_sigmatrisome_ccannot3T$K3vsnot, LUSC_TvNT_sigmatrisome_ccannot3T$secvsnot)

LUSC_TvNT_sigmatrisome_ccannot3T_secK1f<-fisher.test(LUSC_TvNT_sigmatrisome_ccannot3T_secK1) 
LUSC_TvNT_sigmatrisome_ccannot3T_secK3f<-fisher.test(LUSC_TvNT_sigmatrisome_ccannot3T_secK3) 

#basal
#secretory
LUSC_TvNT_sigmatrisome_ccannot3T_basK1<-table(LUSC_TvNT_sigmatrisome_ccannot3T$K1vsnot, LUSC_TvNT_sigmatrisome_ccannot3T$basalvsnot)
LUSC_TvNT_sigmatrisome_ccannot3T_basK3<-table(LUSC_TvNT_sigmatrisome_ccannot3T$K3vsnot, LUSC_TvNT_sigmatrisome_ccannot3T$basalvsnot)

LUSC_TvNT_sigmatrisome_ccannot3T_basK1f<-fisher.test(LUSC_TvNT_sigmatrisome_ccannot3T_basK1) 
LUSC_TvNT_sigmatrisome_ccannot3T_basK3f<-fisher.test(LUSC_TvNT_sigmatrisome_ccannot3T_basK3) 

# Supp Fig 4D LUAD LUSC Cluster DEG ---------------------------------------------------
design_LUSC_T_KmClust<-model.matrix(~0+Km_K3, data = LUSC_OS_cc)
colnames(design_LUSC_T_KmClust)<-c("K1", "K3")

#reorder expression matrix 
LUSC_rnaseq_logcpmTt_ord<-LUSC_rnaseq_logcpmTt
rownames(LUSC_rnaseq_logcpmTt_ord)<-substr(rownames(LUSC_rnaseq_logcpmTt_ord),1,12)
LUSC_rnaseq_logcpmTt_ord<-LUSC_rnaseq_logcpmTt_ord[match(rownames(LUSC_OS_cc), rownames(LUSC_rnaseq_logcpmTt_ord)),]

fit_LUSC_KmClust <- lmFit(t(LUSC_rnaseq_logcpmTt_ord[,2:ncol(LUSC_rnaseq_logcpmTt_ord)]), design_LUSC_T_KmClust)
LUSC_KmClust_contrast<-makeContrasts(coef1 =K3-K1, levels = design_LUSC_T_KmClust)
LUSC_KmClust_contrastfit<-contrasts.fit(fit_LUSC_KmClust, LUSC_KmClust_contrast)
LUSC_KmClust_contrastfit<-eBayes(LUSC_KmClust_contrastfit)
LUSC_KmClust_contrastfitTable<-topTable(LUSC_KmClust_contrastfit, n= 'Inf')

#matrisome only
LUSC_KmClust_contrastfitTable_mat<-LUSC_KmClust_contrastfitTable[which(LUSC_KmClust_contrastfitTable$ID %in% matrisome_master_min$Gene.Symbol),]

#core matrisome only
LUSC_KmClust_contrastfitTable_coremat<-LUSC_KmClust_contrastfitTable[which(LUSC_KmClust_contrastfitTable$ID %in% core_genenames),]

#subset by up and downregulated
LUSC_KmClust_contrastfitTable_up<-LUSC_KmClust_contrastfitTable[which(LUSC_KmClust_contrastfitTable$adj.P.Val<0.05 & LUSC_KmClust_contrastfitTable$logFC>0),]
LUSC_KmClust_contrastfitTable_dn<-LUSC_KmClust_contrastfitTable[which(LUSC_KmClust_contrastfitTable$adj.P.Val<0.05 & LUSC_KmClust_contrastfitTable$logFC<0),]

#all matrisome up and downregulated
LUSC_KmClust_contrastfitTable_mat_up<-LUSC_KmClust_contrastfitTable_mat[which(LUSC_KmClust_contrastfitTable_mat$adj.P.Val<0.05 & LUSC_KmClust_contrastfitTable_mat$logFC>0),]
LUSC_KmClust_contrastfitTable_mat_dn<-LUSC_KmClust_contrastfitTable_mat[which(LUSC_KmClust_contrastfitTable_mat$adj.P.Val<0.05 & LUSC_KmClust_contrastfitTable_mat$logFC<0),]

#core matrisome up and downregulated
LUSC_KmClust_contrastfitTable_coremat_up<-LUSC_KmClust_contrastfitTable_coremat[which(LUSC_KmClust_contrastfitTable_coremat$adj.P.Val<0.05 & LUSC_KmClust_contrastfitTable_coremat$logFC>0),]
LUSC_KmClust_contrastfitTable_coremat_dn<-LUSC_KmClust_contrastfitTable_coremat[which(LUSC_KmClust_contrastfitTable_coremat$adj.P.Val<0.05 & LUSC_KmClust_contrastfitTable_coremat$logFC<0),]
LUSC_TvNT_sigmatrisome_ccannot3_kmord<-LUSC_TvNT_sigmatrisome_ccannot3[match(LUSC_OS_cc$Acc,LUSC_TvNT_sigmatrisome_ccannot3$Acc),]

LUSCcolours2 <- list("TvNT"=c("NT"="darkgrey","T"="black"),
                     "Km_K3"=c("1"="gold","2"="royalblue", "3"="forestgreen"),
                     "Stage" = c( "IA" ="lightblue", "IB" = "navyblue","II" = "orange", "IIA" = "orangered", "IIB" = "red2", "IIIA" = "purple", "IIIB" = "purple4", "IV" = "black"),
                     "LUSC_MolSubtype" = c("primitive" = "gold", "classical" = "red2", "secretory" = "purple4", "basal" = "royalblue"))
# 

LUSC_colAnnclust_kmmolsubtype <- HeatmapAnnotation(df=LUSC_TvNT_sigmatrisome_ccannot3_kmord[,c(2,5,29)], which="col",
                                                   col=LUSCcolours2,
                                                   annotation_width=unit(c(1, 4), "cm"), gap=unit(1, "mm"), na_col = "white")

#narrow the colour range
col_fun <- circlize::colorRamp2(c(-2.5, 0, 2.5), c("blue", "white", "red")) #all nonsig set to zero so will be white

genenumber<-10 #choose number of genes to show as up and downregulated

#Matrisomal genes only
genesup<-LUSC_KmClust_contrastfitTable_mat_up$ID[1:genenumber]
genesdn<-LUSC_KmClust_contrastfitTable_mat_dn$ID[1:genenumber]
LUSC_clustDEGs_matgenes<-LUSC_rnaseq_logcpmTt[,match(c(genesup, genesdn), colnames(LUSC_rnaseq_logcpmTt))]

LUSC_clustDEGs_matgenes<-apply(LUSC_clustDEGs_matgenes,2, scale)
rownames(LUSC_clustDEGs_matgenes)<-paste0(substr(rownames(LUSC_rnaseq_logcpmTt),1,12), "_1")

#order Accs in order of cluster matrix
LUSC_clustDEGs_matgenes<-LUSC_clustDEGs_matgenes[match(LUSC_TvNT_sigmatrisome_ccannot3_kmord$Acc, rownames(LUSC_clustDEGs_matgenes)),]

#Supp Fig 4D
LUSC_KmK3_clustmatgenestop10_hmap<-Heatmap(t(LUSC_clustDEGs_matgenes),
                                           col = col_fun,
                                           name = "z",
                                           border = TRUE,#black border
                                           cluster_rows = FALSE,
                                           cluster_columns = FALSE, #this does not work
                                           column_names_gp = gpar(fontsize = 4),
                                           row_names_gp = gpar(fontsize = 8),
                                           show_row_names = TRUE,
                                           top_annotation=LUSC_colAnnclust_kmmolsubtype,
                                           show_column_names = FALSE)

# Supp Table 4 Demographic info --------------------------------------------------------
LUSC_OS_cc
#Age
mean(LUSC_OS_cc$Age, na.rm =TRUE)
sd(LUSC_OS_cc$Age, na.rm = TRUE)

LUSC_OS_ccK2<-subset(LUSC_OS_cc, Km_K3==2)
LUSC_OS_ccK3<-subset(LUSC_OS_cc, Km_K3==3)
# LUSC_OS_ccnotK2<-subset(LUSC_OS_cc, Km_K3!=2)

mean(LUSC_OS_ccK2$Age, na.rm =TRUE)
sd(LUSC_OS_ccK2$Age, na.rm = TRUE)

mean(LUSC_OS_ccK3$Age, na.rm =TRUE)
sd(LUSC_OS_ccK3$Age, na.rm = TRUE)

wilcox.test(Age~ Km_K3, LUSC_OS_cc)

#Gender
summary(LUSC_OS_cc$Gender, na.rm =TRUE)
summary(LUSC_OS_ccK2$Gender, na.rm =TRUE)
summary(LUSC_OS_ccK3$Gender, na.rm =TRUE)

gendertable<-table(LUSC_OS_cc$Gender,as.factor(LUSC_OS_cc$Km_K3))
fisher.test(gendertable, LUSC_OS_cc)

#Race
summary(LUSC_OS_cc$Race, na.rm =TRUE)
summary(LUSC_OS_ccK2$Race, na.rm =TRUE)
summary(LUSC_OS_ccK3$Race, na.rm =TRUE)

racetable<-table(LUSC_OS_cc$Race,as.factor(LUSC_OS_cc$Km_K3))
fisher.test(racetable, LUSC_OS_cc)

#Stage
summary(as.factor(LUSC_OS_cc$Stage), na.rm =TRUE)
summary(as.factor(LUSC_OS_ccK2$Stage), na.rm =TRUE)
summary(as.factor(LUSC_OS_ccK3$Stage), na.rm =TRUE)

stagetable<-table(LUSC_OS_cc$Stage,as.factor(LUSC_OS_cc$Km_K3))
fisher.test(stagetable, LUSC_OS_cc)

#Smoking status
summary(as.factor(LUSC_OS_cc$Smokingstatus), na.rm =TRUE)
summary(as.factor(LUSC_OS_ccK2$Smokingstatus), na.rm =TRUE)
summary(as.factor(LUSC_OS_ccK3$Smokingstatus), na.rm =TRUE)

smokingtable<-table(LUSC_OS_cc$Smokingstatus,as.factor(LUSC_OS_cc$Km_K3))
fisher.test(smokingtable, LUSC_OS_cc)

#Packyrs
mean(LUSC_OS_cc$Packyrs, na.rm =TRUE)
sd(LUSC_OS_cc$Packyrs, na.rm = TRUE)

mean(LUSC_OS_ccK2$Packyrs, na.rm =TRUE)
sd(LUSC_OS_ccK2$Packyrs, na.rm = TRUE)

mean(LUSC_OS_ccK3$Packyrs, na.rm =TRUE)
sd(LUSC_OS_ccK3$Packyrs, na.rm = TRUE)

wilcox.test(Packyrs~ Km_K3, LUSC_OS_cc)

# Fig 4A-C ESTIMATE score ----------------------------------------------------------

LUSC_estimatev2<-read.table(file = "/home/z3403160/LungTCGA/CodeTest/LUSC_EstimateRNAseqv2_IDStromal.txt", sep = "\t", header = T)

LUSC_estimatev2$ID<-gsub("-", "\\.", as.character(LUSC_estimatev2$ID))
LUSC_estimatev2$ID<-substr(LUSC_estimatev2$ID, 1, 12)

LUSC_OS_cc2<-merge(LUSC_OS_cc, LUSC_estimatev2, by.x = "row.names", by.y = "ID", keep.x = TRUE)

#determine if the scores are significantly different between clusters
LUSC_estimate_K1vsK3stromalwilcox<-wilcox.test(LUSC_OS_cc2$Stromal_score~ LUSC_OS_cc2$Km_K3, data = LUSC_OS_cc2) #p<2.2E-16
LUSC_estimate_K1vsK3immunewilcox<-wilcox.test(LUSC_OS_cc2$Immune_score~ LUSC_OS_cc2$Km_K3, data = LUSC_OS_cc2) #p=8.004E-10
LUSC_estimate_K1vsK3puritywilcox<-wilcox.test(LUSC_OS_cc2$ESTIMATE_score~ LUSC_OS_cc2$Km_K3, data = LUSC_OS_cc2) #p<2.2E-16

#stromal score
LUSC_estimate_cluststromal<-ggplot(LUSC_OS_cc2, aes_string(x=as.factor(LUSC_OS_cc2$Km_K3), y=LUSC_OS_cc2$Stromal_score, fill = as.factor(LUSC_OS_cc2$Km_K3)))+
  theme_bw()+
  geom_violin(trim=FALSE) +
  scale_y_continuous(name = "Stromal Score") +
  scale_x_discrete(name = "ECM Molecular Subtype")+
  scale_fill_manual(values = c("royalblue", "green"))+
  ggtitle("Stromal")+
  geom_boxplot(alpha = 0.7, width = 0.1)+
  labs(fill = "Cluster")+
  theme(axis.text=element_text(size = 10, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 12, face="bold", hjust= 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor= element_blank(),
        legend.position = "none")

#immune score
LUSC_estimate_clustimmune<-ggplot(LUSC_OS_cc2, aes_string(x=as.factor(LUSC_OS_cc2$Km_K3), y=LUSC_OS_cc2$Immune_score, fill = as.factor(LUSC_OS_cc2$Km_K3)))+
  theme_bw()+
  geom_violin(trim=FALSE) +
  scale_y_continuous(name = "Immune Score") +
  scale_x_discrete(name = "ECM Molecular Subtype")+
  scale_fill_manual(values = c("royalblue", "green"))+
  ggtitle("Immune")+
  geom_boxplot(alpha = 0.7, width = 0.1)+
  labs(fill = "Cluster")+
  theme(axis.text=element_text(size = 10, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 12, face="bold", hjust= 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor= element_blank(),
        legend.position = "none")

#tumour purity
LUSC_estimate_clustpurity<-ggplot(LUSC_OS_cc2, aes_string(x=as.factor(LUSC_OS_cc2$Km_K3), y=LUSC_OS_cc2$ESTIMATE_score, fill = as.factor(LUSC_OS_cc2$Km_K3)))+
  theme_bw()+
  geom_violin(trim=FALSE) +
  scale_y_continuous(name = "ESTIMATE Score") +
  scale_x_discrete(name = "ECM Molecular Subtype")+
  scale_fill_manual(values = c("royalblue", "green"))+
  ggtitle("Purity")+
  geom_boxplot(alpha = 0.7, width = 0.1)+
  labs(fill = "Cluster")+
  theme(axis.text=element_text(size = 10, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 12, face="bold", hjust= 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor= element_blank(),
        legend.position = "none")

library("cowplot")
LUSC_estimate_plots<-ggdraw() +
  draw_plot(LUSC_estimate_cluststromal, x = 0, y = 0, width = .33, height = 1)+
  draw_plot(LUSC_estimate_clustimmune, x = 0.33, y = 0, width = .33, height = 1)+
  draw_plot(LUSC_estimate_clustpurity, x = 0.66, y = 0, width = .33, height = 1)

#coxph model of scores as continuous variables
#tumor purity
LUSC_estimate_coxph<-coxph(Surv(Timetoevent_5yrs, vital_status5yrs_numeric) ~ ESTIMATE_score , data = LUSC_OS_cc2)
summary(LUSC_estimate_coxph)$coef #HR = 1.0001 [0.999-1] p=0.2
summary(LUSC_estimate_coxph)$conf.int

#estimate stromal score
LUSC_stromal_coxph<-coxph(Surv(Timetoevent_5yrs, vital_status5yrs_numeric) ~ Stromal_score , data = LUSC_OS_cc2)
summary(LUSC_stromal_coxph)$conf.int #HR = 1.000308 [1.000016-1.000599] p=0.0386

#estimate immune score
LUSC_immune_coxph<-coxph(Surv(Timetoevent_5yrs, vital_status5yrs_numeric) ~ Immune_score , data = LUSC_OS_cc2)
summary(LUSC_immune_coxph)$conf.int #HR = 1.0005 [0.999-1.000322] p=0.717

# Fig 7D Fibrosis Score by Cluster -----------------------------------------------
FibrosisGeneList<-read.csv(file = './McDonough2019_IPFgenelist.csv')

LUSC_fib_cols<-which(colnames(LUSC_rnaseq_logcpmTNTtz) %in% FibrosisGeneList$Gene)
LUSC_fibmatrix<-LUSC_rnaseq_logcpmTNTtz[,LUSC_fib_cols]
#filter genelist for only those that are in the expression matrix
LUSC_fibrosisgenelist<-FibrosisGeneList[which(FibrosisGeneList$Gene %in% colnames(LUSC_fibmatrix)),]
LUSC_fibrosisgenelist<-LUSC_fibrosisgenelist[match(LUSC_fibrosisgenelist$Gene,colnames(LUSC_fibmatrix)),]

LUSC_fibmatrix<-data.frame(t(LUSC_fibmatrix),
                           Fibrosis_Coeff = LUSC_fibrosisgenelist[,2])

LUSC_TNT_fib_score<-c()
genescore<-c()
for (i in 1:(ncol(LUSC_fibmatrix)-1)){
  for (j in 1:nrow(LUSC_fibmatrix)){
    genescore[j]<-LUSC_fibmatrix[j,i]*LUSC_fibmatrix[j,ncol(LUSC_fibmatrix)]
  }
  LUSC_TNT_fib_score[i]<-sum(genescore, rm.na= TRUE)
}
LUSC_TNT_fib_score<-data.frame(Acc = rownames(LUSC_rnaseq_logcpmTNTtz),
                               LUSC_TNT_fib_score)
rownames(LUSC_TNT_fib_score)<-LUSC_TNT_fib_score$Acc
LUSC_TNT_fib_score$Acc<-paste0(substr(LUSC_TNT_fib_score$Acc,1,12), "_",
                               ifelse(substr(LUSC_TNT_fib_score$Acc,14,14)==1,0,1)) 
LUSC_TvNT_sigmatrisome_ccannot3<-merge(LUSC_TvNT_sigmatrisome_ccannot3,
                                       LUSC_TNT_fib_score,
                                       by.x = "Row.names", by.y = "Acc", all.x = TRUE)
rownames(LUSC_TvNT_sigmatrisome_ccannot3)<-LUSC_TvNT_sigmatrisome_ccannot3$Row.names

#check association with Kmeans cluster annotations
LUSC_fib_clustKmeans_TonlyK2vK3_wilcox<-wilcox.test(LUSC_TNT_fib_score~Km_K3, data = subset(LUSC_TvNT_sigmatrisome_ccannot3, Km_K3 !=1))

LUSC_fib_clustKm_violin<-ggplot(LUSC_TvNT_sigmatrisome_ccannot3, aes(x=as.factor(Km_K3), y=LUSC_TNT_fib_score, fill = as.factor(Km_K3)))+
  theme_bw()+
  geom_violin(trim=FALSE) +
  scale_y_continuous(name = "Fibrosis Score") +
  scale_x_discrete(name = "ECM Molecular Subtype")+
  scale_fill_manual(values = c("gold", "royalblue", "forestgreen"))+
  ggtitle("LUSC Fibrosis Score (McDonough)")+
  geom_boxplot(alpha = 0.7, width = 0.1)+
  labs(fill = "Cluster")+
  xlab("ECM Molecular Subtypes")+
  theme(axis.text=element_text(size = 10, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 12, face="bold", hjust= 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor= element_blank(),
        legend.position = "none")
LUSC_fib_clustKm_violin

# Fig 4 D-I Cibersort_LambreschtsscRNAseq -------------------------------------------
LUSC_cs<-read.table(file = "./LUSC_cibersort_LambreschtsscRNAseq_DeconResults.txt", header = TRUE, sep = "\t")
#append a tumour and non-tumour column to these.
LUSC_cs_TNT<-gsub(".*_","", LUSC_cs$Mixture)

LUSC_cs2<-data.frame(LUSC_cs,
                     TNT = as.factor(LUSC_cs_TNT))

#calculate T vs NT differences in scores
LUSC_cs_TNTwilcoxres<-data.frame()
for (i in 2:(ncol(LUSC_cs2)-5)){
  LUSC_cs_TNT_wtest<-wilcox.test(LUSC_cs2[,i] ~ TNT, LUSC_cs2)
  LUSC_cs_TNTwilcoxres[i,1]<-colnames(LUSC_cs2)[i] #name of cell type
  LUSC_cs_TNTwilcoxres[i,2]<-mean(LUSC_cs2[LUSC_cs2$TNT==0,i]) #mean NT
  LUSC_cs_TNTwilcoxres[i,3]<-mean(LUSC_cs2[LUSC_cs2$TNT==1,i]) #mean T
  LUSC_cs_TNTwilcoxres[i,4]<-as.numeric(LUSC_cs_TNTwilcoxres[i,3])/as.numeric(LUSC_cs_TNTwilcoxres[i,2]) #fold change T/NT
  LUSC_cs_TNTwilcoxres[i,5]<-as.numeric(LUSC_cs_TNT_wtest$p.value) #pvalue
}
LUSC_cs_TNTwilcoxres<-LUSC_cs_TNTwilcoxres[-1,]
LUSC_cs_TNTwilcoxres<-data.frame(LUSC_cs_TNTwilcoxres, 
                                 adjPval = p.adjust(LUSC_cs_TNTwilcoxres[,5]))

colnames(LUSC_cs_TNTwilcoxres)<-c("CellType", "Mean_NT", "Mean_T", "FC_TvsNT", "Pval", "adjPval")
rownames(LUSC_cs_TNTwilcoxres)<-LUSC_cs_TNTwilcoxres[,1]

#LUSC cluster annotations
#cluster annotations contained in LUSC_TvNT_sigmatrisome_ccannot2
LUSC_TvNT_sigmatrisome_ccannot3_csord<-LUSC_TvNT_sigmatrisome_ccannot3[order(match(rownames(LUSC_TvNT_sigmatrisome_ccannot3), LUSC_cs2$Mixture)),]
#append cluster info to the cibersort matrix
LUSC_cs2_clust<-data.frame(LUSC_cs2, LUSC_TvNT_sigmatrisome_ccannot3_csord)


#run wilcox test for each cell type
LUSC_cs2_clust_K1vsK2<-subset(LUSC_cs2_clust, Km_K3!=3)
LUSC_cs2_clust_K1vsK3<-subset(LUSC_cs2_clust, Km_K3!=2)
LUSC_cs2_clust_K2vsK3<-subset(LUSC_cs2_clust, Km_K3!=1)


LUSC_cs_TNTwilcoxres_clust<-data.frame()
for (i in 2:53){
  
  #K1vsK2
  LUSC_cs_TNT_wtest_clustK1K2<-wilcox.test(LUSC_cs2_clust_K1vsK2[,i] ~ Km_K3, LUSC_cs2_clust_K1vsK2)
  
  #K1vsK3
  LUSC_cs_TNT_wtest_clustK1K3<-wilcox.test(LUSC_cs2_clust_K1vsK3[,i] ~ Km_K3, LUSC_cs2_clust_K1vsK3)
  
  #K1vsK2
  LUSC_cs_TNT_wtest_clustK2K3<-wilcox.test(LUSC_cs2_clust_K2vsK3[,i] ~ Km_K3, LUSC_cs2_clust_K2vsK3)
  
  #store the pvalues
  LUSC_cs_TNTwilcoxres_clust[i,1]<-as.numeric(LUSC_cs_TNT_wtest_clustK1K2$p.value) #pvalue
  LUSC_cs_TNTwilcoxres_clust[i,2]<-as.numeric(LUSC_cs_TNT_wtest_clustK1K3$p.value) #pvalue
  LUSC_cs_TNTwilcoxres_clust[i,3]<-as.numeric(LUSC_cs_TNT_wtest_clustK2K3$p.value) #pvalue
  
}
colnames(LUSC_cs_TNTwilcoxres_clust)<-c("K1vsK2", "K1vsK3", "K2vsK3")
rownames(LUSC_cs_TNTwilcoxres_clust)<-colnames(LUSC_cs2_clust)[1:53]

#generate violin plots of fibroblast subtypes in different clusters
LUSC_cs_celltype<-lapply(colnames(LUSC_cs2_clust)[2:53], function(x) {
  maintitle<-paste0("LUSC ", x)
  ggplot(LUSC_cs2_clust, aes_string(x=as.factor(LUSC_cs2_clust$Km_K3), y=x, fill = as.factor(LUSC_cs2_clust$Km_K3)))+
    theme_bw()+
    geom_violin(trim=FALSE) +
    scale_y_continuous(name = x) +
    scale_x_discrete(name = "ECM Molecular Subtype")+
    scale_fill_manual(values = c( "gold","royalblue", "green"))+
    ggtitle(maintitle)+
    geom_boxplot(alpha = 0.7, width = 0.1)+
    labs(fill = "Cluster")+
    theme(axis.text=element_text(size = 10, face="bold"),
          axis.title=element_text(size=10, face="bold"),
          plot.title=element_text(size = 12, face="bold", hjust= 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor= element_blank(),
          legend.position = "none")
})

# Fig 4J-K; 5A-B; 6B,H-K; 7A-B,E; Supp6, GSVA for composition ----------------------------------------------------
library(GSVA)
library(GSVAdata)

#read in the manually curated gene signatures
celltypeslist<-read.csv("./celltypes_mastergenelist.csv", header = TRUE)

#convert the mouse genes into human symbols
mousenames<-c(grep("Musmusculus", colnames(celltypeslist)),
              grep("Tsukui", colnames(celltypeslist)))

#use bioMart to convert the gene lists
#convert the mouse genes to human genes (ensembl and hgnc symbol)
library('biomaRt')

convertMouseGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("mgi_symbol"),
                   filters = "mgi_symbol",
                   values = x ,
                   mart = mouse,
                   attributesL = c( "hgnc_symbol"),
                   martL = human,
                   uniqueRows=T)
  
  colnames(genesV2)<-c("mus_symbol", "hgnc_symbol")
  return(genesV2)
}
celltypeslist_musGenes<-celltypeslist[,mousenames]
#convert mouse to human
celltypeslist_hsGenes<-lapply(celltypeslist_musGenes, convertMouseGeneList)
#turn this into a dataframe of only the human symbols with the same dimensions as the original
celltypeslist_hsGenesonly<-lapply(celltypeslist_hsGenes, `[[`, 2)
m1 <-nrow(celltypeslist)
celltypeslist_hsGenesonlydf<- as.data.frame(do.call(cbind, lapply(celltypeslist_hsGenesonly, `length<-`, m1)), stringsAsFactors= FALSE)
colnames(celltypeslist_hsGenesonlydf)<-gsub("Musmusculus", "Hs", colnames(celltypeslist_hsGenesonlydf))
colnames(celltypeslist_hsGenesonlydf)<-gsub("Tsukui", "Hs_Tsukui", colnames(celltypeslist_hsGenesonlydf))

#merge the human back into the original gene list
celltypeslist_Hs<-celltypeslist[,-mousenames]
celltypeslist_Hs<-data.frame(celltypeslist_Hs,celltypeslist_hsGenesonlydf)

library(msigdbr)
hs_msigdbgenesets_hallmarks = msigdbr(species = "Homo sapiens", category = "H")
hs_msigdbgenesets_c2 = msigdbr(species = "Homo sapiens", category = "C2")

#reformat for gsva so each gene set is in a list of its own in HGNC format
hs_msigdbgenesets_hallmarks = split(x = hs_msigdbgenesets_hallmarks$gene_symbol, f = hs_msigdbgenesets_hallmarks$gs_name)
hs_msigdbgenesets_c2 = split(x = hs_msigdbgenesets_c2$gene_symbol, f = hs_msigdbgenesets_c2$gs_name)

msigdb_h<-hs_msigdbgenesets_hallmarks
msigdb_c2<-hs_msigdbgenesets_c2

msigdb_h<-as.list(msigdb_h)

sigs<-c(as.list(celltypeslist_Hs))
#add some random signatures
set.seed(1234)
for (i in 1:length(sigs))
{
  sigs[[sprintf("Rand%d", i)]] = sample(rownames(LUAD_rnaseq_logcpmTNT), length(sigs[[i]]))
}

colnames(LUSC_rnaseq_logcpmTNT)<-gsub("*._", "", colnames(LUSC_rnaseq_logcpmTNT))

colnames(LUSC_rnaseq_logcpmTNT)<-paste0(colnames(LUSC_rnaseq_logcpmTNT),"_", LUSC_rnaseq_logcpmTNT[1,])

set.seed(20131113)
LUSC_gsva_sigs = gsva(LUSC_rnaseq_logcpmTNT, sigs)
LUSC_gsva_testsigs = gsva(LUSC_rnaseq_logcpmTNT[2:nrow(LUSC_rnaseq_logcpmTNT),], test_sigs)

#Msigdb Hallmarks
LUSC_gsva_msigdb_h = gsva(LUSC_rnaseq_logcpmTNT, msigdb_h)

#Msigdb C2
LUSC_gsva_msigdb_c2 = gsva(LUSC_rnaseq_logcpmTNT, msigdb_c2)

#integrin sigs
integringenes<-rownames(LUSC_rnaseq_logcpmTNT)[grep("^ITG.*", rownames(LUSC_rnaseq_logcpmTNT))]
integringenes<-integringenes[-c(19,20,23,30)]

set.seed(20131113)
LUSC_gsva_integrin = gsva(LUSC_rnaseq_logcpmTNT, list(integringenes))

#Core ECM components of FLT4 ligands
#significant ligands for fibrogenic receptors
FLT4lig<-LUSC_ligrecint_clustpval_ann_ECMupfilt$Ligand.ApprovedSymbol[which(LUSC_ligrecint_clustpval_ann_ECMupfilt$Receptor.ApprovedSymbol=="FLT4")]
PDGFRBlig<-LUSC_ligrecint_clustpval_ann_ECMupfilt$Ligand.ApprovedSymbol[which(LUSC_ligrecint_clustpval_ann_ECMupfilt$Receptor.ApprovedSymbol=="PDGFRB")]

fibrogenic_lig<-list(FLT4lig=FLT4lig, PDGFRBlig=PDGFRBlig)

set.seed(20131113)
LUSC_gsva_corematlig4rec = gsva(LUSC_rnaseq_logcpmTNT, fibrogenic_lig)

#integrins
LUSC_TvNT_sigmatrisome_ccannot2_gsvaord_integ<-data.frame(LUSC_TvNT_sigmatrisome_ccannot2_gsvaord,
                                                          Integrins = t(LUSC_gsva_integrin))
LUSC_TvNT_sigmatrisome_ccannot2_gsvaord_integTonly <- subset(LUSC_TvNT_sigmatrisome_ccannot2_gsvaord_integ,Km_K3!=2)
LUSC_gsva_integrins_wilcox<-wilcox.test(LUSC_TvNT_sigmatrisome_ccannot2_gsvaord_integTonly$Integrins~LUSC_TvNT_sigmatrisome_ccannot2_gsvaord_integTonly$Km_K3,
                                        LUSC_TvNT_sigmatrisome_ccannot2_gsvaord_integTonly) 
LUSC_gsva_integrins<-ggplot(LUSC_TvNT_sigmatrisome_ccannot2_gsvaord_integTonly,
                            aes(x=as.factor(Km_K3), y=Integrins, fill = as.factor(Km_K3)))+
  theme_bw()+
  geom_violin(trim=FALSE) +
  scale_y_continuous(name = "Integrin Score") +
  scale_fill_manual(values = c("royalblue", "forestgreen"))+
  ggtitle(sprintf("%.3e", LUSC_gsva_integrins_wilcox$p.value))+
  xlab("Subtype")+
  geom_boxplot(width = 0.1)+
   theme(axis.text=element_text(size = 10, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 10, face="bold", hjust= 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor= element_blank(),
        legend.position = "none")

#corematrisomal fibrogenic mediators
LUSC_TvNT_sigmatrisome_ccannot2_gsvaord_fibrogen<-data.frame(LUSC_TvNT_sigmatrisome_ccannot2_gsvaord,
                                                             Fibrogenic = t(LUSC_gsva_corematlig4rec))
LUSC_TvNT_sigmatrisome_ccannot2_gsvaord_fibrogenTonly <- subset(LUSC_TvNT_sigmatrisome_ccannot2_gsvaord_fibrogen,Km_K3!=2)
LUSC_gsva_fibrogenFLT4_wilcox<-wilcox.test(LUSC_TvNT_sigmatrisome_ccannot2_gsvaord_fibrogenTonly$Fibrogenic.FLT4lig~LUSC_TvNT_sigmatrisome_ccannot2_gsvaord_fibrogenTonly$Km_K3,
                                           LUSC_TvNT_sigmatrisome_ccannot2_gsvaord_fibrogenTonly) 
LUSC_gsva_fibrogenFLT4<-ggplot(LUSC_TvNT_sigmatrisome_ccannot2_gsvaord_fibrogenTonly,
                               aes(x=as.factor(Km_K3), y=Fibrogenic.FLT4lig, fill = as.factor(Km_K3)))+
  theme_bw()+
  geom_violin(trim=FALSE) +
  scale_y_continuous(name = "FLT4 Ligand Score") +
  scale_fill_manual(values = c("royalblue", "forestgreen"))+
  ggtitle(sprintf("%.3e", LUSC_gsva_fibrogenFLT4_wilcox$p.value))+
  xlab("Subtype")+
  geom_boxplot(width = 0.1)+
  theme(axis.text=element_text(size = 10, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 10, face="bold", hjust= 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor= element_blank(),
        legend.position = "none")

LUSC_gsva_fibrogenPDGFRB_wilcox<-wilcox.test(LUSC_TvNT_sigmatrisome_ccannot2_gsvaord_fibrogenTonly$Fibrogenic.PDGFRBlig~LUSC_TvNT_sigmatrisome_ccannot2_gsvaord_fibrogenTonly$Km_K3,
                                             LUSC_TvNT_sigmatrisome_ccannot2_gsvaord_fibrogenTonly) 
LUSC_gsva_fibrogenPDGFRB<-ggplot(LUSC_TvNT_sigmatrisome_ccannot2_gsvaord_fibrogenTonly,
                                 aes(x=as.factor(Km_K3), y=Fibrogenic.PDGFRBlig, fill = as.factor(Km_K3)))+
  theme_bw()+
  geom_violin(trim=FALSE) +
  scale_y_continuous(name = "PDGFRB Ligand Score") +
  scale_fill_manual(values = c("royalblue", "forestgreen"))+
  ggtitle(sprintf("%.3e", LUSC_gsva_fibrogenPDGFRB_wilcox$p.value))+
  xlab("Subtype")+
  geom_boxplot(width = 0.1)+
  theme(axis.text=element_text(size = 10, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 10, face="bold", hjust= 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor= element_blank(),
        legend.position = "none")


LUSC_gsva_design<-model.matrix(~0+LUSC_TvNT_sigmatrisome_ccannot2_gsvaord$Km_K3)
colnames(LUSC_gsva_design)<-c("K1", "K2", "K3")

LUSC_gsva_fit <- lmFit(LUSC_gsva_sigs, LUSC_gsva_design) #exp matrix with samples in columns
LUSC_gsva_contrast<-makeContrasts(coef1 =K1-K3, 
                                  levels = LUSC_gsva_design)
LUSC_gsva_contrastfit<-contrasts.fit(LUSC_gsva_fit, LUSC_gsva_contrast)
LUSC_gsva_contrastfit<-eBayes(LUSC_gsva_contrastfit)
LUSC_gsva_contrastfitTableK1K3<-topTable(LUSC_gsva_contrastfit, n= 'Inf', coef = 1)

#LUSC Hallmarks

LUSC_gsva_h_fit <- lmFit(LUSC_gsva_msigdb_h, LUSC_gsva_design) #exp matrix with samples in columns
LUSC_gsva_contrast<-makeContrasts(coef1 =K1-K3, 
                                  levels = LUSC_gsva_design)
LUSC_gsva_h_contrastfit<-contrasts.fit(LUSC_gsva_h_fit, LUSC_gsva_contrast)
LUSC_gsva_h_contrastfit<-eBayes(LUSC_gsva_h_contrastfit)
LUSC_gsva_h_contrastfitTableK1K3<-topTable(LUSC_gsva_h_contrastfit, n= 'Inf', coef = 1)

#LUSC C2
LUSC_gsva_c2_fit <- lmFit(LUSC_gsva_msigdb_c2, LUSC_gsva_design) #exp matrix with samples in columns
LUSC_gsva_contrast<-makeContrasts(coef1 =K3-K1, 
                                  levels = LUSC_gsva_design)
LUSC_gsva_c2_contrastfit<-contrasts.fit(LUSC_gsva_c2_fit, LUSC_gsva_contrast)
LUSC_gsva_c2_contrastfit<-eBayes(LUSC_gsva_c2_contrastfit)
LUSC_gsva_c2_contrastfitTableK3K1<-topTable(LUSC_gsva_c2_contrastfit, n= 'Inf', coef = 1)

#LUSC identify which cell types are significantly different between good and poor prognosis
LUSC_gsva_contrastfitTableK1K3_sig<-LUSC_gsva_contrastfitTableK1K3[which(LUSC_gsva_contrastfitTableK1K3$adj.P.Val<0.05),]

LUSC_gsva_h_contrastfitTableK1K3_sig<-LUSC_gsva_h_contrastfitTableK1K3[which(LUSC_gsva_h_contrastfitTableK1K3$adj.P.Val<0.05),]

LUSC_gsva_sigst_clust<-cbind(t(LUSC_gsva_sigs), LUSC_TvNT_sigmatrisome_ccannot2_gsvaord)

#compare clusters for specific signatures
LUSC_gsva_sigst_clustK1K3<-subset(LUSC_gsva_sigst_clust, Km_K3 != 2)
LUSC_gsva_sigst_clustK2K3<-subset(LUSC_gsva_sigst_clust, Km_K3 != 1)
LUSC_gsva_sigst_clustK1K2<-subset(LUSC_gsva_sigst_clust, Km_K3 != 3)

LUSC_gsva_sigst_clust_KeggECM_K1vsK3<-wilcox.test(LUSC_gsva_sigst_clustK1K3$kegg_ecm~LUSC_gsva_sigst_clustK1K3$Km_K3) #p<2.2E-16
LUSC_gsva_sigst_clust_KeggECM_K2vsK3<-wilcox.test(LUSC_gsva_sigst_clustK2K3$kegg_ecm~LUSC_gsva_sigst_clustK2K3$Km_K3) #p=0.2668
LUSC_gsva_sigst_clust_KeggECM_K1vsK2<-wilcox.test(LUSC_gsva_sigst_clustK1K2$kegg_ecm~LUSC_gsva_sigst_clustK1K2$Km_K3) #p=1.33E-10

LUSC_gsva_sigst_clust_reactomeECM_K1vsK3<-wilcox.test(LUSC_gsva_sigst_clustK1K3$reactome_ecm_degradation~LUSC_gsva_sigst_clustK1K3$Km_K3) #p<2.2E-16
LUSC_gsva_sigst_clust_reactomeECM_K2vsK3<-wilcox.test(LUSC_gsva_sigst_clustK2K3$reactome_ecm_degradation~LUSC_gsva_sigst_clustK2K3$Km_K3) #p=0.01115
LUSC_gsva_sigst_clust_reactomeECM_K1vsK2<-wilcox.test(LUSC_gsva_sigst_clustK1K2$reactome_ecm_degradation~LUSC_gsva_sigst_clustK1K2$Km_K3) #p=1.311E-7

LUSC_gsva_sigst_clust_wilcox<-data.frame()
for (i in 1:73){
  LUSC_gsva_sigst_clust_wilcox[i,1]<-wilcox.test(LUSC_gsva_sigst_clustK1K3[,i]~LUSC_gsva_sigst_clustK1K3$Km_K3)$p.value 
  LUSC_gsva_sigst_clust_wilcox[i,2]<-wilcox.test(LUSC_gsva_sigst_clustK2K3[,i]~LUSC_gsva_sigst_clustK2K3$Km_K3)$p.value 
  LUSC_gsva_sigst_clust_wilcox[i,3]<-wilcox.test(LUSC_gsva_sigst_clustK1K2[,i]~LUSC_gsva_sigst_clustK1K2$Km_K3)$p.value 
}
colnames(LUSC_gsva_sigst_clust_wilcox)<-c("K1vsK3_pval", "K2vsK3_pval", "K1vsK2_pval")
rownames(LUSC_gsva_sigst_clust_wilcox)<-colnames(LUSC_gsva_sigst_clust)[1:73]
LUSC_gsva_sigst_clust_pvaladj<-apply(LUSC_gsva_sigst_clust_wilcox, 2, p.adjust)
colnames(LUSC_gsva_sigst_clust_pvaladj)<-c("K1vsK3_pvaladj", "K2vsK3_pvaladj", "K1vsK2_pvaladj")

LUSC_gsva_sigst_clust_wilcox<-cbind(LUSC_gsva_sigst_clust_wilcox,LUSC_gsva_sigst_clust_pvaladj)

LUSC_gsva_contrastfitTable_sig_vio<-lapply(rownames(LUSC_gsva_contrastfitTableK1K3_sig), function(x) {
  maintitle<-paste0("LUSC ", x)
  ggplot(LUSC_gsva_sigst_clust, aes_string(x=as.factor(LUSC_gsva_sigst_clust$Km_K3), y=x, fill = as.factor(LUSC_gsva_sigst_clust$Km_K3)))+
    theme_bw()+
    geom_violin(trim=FALSE) +
    scale_y_continuous(name = x) +
    scale_x_discrete(name = "ECM Molecular Subtype")+
    scale_fill_manual(values = c("royalblue", "gold", "green"))+
    ggtitle(maintitle)+
    geom_boxplot(alpha = 0.7, width = 0.1)+
    labs(fill = "Cluster")+
    theme(axis.text=element_text(size = 10, face="bold"),
          axis.title=element_text(size=10, face="bold"),
          plot.title=element_text(size = 12, face="bold", hjust= 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor= element_blank(),
          legend.position = "none")
})
LUSC_gsva_contrastfitTable_sig_vio[[1]]

#group cells into stacked bar plot
library(reshape2)
LUSC_gsva_sigst_clust_long<-data.frame(Acc = rownames(LUSC_gsva_sigst_clust),
                                       LUSC_gsva_sigst_clust)
LUSC_gsva_sigst_clust_long<-aggregate(LUSC_gsva_sigst_clust_long, by = list(LUSC_gsva_sigst_clust_long$Km_K3), FUN=mean)
colnames(LUSC_gsva_sigst_clust_long)[1]<-"Km_K3"

LUSC_gsva_sigst_clust_long_immune<-melt(LUSC_gsva_sigst_clust_long[,c(1,11:35,77)], id.vars = "Km_K3")

LUSC_gsva_sigst_clust_long_immuneTcell<-LUSC_gsva_sigst_clust_long_immune[grep("T",LUSC_gsva_sigst_clust_long_immune$variable),]
LUSC_gsva_sigst_clust_long_immuneTcell<-droplevels(LUSC_gsva_sigst_clust_long_immuneTcell)

LUSC_gsva_clust_immuneTcell_stacked<-ggplot(LUSC_gsva_sigst_clust_long_immuneTcell, aes(x=as.factor(LUSC_gsva_sigst_clust_long_immuneTcell$Km_K3), y=value, fill = variable))+
  theme_bw()+
  geom_bar(position = "stack", stat = "identity") +
  ggtitle("T Cell Subsets")+
  xlab("Matreotype")+
  ylab("Score")+
  theme(axis.text=element_text(size = 10, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 12, face="bold", hjust= 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor= element_blank(),
        legend.position = "right")

dev.off()

LUSC_gsva_sigst_clust_long_immuneNKcell<-LUSC_gsva_sigst_clust_long_immune[grep("NK",LUSC_gsva_sigst_clust_long_immune$variable),]
LUSC_gsva_sigst_clust_long_immuneNKcell<-droplevels(LUSC_gsva_sigst_clust_long_immuneNKcell)

LUSC_gsva_clust_immuneNKcell_stacked<-ggplot(LUSC_gsva_sigst_clust_long_immuneNKcell, aes(x=as.factor(LUSC_gsva_sigst_clust_long_immuneNKcell$Km_K3), y=value, fill = variable))+
  theme_bw()+
  geom_bar(position = "stack", stat = "identity") +
  ggtitle("NK cell Subsets")+
  xlab("Matreotype")+
  ylab("Score")+
  theme(axis.text=element_text(size = 10, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 12, face="bold", hjust= 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor= element_blank(),
        legend.position = "right")
LUSC_gsva_sigst_clust_long_immuneDCcell<-LUSC_gsva_sigst_clust_long_immune[grep("DC",LUSC_gsva_sigst_clust_long_immune$variable),]
LUSC_gsva_sigst_clust_long_immuneDCcell<-droplevels(LUSC_gsva_sigst_clust_long_immuneDCcell)

LUSC_gsva_clust_immuneDCcell_stacked<-ggplot(LUSC_gsva_sigst_clust_long_immuneDCcell, aes(x=as.factor(LUSC_gsva_sigst_clust_long_immuneDCcell$Km_K3), y=value, fill = variable))+
  theme_bw()+
  geom_bar(position = "stack", stat = "identity") +
  ggtitle("DC Cell Subsets")+
  xlab("Matreotype")+
  ylab("Score")+
  theme(axis.text=element_text(size = 10, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 12, face="bold", hjust= 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor= element_blank(),
        legend.position = "right")

myeloidcells<-c("Eosinophils", "Macrophages", "Mast_cells", "Neutrophils")
LUSC_gsva_sigst_clust_long_immunemyeloidcell<-LUSC_gsva_sigst_clust_long_immune[which(LUSC_gsva_sigst_clust_long_immune$variable %in% myeloidcells),]
LUSC_gsva_sigst_clust_long_immunemyeloidcell<-droplevels(LUSC_gsva_sigst_clust_long_immunemyeloidcell)

LUSC_gsva_clust_immunemyeloidcell_stacked<-ggplot(LUSC_gsva_sigst_clust_long_immunemyeloidcell, aes(x=as.factor(LUSC_gsva_sigst_clust_long_immunemyeloidcell$Km_K3), y=value, fill = variable))+
  theme_bw()+
  geom_bar(position = "stack", stat = "identity") +
  ggtitle("Myeloid Lineage Cells")+
  xlab("Matreotype")+
  ylab("Score")+
  theme(axis.text=element_text(size = 10, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 12, face="bold", hjust= 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor= element_blank(),
        legend.position = "right")

#Bubble plot of the hallmarks comparing K1 and K3
#order the significant hallmarks by logFC high to low
LUSC_gsva_h_contrastfitTableK1K3_sigord<-LUSC_gsva_h_contrastfitTableK1K3_sig[order(LUSC_gsva_h_contrastfitTableK1K3_sig$logFC, decreasing = TRUE),]
LUSC_gsva_h_contrastfitTableK1K3_sigord$X<-gsub("HALLMARK_", "", LUSC_gsva_h_contrastfitTableK1K3_sigord$X)

LUSC_gsva_h_contrastfitTableK3K1_sigord<-LUSC_gsva_h_contrastfitTableK1K3_sigord
LUSC_gsva_h_contrastfitTableK3K1_sigord$logFC<--LUSC_gsva_h_contrastfitTableK3K1_sigord$logFC
rownames(LUSC_gsva_h_contrastfitTableK3K1_sigord)<-LUSC_gsva_h_contrastfitTableK3K1_sigord$X

LUSC_gsva_h_contrastfitTableK3K1_sigord_bubble<-ggplot(LUSC_gsva_h_contrastfitTableK3K1_sigord, aes(x = 1, y = rownames(LUSC_gsva_h_contrastfitTableK1K3_sigord))) + 
  geom_point(aes( color  = logFC, size = -log10(adj.P.Val)), alpha = 0.7, stroke = 0.1) +
  scale_y_discrete(limits = LUSC_gsva_h_contrastfitTableK1K3_sigord$X, name = "Hallmarks", expand=c(0,2))+ #expand=c(multiple, additive)
  scale_color_gradient(low="blue", high="red")+
  scale_size(range = c(0.5, 12)) +  
  theme(axis.text=element_text(size = 10, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 12, face="bold", hjust= 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor= element_blank(),
        legend.position = "right",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        panel.background = element_blank(),
        panel.grid.major.y = element_line( size=.1, color="grey" ))#,

#C2 MsigDb 
#there were a lot of significant hits in the C2 comparison - need to manually choose which to include
LUSC_ECMrows<-grep("ECM", rownames(LUSC_gsva_c2_contrastfitTableK3K1))
LUSC_colrows<-grep("COLLAGEN*", rownames(LUSC_gsva_c2_contrastfitTableK3K1))
LUSC_fibrows<-grep("FIBRONECTIN", rownames(LUSC_gsva_c2_contrastfitTableK3K1))
LUSC_integrinrows<-grep("INTEGRIN*", rownames(LUSC_gsva_c2_contrastfitTableK3K1))
LUSC_matrisomerows<-grep("MATRISOM*", rownames(LUSC_gsva_c2_contrastfitTableK3K1))
LUSC_slrprows<-grep("SLRP*", rownames(LUSC_gsva_c2_contrastfitTableK3K1))
LUSC_proteoglycanrows<-grep("PROTEOGLYCAN*", rownames(LUSC_gsva_c2_contrastfitTableK3K1))

LUSC_gsva_c2_contrastfitTableK3K1_ecm<-LUSC_gsva_c2_contrastfitTableK3K1[c(LUSC_ECMrows,
                                                                           LUSC_colrows,
                                                                           LUSC_fibrows,
                                                                           LUSC_integrinrows,
                                                                           LUSC_matrisomerows,
                                                                           LUSC_slrprows,
                                                                           LUSC_proteoglycanrows),]

#these all look significant so don't need to filter further but do need to order them
LUSC_gsva_c2_contrastfitTableK3K1_ecmord<-LUSC_gsva_c2_contrastfitTableK3K1_ecm[order(LUSC_gsva_c2_contrastfitTableK3K1_ecm$logFC, decreasing = TRUE),]

#look at non-canonical ECM pathways
LUSC_gsva_c2_contrastfitTableK3K1_nonECM<-LUSC_gsva_c2_contrastfitTableK3K1[-c(LUSC_ECMrows,
                                                                               LUSC_colrows,
                                                                               LUSC_fibrows,
                                                                               LUSC_integrinrows,
                                                                               LUSC_matrisomerows,
                                                                               LUSC_slrprows,
                                                                               LUSC_proteoglycanrows),]

#take top 25 up and top 25 down that are significant
LUSC_gsva_c2_contrastfitTableK3K1_sig<-LUSC_gsva_c2_contrastfitTableK3K1_nonECM[which(LUSC_gsva_c2_contrastfitTableK3K1$adj.P.Val<0.01),]

LUSC_gsva_c2_contrastfitTableK3K1_sigup<-LUSC_gsva_c2_contrastfitTableK3K1_sig[which(LUSC_gsva_c2_contrastfitTableK3K1_sig$logFC>0),]
LUSC_gsva_c2_contrastfitTableK3K1_sigdn<-LUSC_gsva_c2_contrastfitTableK3K1_sig[which(LUSC_gsva_c2_contrastfitTableK3K1_sig$logFC<0),]

#order by FC and combine into each
LUSC_gsva_c2_contrastfitTableK3K1_sigupord<-LUSC_gsva_c2_contrastfitTableK3K1_sigup[order(LUSC_gsva_c2_contrastfitTableK3K1_sigup$logFC, decreasing = TRUE),]
LUSC_gsva_c2_contrastfitTableK3K1_sigdnord<-LUSC_gsva_c2_contrastfitTableK3K1_sigdn[order(LUSC_gsva_c2_contrastfitTableK3K1_sigdn$logFC, decreasing = FALSE),]

LUSC_gsva_c2_contrastfitTableK3K1_sigdnord25<-LUSC_gsva_c2_contrastfitTableK3K1_sigdnord[1:25,]

LUSC_gsva_c2_contrastfitTableK3K1_sig25ord<-rbind(LUSC_gsva_c2_contrastfitTableK3K1_sigupord[1:25,],
                                                  LUSC_gsva_c2_contrastfitTableK3K1_sigdnord25[order(LUSC_gsva_c2_contrastfitTableK3K1_sigdnord25$logFC, decreasing = TRUE),])

#TKI responsiveness
imatinibpathways<-c("HAEGERSTRAND_RESPONSE_TO_IMATINIB","WP_IMATINIB_AND_CHRONIC_MYELOID_LEUKEMIA","APPEL_IMATINIB_RESPONSE")
LUSC_gsva_c2_contrastfitTableK3K1_imatinib<-LUSC_gsva_c2_contrastfitTableK3K1[which(rownames(LUSC_gsva_c2_contrastfitTableK3K1) %in% imatinibpathways),]
LUSC_gsva_c2_contrastfitTableK3K1_imatinib_bubble<-ggplot(LUSC_gsva_c2_contrastfitTableK3K1_imatinib, aes(x = 1, y = rownames(LUSC_gsva_c2_contrastfitTableK3K1_imatinib))) + 
  #geom_point(aes(color = logFC, size = -log10(adj.P.Val)), alpha = 0.5) +
  geom_point(aes( color  = logFC, size = -log10(adj.P.Val)), alpha = 0.7, stroke = 0.1) +
  scale_y_discrete(limits = rev(rownames(LUSC_gsva_c2_contrastfitTableK3K1_imatinib)), name = "Imatinib Signatures", expand=c(0,2))+ #expand=c(multiple, additive)
  scale_color_gradient(low="pink", high="red")+
  scale_size(range = c(10, 15)) +  # Adjust the range of points size
  #expand_limits(y=50)+
  theme(axis.text=element_text(size = 10, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 12, face="bold", hjust= 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor= element_blank(),
        legend.position = "right",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major.y = element_line( size=.1, color="grey" ))#,

#instead choose the pathways based on specific parameters from the hallmarks
#cisplatin
LUSC_gsva_c2_contrastfitTableK3K1_cisplatin<-LUSC_gsva_c2_contrastfitTableK3K1[grep("*CISPLATIN*", rownames(LUSC_gsva_c2_contrastfitTableK3K1)),]
LUSC_gsva_c2_contrastfitTableK3K1_cisplatin_up<-LUSC_gsva_c2_contrastfitTableK3K1_cisplatin[grep("*_UP",rownames(LUSC_gsva_c2_contrastfitTableK3K1_cisplatin)),]

LUSC_gsva_c2_contrastfitTableK3K1_cisplatinup_bubble<-ggplot(LUSC_gsva_c2_contrastfitTableK3K1_cisplatin_up, aes(x = 1, y = rownames(LUSC_gsva_c2_contrastfitTableK3K1_cisplatin_up))) + 
  geom_point(aes( color  = logFC, size = -log10(adj.P.Val)), alpha = 0.7, stroke = 0.1) +
  scale_y_discrete(limits = rev(rownames(LUSC_gsva_c2_contrastfitTableK3K1_cisplatin_up)), name = "Cisplatin-Related Signatures", expand=c(0,2))+ #expand=c(multiple, additive)
  scale_color_gradient(low="pink", high="red")+
  scale_size(range = c(10, 15)) + 
  theme(axis.text=element_text(size = 10, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 12, face="bold", hjust= 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor= element_blank(),
        legend.position = "right",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major.y = element_line( size=.1, color="grey" ))

#ECM related signatures  
LUSC_gsva_c2_contrastfitTableK3K1_sigecm_bubble<-ggplot(LUSC_gsva_c2_contrastfitTableK3K1_ecmord, aes(x = 1, y = rownames(LUSC_gsva_c2_contrastfitTableK3K1_ecmord))) + 
  geom_point(aes( color  = logFC, size = -log10(adj.P.Val)), alpha = 0.7, stroke = 0) +
  scale_y_discrete(limits = rev(rownames(LUSC_gsva_c2_contrastfitTableK3K1_ecmord)), name = "ECM-related C2", expand=c(0,2))+ #expand=c(multiple, additive)
  scale_color_gradient(low="purple", high="red")+
  scale_size(range = c(0.5, 12)) +  # Adjust the range of points size
  theme(axis.text=element_text(size = 10, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 12, face="bold", hjust= 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor= element_blank(),
        legend.position = "right",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major.y = element_line( size=.1, color="grey" ))#,

#non-ECM related signatures top 25 up and down
LUSC_gsva_c2_contrastfitTableK3K1_sig25_bubble<-ggplot(LUSC_gsva_c2_contrastfitTableK3K1_sig25ord, aes(x = 1, y = rownames(LUSC_gsva_c2_contrastfitTableK3K1_sig25ord))) + 
  geom_point(aes( color  = logFC, size = -log10(adj.P.Val)), alpha = 0.7, stroke = 0.1) +
  scale_y_discrete(limits = rev(rownames(LUSC_gsva_c2_contrastfitTableK3K1_sig25ord)), name = "top 50 C2", expand=c(0,2))+ #expand=c(multiple, additive)
  scale_color_gradient(low="blue", high="red")+
  scale_size(range = c(0.5, 12)) +
  theme(axis.text=element_text(size = 10, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 12, face="bold", hjust= 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor= element_blank(),
        legend.position = "right",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major.y = element_line( size=.1, color="grey" ))#,

#transcriptional genes of interest
genesofinterest<-c("ZEB1", "ZEB2", "SNAI1", "SNAI2", "SNAI3", "TWIST1", "TWIST2", "CDKN1A", "CDC25A",
                   "FLT4", "PDGFRB","TGFBR2", "FGFR1")
LUSC_rnaseq_logcpmTNT_goi<-LUSC_rnaseq_logcpmTNT[which(rownames(LUSC_rnaseq_logcpmTNT) %in% genesofinterest),]

LUSC_rnaseq_logcpmTNT_goiclust<-data.frame(t(LUSC_rnaseq_logcpmTNT_goi),
                                           LUSC_TvNT_sigmatrisome_ccannot2_gsvaord)
LUSC_rnaseq_logcpmTNT_goiclustTonly<-subset(LUSC_rnaseq_logcpmTNT_goiclust, Km_K3!=2)

LUSC_clust_ZEB1_wilcox<-wilcox.test(LUSC_rnaseq_logcpmTNT_goiclustTonly$ZEB1~LUSC_rnaseq_logcpmTNT_goiclustTonly$Km_K3, LUSC_rnaseq_logcpmTNT_goiclustTonly)
LUSC_clust_ZEB1<-ggplot(LUSC_rnaseq_logcpmTNT_goiclustTonly,
                        aes(x=as.factor(Km_K3), y=ZEB1, fill = as.factor(Km_K3)))+
  theme_bw()+
  geom_violin(trim=FALSE) +
  scale_y_continuous(name = "Relative ZEB1 Expression") +
  scale_fill_manual(values = c("royalblue", "forestgreen"))+
  ggtitle(sprintf("%.3e", LUSC_clust_ZEB1_wilcox$p.value))+
  xlab("Subtype")+
  geom_boxplot(width = 0.1)+
  theme(axis.text=element_text(size = 10, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 10, face="bold", hjust= 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor= element_blank(),
        legend.position = "none")

LUSC_clust_ZEB2_wilcox<-wilcox.test(LUSC_rnaseq_logcpmTNT_goiclustTonly$ZEB2~LUSC_rnaseq_logcpmTNT_goiclustTonly$Km_K3, LUSC_rnaseq_logcpmTNT_goiclustTonly)
LUSC_clust_ZEB2<-ggplot(LUSC_rnaseq_logcpmTNT_goiclustTonly,
                        aes(x=as.factor(Km_K3), y=ZEB2, fill = as.factor(Km_K3)))+
  theme_bw()+
  geom_violin(trim=FALSE) +
  scale_y_continuous(name = "Relative ZEB2 Expression") +
  scale_fill_manual(values = c("royalblue", "forestgreen"))+
  ggtitle(sprintf("%.3e", LUSC_clust_ZEB2_wilcox$p.value))+
  xlab("Subtype")+
  geom_boxplot(width = 0.1)+
  theme(axis.text=element_text(size = 10, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 10, face="bold", hjust= 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor= element_blank(),
        legend.position = "none")

LUSC_clust_TWIST1_wilcox<-wilcox.test(LUSC_rnaseq_logcpmTNT_goiclustTonly$TWIST1~LUSC_rnaseq_logcpmTNT_goiclustTonly$Km_K3, LUSC_rnaseq_logcpmTNT_goiclustTonly)
LUSC_clust_TWIST1<-ggplot(LUSC_rnaseq_logcpmTNT_goiclustTonly,
                          aes(x=as.factor(Km_K3), y=TWIST1, fill = as.factor(Km_K3)))+
  theme_bw()+
  geom_violin(trim=FALSE) +
  scale_y_continuous(name = "Relative ZEB2 Expression") +
  scale_fill_manual(values = c("royalblue", "forestgreen"))+
  ggtitle(sprintf("%.3e", LUSC_clust_TWIST1_wilcox$p.value))+
  xlab("Subtype")+
  geom_boxplot(width = 0.1)+
  theme(axis.text=element_text(size = 10, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 10, face="bold", hjust= 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor= element_blank(),
        legend.position = "none")

LUSC_clust_TWIST2_wilcox<-wilcox.test(LUSC_rnaseq_logcpmTNT_goiclustTonly$TWIST2~LUSC_rnaseq_logcpmTNT_goiclustTonly$Km_K3, LUSC_rnaseq_logcpmTNT_goiclustTonly)
LUSC_clust_TWIST2<-ggplot(LUSC_rnaseq_logcpmTNT_goiclustTonly,
                          aes(x=as.factor(Km_K3), y=TWIST2, fill = as.factor(Km_K3)))+
  theme_bw()+
  geom_violin(trim=FALSE) +
  scale_y_continuous(name = "Relative TWIST2 Expression") +
  scale_fill_manual(values = c("royalblue", "forestgreen"))+
  ggtitle(sprintf("%.3e", LUSC_clust_TWIST2_wilcox$p.value))+
  xlab("Subtype")+
  geom_boxplot(width = 0.1)+
  theme(axis.text=element_text(size = 10, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 10, face="bold", hjust= 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor= element_blank(),
        legend.position = "none")

LUSC_clust_SNAI1_wilcox<-wilcox.test(LUSC_rnaseq_logcpmTNT_goiclustTonly$SNAI1~LUSC_rnaseq_logcpmTNT_goiclustTonly$Km_K3, LUSC_rnaseq_logcpmTNT_goiclustTonly)
LUSC_clust_SNAI1<-ggplot(LUSC_rnaseq_logcpmTNT_goiclustTonly,
                         aes(x=as.factor(Km_K3), y=SNAI1, fill = as.factor(Km_K3)))+
  theme_bw()+
  geom_violin(trim=FALSE) +
  scale_y_continuous(name = "Relative SNAI1 Expression") +
  scale_fill_manual(values = c("royalblue", "forestgreen"))+
  ggtitle(sprintf("%.3e", LUSC_clust_SNAI1_wilcox$p.value))+
  xlab("Subtype")+
  geom_boxplot(width = 0.1)+
  theme(axis.text=element_text(size = 10, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 10, face="bold", hjust= 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor= element_blank(),
        legend.position = "none")

LUSC_clust_SNAI2_wilcox<-wilcox.test(LUSC_rnaseq_logcpmTNT_goiclustTonly$SNAI2~LUSC_rnaseq_logcpmTNT_goiclustTonly$Km_K3, LUSC_rnaseq_logcpmTNT_goiclustTonly)
LUSC_clust_SNAI2<-ggplot(LUSC_rnaseq_logcpmTNT_goiclustTonly,
                         aes(x=as.factor(Km_K3), y=SNAI2, fill = as.factor(Km_K3)))+
  theme_bw()+
  geom_violin(trim=FALSE) +
  scale_y_continuous(name = "Relative SNAI2 Expression") +
  scale_fill_manual(values = c("royalblue", "forestgreen"))+
  ggtitle(sprintf("%.3e", LUSC_clust_SNAI2_wilcox$p.value))+
  xlab("Subtype")+
  geom_boxplot(width = 0.1)+
  theme(axis.text=element_text(size = 10, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 10, face="bold", hjust= 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor= element_blank(),
        legend.position = "none")

LUSC_clust_SNAI3_wilcox<-wilcox.test(LUSC_rnaseq_logcpmTNT_goiclustTonly$SNAI3~LUSC_rnaseq_logcpmTNT_goiclustTonly$Km_K3, LUSC_rnaseq_logcpmTNT_goiclustTonly)
LUSC_clust_SNAI3<-ggplot(LUSC_rnaseq_logcpmTNT_goiclustTonly,
                         aes(x=as.factor(Km_K3), y=SNAI3, fill = as.factor(Km_K3)))+
  theme_bw()+
  geom_violin(trim=FALSE) +
  scale_y_continuous(name = "Relative SNAI3 Expression") +
  scale_fill_manual(values = c("royalblue", "forestgreen"))+
  ggtitle(sprintf("%.3e", LUSC_clust_SNAI3_wilcox$p.value))+
  xlab("Subtype")+
  geom_boxplot(width = 0.1)+
  theme(axis.text=element_text(size = 10, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 10, face="bold", hjust= 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor= element_blank(),
        legend.position = "none")

#CDKN1A
LUSC_clust_CDKN1A_wilcox<-wilcox.test(LUSC_rnaseq_logcpmTNT_goiclustTonly$CDKN1A~LUSC_rnaseq_logcpmTNT_goiclustTonly$Km_K3, LUSC_rnaseq_logcpmTNT_goiclustTonly)
LUSC_clust_CDKN1A<-ggplot(LUSC_rnaseq_logcpmTNT_goiclustTonly,
                          aes(x=as.factor(Km_K3), y=CDKN1A, fill = as.factor(Km_K3)))+
  theme_bw()+
  geom_violin(trim=FALSE) +
  scale_y_continuous(name = "Relative CDKN1A Expression") +
  scale_fill_manual(values = c("royalblue", "forestgreen"))+
  ggtitle(sprintf("%.3e", LUSC_clust_CDKN1A_wilcox$p.value))+
  xlab("Subtype")+
  geom_boxplot(width = 0.1)+
  theme(axis.text=element_text(size = 10, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 10, face="bold", hjust= 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor= element_blank(),
        legend.position = "none")

#CDC25A
LUSC_clust_CDC25A_wilcox<-wilcox.test(LUSC_rnaseq_logcpmTNT_goiclustTonly$CDC25A~LUSC_rnaseq_logcpmTNT_goiclustTonly$Km_K3, LUSC_rnaseq_logcpmTNT_goiclustTonly)
LUSC_clust_CDC25A<-ggplot(LUSC_rnaseq_logcpmTNT_goiclustTonly,
                          aes(x=as.factor(Km_K3), y=CDC25A, fill = as.factor(Km_K3)))+
  theme_bw()+
  geom_violin(trim=FALSE) +
  scale_y_continuous(name = "Relative CDC25A Expression") +
  scale_fill_manual(values = c("royalblue", "forestgreen"))+
  ggtitle(sprintf("%.3e", LUSC_clust_CDC25A_wilcox$p.value))+
  xlab("Subtype")+
  geom_boxplot(width = 0.1)+
  theme(axis.text=element_text(size = 10, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 10, face="bold", hjust= 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor= element_blank(),
        legend.position = "none")

#FLT4
LUSC_clust_FLT4_wilcox<-wilcox.test(LUSC_rnaseq_logcpmTNT_goiclustTonly$FLT4~LUSC_rnaseq_logcpmTNT_goiclustTonly$Km_K3, LUSC_rnaseq_logcpmTNT_goiclustTonly)
LUSC_clust_FLT4<-ggplot(LUSC_rnaseq_logcpmTNT_goiclustTonly,
                        aes(x=as.factor(Km_K3), y=FLT4, fill = as.factor(Km_K3)))+
  theme_bw()+
  geom_violin(trim=FALSE) +
  scale_y_continuous(name = "Relative FLT4 Expression") +
  scale_fill_manual(values = c("royalblue", "forestgreen"))+
  ggtitle(sprintf("%.3e", LUSC_clust_FLT4_wilcox$p.value))+
  xlab("Subtype")+
  geom_boxplot(width = 0.1)+
  theme(axis.text=element_text(size = 10, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 10, face="bold", hjust= 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor= element_blank(),
        legend.position = "none")

#PDGFRB
LUSC_clust_PDGFRB_wilcox<-wilcox.test(LUSC_rnaseq_logcpmTNT_goiclustTonly$PDGFRB~LUSC_rnaseq_logcpmTNT_goiclustTonly$Km_K3, LUSC_rnaseq_logcpmTNT_goiclustTonly)
LUSC_clust_PDGFRB<-ggplot(LUSC_rnaseq_logcpmTNT_goiclustTonly,
                          aes(x=as.factor(Km_K3), y=PDGFRB, fill = as.factor(Km_K3)))+
  theme_bw()+
  geom_violin(trim=FALSE) +
  scale_y_continuous(name = "Relative PDGFRB Expression") +
  scale_fill_manual(values = c("royalblue", "forestgreen"))+
  ggtitle(sprintf("%.3e", LUSC_clust_PDGFRB_wilcox$p.value))+
  xlab("Subtype")+
  geom_boxplot(width = 0.1)+
  theme(axis.text=element_text(size = 10, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 10, face="bold", hjust= 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor= element_blank(),
        legend.position = "none")

#integrins
LUSC_TvNT_sigmatrisome_ccannot2_gsvaord_integ<-data.frame(LUSC_TvNT_sigmatrisome_ccannot2_gsvaord,
                                                          Integrins = t(LUSC_gsva_integrin))
LUSC_TvNT_sigmatrisome_ccannot2_gsvaord_integTonly <- subset(LUSC_TvNT_sigmatrisome_ccannot2_gsvaord_integ,Km_K3!=2)
LUSC_gsva_integrins_wilcox<-wilcox.test(LUSC_TvNT_sigmatrisome_ccannot2_gsvaord_integTonly$Integrins~LUSC_TvNT_sigmatrisome_ccannot2_gsvaord_integTonly$Km_K3,
                                        LUSC_TvNT_sigmatrisome_ccannot2_gsvaord_integTonly) 
LUSC_gsva_integrins<-ggplot(LUSC_TvNT_sigmatrisome_ccannot2_gsvaord_integTonly,
                            aes(x=as.factor(Km_K3), y=Integrins, fill = as.factor(Km_K3)))+
  theme_bw()+
  geom_violin(trim=FALSE) +
  scale_y_continuous(name = "Integrin Score") +
  scale_fill_manual(values = c("royalblue", "forestgreen"))+
  ggtitle(sprintf("%.3e", LUSC_gsva_integrins_wilcox$p.value))+
  xlab("Subtype")+
  geom_boxplot(width = 0.1)+
  theme(axis.text=element_text(size = 10, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 10, face="bold", hjust= 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor= element_blank(),
        legend.position = "none")

# Fig 5E; 6A; 7C; Ligand-receptor interaction ---------------------------------------------
#read n the receptors manually annotated into groups
LigandRec<-read.csv("./LigandReceptor_matrisomeannotated.csv", header = TRUE, fill = TRUE, stringsAsFactors = FALSE)

LigandRec$Receptor.Group<-as.factor(LigandRec$Receptor.Group)

LUSC_Ligandrowindex<-rep(NA, length = nrow(LigandRec))
LUSC_Recrowindex<-rep(NA, length = nrow(LigandRec))

LUSC_logcpm_min<-min(apply(LUSC_rnaseq_logcpmTNT[2:nrow(LUSC_rnaseq_logcpmTNT),],1, min)) #indicates -5.1285 is min value

LUSC_rnaseq_logcpmTNT_offset<-rbind(TvNT = LUSC_rnaseq_logcpmTNT[1,],
                                    LUSC_rnaseq_logcpmTNT[2:nrow(LUSC_rnaseq_logcpmTNT),]+abs(LUSC_logcpm_min)+0.01) #where 0.01 is the offset

for (i in 1:nrow(LigandRec)) {
  LUSC_Ligandrowindex[i]<- ifelse(length(which(LigandRec$Ligand.ApprovedSymbol[i] == rownames(LUSC_rnaseq_logcpmTNT_offset))>0), which(LigandRec$Ligand.ApprovedSymbol[i] == rownames(LUSC_rnaseq_logcpmTNT_offset)), NA)
  LUSC_Recrowindex[i]<- ifelse(length(which(LigandRec$Receptor.ApprovedSymbol[i] == rownames(LUSC_rnaseq_logcpmTNT_offset))>0), which(LigandRec$Receptor.ApprovedSymbol[i] == rownames(LUSC_rnaseq_logcpmTNT_offset)), NA)
}

LigandRec_form<-data.frame(LigandRec,
                           LUSC_ligandrow = LUSC_Ligandrowindex,
                           LUSC_receptorrow = LUSC_Recrowindex)
LUSC_ligrecint<-matrix(NA, nrow = nrow(LigandRec_form), ncol = ncol(LUSC_rnaseq_logcpmTNT_offset))
LUSC_flag<-matrix(NA, nrow = nrow(LigandRec_form), ncol = ncol(LUSC_rnaseq_logcpmTNT_offset))
for (i in 1:nrow(LigandRec_form)) {
  for (j in 1: ncol(LUSC_rnaseq_logcpmTNT_offset)){
    LUSC_ligrecint[i,j]<-LUSC_rnaseq_logcpmTNT_offset[LigandRec_form[i,"LUSC_ligandrow"],j]*LUSC_rnaseq_logcpmTNT_offset[LigandRec_form[i,"LUSC_receptorrow"],j]
    LUSC_flag[i,j]<-ifelse(LUSC_rnaseq_logcpmTNT_offset[LigandRec_form[i,"LUSC_ligandrow"],j]<0 && LUSC_rnaseq_logcpmTNT_offset[LigandRec_form[i,"LUSC_receptorrow"],j]<0, 1,0)
  }
} 
rownames(LUSC_ligrecint)<- LigandRec_form$Pair.Name
colnames(LUSC_ligrecint)<- colnames(LUSC_rnaseq_logcpmTNT_offset)
LUSC_ligrecint<-rbind(LUSC_ligrecint, TvNT = LUSC_rnaseq_logcpmTNT_offset[1,])

#LUSC
LUSC_ligrecintt<-t(LUSC_ligrecint)
Acc<-paste0(substr(rownames(LUSC_ligrecintt),1,12), "_", LUSC_ligrecintt[,ncol(LUSC_ligrecintt)])
LUSC_ligrecintt<-data.frame(Row.names = Acc,
                            LUSC_ligrecintt)
#merge ligand receptor data with cluster info
LUSC_ligrecintt_clust<-merge(LUSC_TvNT_sigmatrisome_ccannot3, LUSC_ligrecintt, by = "Row.names")
rownames(LUSC_ligrecintt_clust)<-LUSC_ligrecintt_clust$Row.names
LUSC_ligrecintt_clust<-LUSC_ligrecintt_clust[colSums(!is.na(LUSC_ligrecintt_clust))>0]

LUSC_ligrecint_clustpval<-data.frame(K1else_pval = rep(NA, length = ncol(LUSC_ligrecintt_clust)),
                                     K2else_pval = rep(NA, length = ncol(LUSC_ligrecintt_clust)),
                                     K3else_pval = rep(NA, length = ncol(LUSC_ligrecintt_clust)),
                                     K2vsK3_pval = rep(NA, length = ncol(LUSC_ligrecintt_clust)),
                                     K1else_FC = rep(NA, length = ncol(LUSC_ligrecintt_clust)),
                                     K2else_FC = rep(NA, length = ncol(LUSC_ligrecintt_clust)),
                                     K3else_FC = rep(NA, length = ncol(LUSC_ligrecintt_clust)),
                                     K2vsK3_FC = rep(NA, length = ncol(LUSC_ligrecintt_clust)))
rownames(LUSC_ligrecint_clustpval)<-colnames(LUSC_ligrecintt_clust)
K2vsK3<-subset(LUSC_ligrecintt_clust, Km_K3 == "2" | Km_K3 == "3")

for (i in 37: (ncol(LUSC_ligrecintt_clust)-1)){
  LUSC_ligrecint_clustpval[i,1]<-t.test(LUSC_ligrecintt_clust[,i]~LUSC_ligrecintt_clust$K1vsnot)$p.value #K1vselse
  LUSC_ligrecint_clustpval[i,2]<-t.test(LUSC_ligrecintt_clust[,i]~LUSC_ligrecintt_clust$K2vsnot)$p.value #K2vselse
  LUSC_ligrecint_clustpval[i,3]<-t.test(LUSC_ligrecintt_clust[,i]~LUSC_ligrecintt_clust$K3vsnot)$p.value #K3vselse
  
  LUSC_ligrecint_clustpval[i,4]<-t.test(K2vsK3[,i]~as.factor(K2vsK3$K3vsnot))$p.value #K2vsK3
  
  LUSC_ligrecint_clustpval[i,5]<-mean(LUSC_ligrecintt_clust[LUSC_ligrecintt_clust$K1vsnot==1,i])/mean(LUSC_ligrecintt_clust[LUSC_ligrecintt_clust$K1vsnot==0,i])#K1vselse
  LUSC_ligrecint_clustpval[i,6]<-mean(LUSC_ligrecintt_clust[LUSC_ligrecintt_clust$K2vsnot==1,i])/mean(LUSC_ligrecintt_clust[LUSC_ligrecintt_clust$K2vsnot==0,i]) #K2vselse
  LUSC_ligrecint_clustpval[i,7]<-mean(LUSC_ligrecintt_clust[LUSC_ligrecintt_clust$K3vsnot==1,i])/mean(LUSC_ligrecintt_clust[LUSC_ligrecintt_clust$K3vsnot==0,i]) #K3vselse
  LUSC_ligrecint_clustpval[i,8]<-mean(K2vsK3[K2vsK3$K3vsnot==1,i])/mean(K2vsK3[K2vsK3$K3vsnot==0,i])  #K2vsK3
  
}
#adjust the pvalues for multiple comparisons
LUSC_ligrecint_clustpval<-data.frame(LUSC_ligrecint_clustpval,
                                     apply(LUSC_ligrecint_clustpval[,1:4], 2, p.adjust))
colnames(LUSC_ligrecint_clustpval)[9:12]<-paste0(colnames(LUSC_ligrecint_clustpval)[1:4], "adj")

#now recombine result with annotated input for visual representation
LUSC_ligrecint_clustpval_ann<-merge(LUSC_ligrecint_clustpval, LigandRec_form, by.x = "row.names", by.y = 4, all=TRUE)

#make the receptor group names more plot-friendly by shortening the long ones
levels(LUSC_ligrecint_clustpval_ann$Receptor.Group)[levels(LUSC_ligrecint_clustpval_ann$Receptor.Group)=="Thyroid Stimulating Hormone Receptor"]<-"TSH Receptor"
levels(LUSC_ligrecint_clustpval_ann$Receptor.Group)[levels(LUSC_ligrecint_clustpval_ann$Receptor.Group)=="Thyrotropin-Releasing Hormone Receptor"]<-"TRH Receptor"

levels(LUSC_ligrecint_clustpval_ann$Receptor.Group)[grep("Receptor", levels(LUSC_ligrecint_clustpval_ann$Receptor.Group))]<-gsub("Receptor","Rec.", levels(LUSC_ligrecint_clustpval_ann$Receptor.Group)[grep("Receptor", levels(LUSC_ligrecint_clustpval_ann$Receptor.Group))])
levels(LUSC_ligrecint_clustpval_ann$Receptor.Group)[grep("Lipoprotein ", levels(LUSC_ligrecint_clustpval_ann$Receptor.Group))]<-"Lipoprotein"

LUSC_ligrecint_clustpval_ann$Category_Lig<-gsub("ECM Glycoproteins", "Glycoproteins", LUSC_ligrecint_clustpval_ann$Category_Lig)

#subset for only those that contain significantly different lig-receptor scores
LUSC_ligrecint_clustpval_ann_K2vsK3sig<-LUSC_ligrecint_clustpval_ann[which(LUSC_ligrecint_clustpval_ann$K2vsK3_pvaladj<0.05),]
LUSC_ligrecint_clustpval_ann_K2vsK3sigup<-LUSC_ligrecint_clustpval_ann[which(LUSC_ligrecint_clustpval_ann$K2vsK3_FC>1 & LUSC_ligrecint_clustpval_ann$K2vsK3_pvaladj<0.05),]
LUSC_ligrecint_clustpval_ann_K2vsK3sigdown<-LUSC_ligrecint_clustpval_ann[which(LUSC_ligrecint_clustpval_ann$K2vsK3_FC<1 & LUSC_ligrecint_clustpval_ann$K2vsK3_pvaladj<0.05),]
#order as descending for plotting purposes
LUSC_ligrecint_clustpval_ann_K2vsK3sigupord<-LUSC_ligrecint_clustpval_ann_K2vsK3sigup[order(LUSC_ligrecint_clustpval_ann_K2vsK3sigup$K2vsK3_FC, decreasing = TRUE),]
LUSC_ligrecint_clustpval_ann_K2vsK3sigupord$Receptor.Group<-reorder(LUSC_ligrecint_clustpval_ann_K2vsK3sigupord$Receptor.Group, LUSC_ligrecint_clustpval_ann_K2vsK3sigupord$K2vsK3_FC)

LUSC_ligrecint_clustpval_ann_K2vsK3sigdownord<-LUSC_ligrecint_clustpval_ann_K2vsK3sigdown[order(LUSC_ligrecint_clustpval_ann_K2vsK3sigdown$K2vsK3_FC, decreasing = TRUE),]
LUSC_ligrecint_clustpval_ann_K2vsK3sigdownord$Receptor.Group<-reorder(LUSC_ligrecint_clustpval_ann_K2vsK3sigdownord$Receptor.Group, LUSC_ligrecint_clustpval_ann_K2vsK3sigdownord$K2vsK3_FC)

library(circlize) #for plotting the circos plot

#generate a dataframe where the links are the mean or max fold change values for that category of pairs
LUSC_ligrec_recgrp_clustK2vsK3up_FC<-data.frame(from = LUSC_ligrecint_clustpval_ann_K2vsK3sigupord$Category_Lig,
                                                to = LUSC_ligrecint_clustpval_ann_K2vsK3sigupord$Receptor.Group,
                                                value = LUSC_ligrecint_clustpval_ann_K2vsK3sigupord$K2vsK3_FC) #value based on avg FC
#if you collapse each group to a single value
LUSC_ligrec_recgrp_clustK2vsK3up_FCmn<-aggregate(value ~ from + to, data = LUSC_ligrec_recgrp_clustK2vsK3up_FC ,FUN = mean, na.rm = TRUE)
LUSC_ligrec_recgrp_clustK2vsK3up_FCmx<-aggregate(value ~ from + to, data = LUSC_ligrec_recgrp_clustK2vsK3up_FC ,FUN = max, na.rm = TRUE)

#define colours based on what you are plotting
grid.col <- c(setNames(rainbow(length(levels(LUSC_ligrec_recgrp_clustK2vsK3up_FC$to))), levels(LUSC_ligrec_recgrp_clustK2vsK3up_FC$to)),#assign rainbow to receptors
              "Collagens" = "royalblue", "Secreted Factors" = "forestgreen", "ECM Glycoproteins" = "gold",
              "ECM Regulators" = "red2", "ECM-affiliated Proteins" = "purple", "Proteoglycans" = "black")

#define colours as universal colours so all plots have the same colour assignments
recgroups<-unique(LUSC_ligrecint_clustpval_ann$Receptor.Group)
LUSC_universal_grid.col <- c(setNames(rainbow(length(levels(recgroups))), levels(recgroups)),#assign rainbow to receptors
                             "Collagens" = "royalblue", "Secreted Factors" = "forestgreen", "Glycoproteins" = "gold",
                             "ECM Regulators" = "red2", "ECM-affiliated Proteins" = "purple", "Proteoglycans" = "black")


LUSC_ligrec_recgrp_clustK2vsK3upcoremat<-subset(LUSC_ligrecint_clustpval_ann_K2vsK3sigupord, Division_Lig=="Core matrisome")
LUSC_ligrec_recgrp_clustK2vsK3upcoremat_FC<-data.frame(from = LUSC_ligrec_recgrp_clustK2vsK3upcoremat$Category_Lig,
                                                       to = LUSC_ligrec_recgrp_clustK2vsK3upcoremat$Receptor.Group,
                                                       value = LUSC_ligrec_recgrp_clustK2vsK3upcoremat$K2vsK3_FC)
LUSC_ligrec_recgrp_clustK2vsK3upcoremat_FCmx<-aggregate(value ~ from + to, data = LUSC_ligrec_recgrp_clustK2vsK3upcoremat_FC ,FUN = max, na.rm = TRUE)


#plot based on max
chordDiagram(LUSC_ligrec_recgrp_clustK2vsK3upcoremat_FCmx,  grid.col = grid.col, annotationTrack = "grid", preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), col = "black", cex = 0.85)
}, bg.border = NA)

LUSC_ligrecint_clustpval_ann_ECMupfilt<-LUSC_ligrecint_clustpval_ann[which(LUSC_ligrecint_clustpval_ann$Ligand.ApprovedSymbol %in% LUSC_KmClust_contrastfitTable_coremat_up$ID),]
LUSC_ligrecint_clustpval_ann_ECMdnfilt<-LUSC_ligrecint_clustpval_ann[which(LUSC_ligrecint_clustpval_ann$Ligand.ApprovedSymbol %in% LUSC_KmClust_contrastfitTable_coremat_dn$ID),]

LUSC_ligrecint_clustpval_ann_ECMupfilt_FC<-data.frame(from = LUSC_ligrecint_clustpval_ann_ECMupfilt$Category_Lig,
                                                      to = LUSC_ligrecint_clustpval_ann_ECMupfilt$Receptor.Group,
                                                      value = LUSC_ligrecint_clustpval_ann_ECMupfilt$K2vsK3_FC)
LUSC_ligrecint_clustpval_ann_ECMdnfilt_FC<-data.frame(from = LUSC_ligrecint_clustpval_ann_ECMdnfilt$Category_Lig,
                                                      to = LUSC_ligrecint_clustpval_ann_ECMdnfilt$Receptor.Group,
                                                      value = LUSC_ligrecint_clustpval_ann_ECMdnfilt$K2vsK3_FC)

LUSC_ligrecint_clustpval_ann_ECMupfilt_FCmx<-aggregate(value ~ from + to, data = LUSC_ligrecint_clustpval_ann_ECMupfilt_FC ,FUN = max, na.rm = TRUE)

LUSC_ligrecint_clustpval_ann_ECMdnfilt_FCmx<-aggregate(value ~ from + to, data = LUSC_ligrecint_clustpval_ann_ECMdnfilt_FC ,FUN = max, na.rm = TRUE)

#Fig 5E
#plot based on max - ECM upregulated
circos.par(circle.margin = 0.2)
chordDiagram(LUSC_ligrecint_clustpval_ann_ECMupfilt_FCmx,  grid.col = LUSC_universal_grid.col, annotationTrack = "grid", preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), col = "black", cex = 0.85)
}, bg.border = NA)

#plot based on max - ECM downregulated
circos.par(circle.margin = 0.2)
chordDiagram(LUSC_ligrecint_clustpval_ann_ECMdnfilt_FCmx,  grid.col = LUSC_universal_grid.col, annotationTrack = "grid", preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), col = "black", cex = 0.85)
}, bg.border = NA)

#Fig 7C plot fibrogenic 
LUSC_ligrecint_clustpval_ann_ECMupfilt_FC_fibrogen<-LUSC_ligrecint_clustpval_ann_ECMupfilt[which(LUSC_ligrecint_clustpval_ann_ECMupfilt$Receptor.ApprovedSymbol %in% c("FLT4", "PDGFRB")),]

LUSC_fibrogenic_grid.col <- c(setNames(rainbow(2*nrow(LUSC_ligrecint_clustpval_ann_ECMupfilt_FC_fibrogen)),
                                       c(LUSC_ligrecint_clustpval_ann_ECMupfilt_FC_fibrogen$Receptor.ApprovedSymbol,
                                         LUSC_ligrecint_clustpval_ann_ECMupfilt_FC_fibrogen$Ligand.ApprovedSymbol)))

LUSC_fibrogenic_grid.colv2 <- c("royalblue", "sky2", "navy", "forestgreen",
                                "red2", "orange", "purple2")
names(LUSC_fibrogenic_grid.colv2)<-c(LUSC_ligrecint_clustpval_ann_ECMupfilt_FC_fibrogen$Ligand.ApprovedSymbol,
                                     unique(LUSC_ligrecint_clustpval_ann_ECMupfilt_FC_fibrogen$Receptor.ApprovedSymbol))

chordDiagram(LUSC_ligrecint_clustpval_ann_ECMupfilt_FC_fibrogen[,c(16,15,9)],  grid.col = LUSC_fibrogenic_grid.col, annotationTrack = "grid", preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), col = "black", cex = 0.85)
}, bg.border = NA)

#integrins rows
integringenes<-rownames(LUSC_rnaseq_logcpmTNT)[grep("^ITG.*", rownames(LUSC_rnaseq_logcpmTNT))]
integringenes<-integringenes[-c(19,20,23,30)]

LUSC_ligrecint_clustpval_ann_ECMupfilt_FC_integrins<-LUSC_ligrecint_clustpval_ann_ECMupfilt[which(LUSC_ligrecint_clustpval_ann_ECMupfilt$Receptor.ApprovedSymbol %in% integringenes),]

#filter by significant interactions only
LUSC_ligrecint_clustpval_ann_ECMupfilt_FC_integrins<-LUSC_ligrecint_clustpval_ann_ECMupfilt_FC_integrins[which(LUSC_ligrecint_clustpval_ann_ECMupfilt_FC_integrins$K2vsK3_pvaladj <0.05),]


LUSC_integrin_grid.col <- c(setNames(rainbow(2*nrow(LUSC_ligrecint_clustpval_ann_ECMupfilt_FC_integrins)),
                                     c(LUSC_ligrecint_clustpval_ann_ECMupfilt_FC_integrins$Receptor.ApprovedSymbol,
                                       LUSC_ligrecint_clustpval_ann_ECMupfilt_FC_integrins$Ligand.ApprovedSymbol)))

#Fig 6A
chordDiagram(LUSC_ligrecint_clustpval_ann_ECMupfilt_FC_integrins[,c(16,15,9)],  grid.col = LUSC_integrin_grid.col, annotationTrack = "grid", preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), col = "black", cex = 0.85)
}, bg.border = NA)

# Fig 5F Ligand-Receptor Pathway Level Analysis ----------------------------------

library(clusterProfiler)

library(org.Hs.eg.db)
library(ensembldb)
require(biomaRt)
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

receptors<-LUSC_ligrecint_clustpval[37:nrow(LUSC_ligrecint_clustpval),1]
receptors<-gsub(".*_", "", receptors)

siggene_up<-unique(LUSC_ligrecint_clustpval_ann_K2vsK3sigup$Receptor.ApprovedSymbol[which(LUSC_ligrecint_clustpval_ann_K2vsK3sigup$K2vsK3_pvaladj<0.05)]) #is a list of genes with a value assigned to them (looks like pvals)
siggene_dn<-unique(LUSC_ligrecint_clustpval_ann_K2vsK3sigdown$Receptor.ApprovedSymbol[which(LUSC_ligrecint_clustpval_ann_K2vsK3sigdown$K2vsK3_pvaladj<0.05)]) #is a list of genes with a value assigned to them (looks like pvals)

siggene_up_entrez = getBM(attributes = c("hgnc_symbol", "entrezgene_id"),
                          filters = "hgnc_symbol",
                          values = siggene_up,
                          mart = human,
                          uniqueRows=T,
                          useCache = FALSE)

siggene_dn_entrez = getBM(attributes = c("hgnc_symbol", "entrezgene_id"),
                          filters = "hgnc_symbol",
                          values = siggene_dn,
                          mart = human,
                          uniqueRows=T,
                          useCache = FALSE)

#using the clusterprofiler package with msigdb
library(msigdbr)

#map the MSigDB Hallmark signatures (gs_name) to genes (entrez and gene symbol)
msigdb_hallmark_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene, human_gene_symbol)
head(msigdb_hallmark_t2g)

#perform overrepresentation analysis with the hallmarks
LUSC_ligrecsigup_msigdbHallmark <- enricher(gene = siggene_up,
                                            #universe = universeofgenes,
                                            TERM2GENE=msigdb_hallmark_t2g[,c(1,3)],
                                            minGSSize = 4,
                                            pvalueCutoff = Inf,
                                            qvalueCutoff = Inf)

LUSC_ligrecsigup_msigdbHallmarkdf<-data.frame(ID = LUSC_ligrecsigup_msigdbHallmark$ID, 
                                              GeneRatio = LUSC_ligrecsigup_msigdbHallmark$GeneRatio,
                                              BgRatio = LUSC_ligrecsigup_msigdbHallmark$BgRatio,
                                              Pvalue = LUSC_ligrecsigup_msigdbHallmark$pvalue,
                                              qvalue = LUSC_ligrecsigup_msigdbHallmark$qvalue,
                                              GeneID = LUSC_ligrecsigup_msigdbHallmark$geneID)


mappedHallmarkID<-vector(mode = "list",nrow(LUSC_ligrecint_clustpval_ann_K2vsK3sigup)) 
mappedHallmarkdescription<-vector(mode = "list",nrow(LUSC_ligrecint_clustpval_ann_K2vsK3sigup)) 

for (i in 1:nrow(LUSC_ligrecint_clustpval_ann_K2vsK3sigup)){
  receptorgene<-LUSC_ligrecint_clustpval_ann_K2vsK3sigup$Receptor.ApprovedSymbol[i]
  
  for (j in 1:nrow(LUSC_ligrecsigup_msigdbHallmark)){
    mappedHallmarkID_yesno<-grep(receptorgene,LUSC_ligrecsigup_msigdbHallmark$geneID[j]) 
    if (length(mappedHallmarkID_yesno)>=1){
      HallmarkID<-LUSC_ligrecsigup_msigdbHallmark$ID[j] 
      mappedHallmarkIDvec<-c(unlist(mappedHallmarkID[[i]]), HallmarkID)
      mappedHallmarkID[[i]]<-mappedHallmarkIDvec
      Hallmark_description<-LUSC_ligrecsigup_msigdbHallmark$Description[j] 
      mappedHallmarkdescriptionvec<-c(unlist(mappedHallmarkdescription[[i]]), Hallmark_description)
      mappedHallmarkdescription[[i]]<-mappedHallmarkdescriptionvec
      
    }
  }
}

#now duplicate each row and assign it to the GOID so that the dataframe is in a long format
LUSC_ligrecint_clustpval_ann_K2vsK3sigup_hallmark<-LUSC_ligrecint_clustpval_ann_K2vsK3sigup

#create the dataframe of the receptor ligand itneractions that you can append the GO terms to
temp<-c()
LUSC_ligrecint_hallmarkcomplete<-c()

for (i in 1:length(mappedHallmarkID)){
  temp<-LUSC_ligrecint_clustpval_ann_K2vsK3sigup_hallmark[rep(i, each = length(mappedHallmarkID[[i]])), ]
  LUSC_ligrecint_hallmarkcomplete<-rbind(LUSC_ligrecint_hallmarkcomplete, temp)
}

#unlist the GO terms to append
mappedHallmarkID_long<-unlist(mappedHallmarkID)
mappedHallmarkdescription_long<-unlist(mappedHallmarkdescription)

#now merge all together
LUSC_ligrecint_hallmarkcomplete$mappedHallmarkID<-mappedHallmarkID_long
LUSC_ligrecint_hallmarkcomplete$mappedHallmarkdescription<-mappedHallmarkdescription_long

#now collapse by hallmarks
#filter first so it's just core ECM
LUSC_ligrecint_hallmarkcomplete_coremat<-subset(LUSC_ligrecint_hallmarkcomplete, Division_Lig == "Core matrisome")

#create simpler dataframe for this
LUSC_ligrec_recgrp_clustK2vsK3upcoremat_FChallmark<-data.frame(from = LUSC_ligrecint_hallmarkcomplete_coremat$Category_Lig,
                                                               to = LUSC_ligrecint_hallmarkcomplete_coremat$mappedHallmarkdescription,
                                                               value = LUSC_ligrecint_hallmarkcomplete_coremat$K2vsK3_FC) #value based on avg FC
#if you collapse each group to a single value
LUSC_ligrec_recgrp_clustK2vsK3upcoremat_FCmnhallmark<-aggregate(value ~ from + to, data = LUSC_ligrec_recgrp_clustK2vsK3upcoremat_FChallmark ,FUN = mean, na.rm = TRUE)
LUSC_ligrec_recgrp_clustK2vsK3upcoremat_FCmxhallmark<-aggregate(value ~ from + to, data = LUSC_ligrec_recgrp_clustK2vsK3upcoremat_FChallmark ,FUN = max, na.rm = TRUE)

LUSC_universal_grid.colhallmark <- c(setNames(rainbow(length(levels(LUSC_ligrec_recgrp_clustK2vsK3upcoremat_FCmnhallmark$to))), levels(LUSC_ligrec_recgrp_clustK2vsK3upcoremat_FCmnhallmark$to)),#assign rainbow to receptors
                                     "Collagens" = "royalblue", "Secreted Factors" = "forestgreen", "ECM Glycoproteins" = "gold","Glycoproteins" = "gold",
                                     "ECM Regulators" = "red2", "ECM-affiliated Proteins" = "purple", "Proteoglycans" = "black")

LUSC_ligrec_recgrp_clustK2vsK3upcoremat_FCmnhallmarkv2<-LUSC_ligrec_recgrp_clustK2vsK3upcoremat_FCmnhallmark
LUSC_ligrec_recgrp_clustK2vsK3upcoremat_FCmnhallmarkv2$to<-as.factor(gsub(".*HALLMARK_","",LUSC_ligrec_recgrp_clustK2vsK3upcoremat_FCmnhallmark$to))
LUSC_universal_grid.colhallmarkv2 <- c(setNames(rainbow(length(levels(LUSC_ligrec_recgrp_clustK2vsK3upcoremat_FCmnhallmarkv2$to))), levels(LUSC_ligrec_recgrp_clustK2vsK3upcoremat_FCmnhallmarkv2$to)),#assign rainbow to receptors
                                       "Collagens" = "royalblue", "Secreted Factors" = "forestgreen", "ECM Glycoproteins" = "gold","Glycoproteins" = "gold",
                                       "ECM Regulators" = "red2", "ECM-affiliated Proteins" = "purple", "Proteoglycans" = "black")

#filter now for significant hallmarks
sighallmarks<-LUSC_ligrecsigup_msigdbHallmarkdf[which(LUSC_ligrecsigup_msigdbHallmarkdf$qvalue<0.2),]
LUSC_ligrec_recgrp_clustK2vsK3upcoremat_FCmnhallmarksig<-LUSC_ligrec_recgrp_clustK2vsK3upcoremat_FCmnhallmark[which(LUSC_ligrec_recgrp_clustK2vsK3upcoremat_FCmnhallmark$to %in% sighallmarks$ID),]
LUSC_ligrec_recgrp_clustK2vsK3upcoremat_FCmxhallmarksig<-LUSC_ligrec_recgrp_clustK2vsK3upcoremat_FCmxhallmark[which(LUSC_ligrec_recgrp_clustK2vsK3upcoremat_FCmxhallmark$to %in% sighallmarks$ID),]

#Figure 5F
#remove"hallmark from the names so they're shorter
LUSC_ligrec_recgrp_clustK2vsK3upcoremat_FCmxhallmarksig$to<-gsub(".*HALLMARK_", "", LUSC_ligrec_recgrp_clustK2vsK3upcoremat_FCmxhallmarksig$to)

circos.par(circle.margin = 0.35)
chordDiagram(LUSC_ligrec_recgrp_clustK2vsK3upcoremat_FCmxhallmarksig,  grid.col = LUSC_universal_grid.colhallmarkv2, annotationTrack = "grid", preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), col = "black", cex = 2)
}, bg.border = NA)
dev.off()
circos.clear()

# Fig 5C,D; 6D; Supp 8B; RPPA Data ---------------------------------------------------------------
#TCGA RPPA data obtained from https://tcpaportal.org/tcpa/download.html
#Level4 data has already been batch corrected
LUSC_rppa<-read.csv("./TCGA-LUSC-L4.csv", sep = ",", header = TRUE)
LUSC_rppa$Sample_ID_format<-gsub("-", "\\.", LUSC_rppa$Sample_ID)
LUSC_rppa$Sample_ID_format<-substr(LUSC_rppa$Sample_ID_format, 1, 12)

LUSC_TvNT_sigmatrisome_ccannot3T$Row.names_format<-gsub("_.*", "", LUSC_TvNT_sigmatrisome_ccannot3T$Row.names)

LUSC_rppa_clust<-merge(LUSC_rppa, LUSC_TvNT_sigmatrisome_ccannot3T, by.x = "Sample_ID_format", by.y = "Row.names_format")

#Differential RPPA data by clusters
ttestout<-matrix(NA, nrow = ncol(LUSC_rppa_clust), ncol = 5)
colnames(ttestout)<-c("K3vsK2_pval", "K2_mean", "K3_mean", "K3minusK2_FC", "K3vsK2_pvaladj")
rownames(ttestout)<-colnames(LUSC_rppa_clust)

for (i in 6:242){
  ttest<-t.test(LUSC_rppa_clust[,i]~ LUSC_rppa_clust$Km_K3)
  ttestout[i,1]<-ttest$p.value #significance test
  ttestout[i,2]<-ttest$estimate[1] #mean for K2
  ttestout[i,3]<-ttest$estimate[2] #mean for K3
  ttestout[i,4]<-mean(LUSC_rppa_clust[which(LUSC_rppa_clust$Km_K3==3),i], na.rm = TRUE)-mean(LUSC_rppa_clust[which(LUSC_rppa_clust$Km_K3==2),i], na.rm = TRUE)
}
ttestout<-as.data.frame(ttestout)
ttestout$K3vsK1_pvaladj<-p.adjust(ttestout$K3vsK2_pval, method = "BH")

#filter out only those that are significantly different
ttestout_sig<-ttestout[which(ttestout$K3vsK2_pvaladj<0.1),]

#p21
LUSC_RPPA_p21_wilcox<-wilcox.test(LUSC_rppa_clust$P21~LUSC_rppa_clust$Km_K3, LUSC_rppa_clust) 
LUSC_RPPA_p21<-ggplot(LUSC_rppa_clust, aes(x=as.factor(Km_K3), y=P21, fill = as.factor(Km_K3)))+
  theme_bw()+
  geom_violin(trim=FALSE) +
  scale_y_continuous(name = "Relative p21 Expression") +
  scale_fill_manual(values = c("royalblue", "forestgreen"))+
  ggtitle(sprintf("%.3e", LUSC_RPPA_p21_wilcox$p.value))+
  xlab("Subtype")+
  geom_boxplot(width = 0.1)+
  theme(axis.text=element_text(size = 10, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 10, face="bold", hjust= 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor= element_blank(),
        legend.position = "none")

#myosinIIA
LUSC_RPPA_myosinIIA_wilcox<-wilcox.test(LUSC_rppa_clust$MYOSINIIA_pS1943~LUSC_rppa_clust$Km_K3, LUSC_rppa_clust) 
LUSC_RPPA_myosinIIA<-ggplot(LUSC_rppa_clust, aes(x=as.factor(Km_K3), y=MYOSINIIA_pS1943, fill = as.factor(Km_K3)))+
  theme_bw()+
  geom_violin(trim=FALSE) +
  scale_y_continuous(name = "Myosin-IIA pS1943") +
  scale_fill_manual(values = c("royalblue", "forestgreen"))+
  ggtitle(sprintf("%.3e", LUSC_RPPA_myosinIIA_wilcox$p.value))+
  xlab("Subtype")+
  geom_boxplot(width = 0.1)+
  theme(axis.text=element_text(size = 10, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 10, face="bold", hjust= 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor= element_blank(),
        legend.position = "none")

#ERK activation
LUSC_RPPA_MAPK_pT202Y204_wilcox<-wilcox.test(LUSC_rppa_clust$MAPK_pT202Y204~LUSC_rppa_clust$Km_K3, LUSC_rppa_clust) 
LUSC_RPPA_MAPK_pT202Y204<-ggplot(LUSC_rppa_clust, aes(x=as.factor(Km_K3), y=MAPK_pT202Y204, fill = as.factor(Km_K3)))+
  theme_bw()+
  geom_violin(trim=FALSE) +
  scale_y_continuous(name = "MAPK pT202/pY204") +
  scale_fill_manual(values = c("royalblue", "forestgreen"))+
  ggtitle(sprintf("%.3e", LUSC_RPPA_MAPK_pT202Y204_wilcox$p.value))+
  xlab("Subtype")+
  geom_boxplot(width = 0.1)+
  theme(axis.text=element_text(size = 10, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 10, face="bold", hjust= 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor= element_blank(),
        legend.position = "none")

LUSC_RPPA_MEK1_pS217S221_wilcox<-wilcox.test(LUSC_rppa_clust$MEK1_pS217S221~LUSC_rppa_clust$Km_K3, LUSC_rppa_clust) 

LUSC_RPPA_MEK1_pS217S221<-ggplot(LUSC_rppa_clust, aes(x=as.factor(Km_K3), y=MEK1_pS217S221, fill = as.factor(Km_K3)))+
  theme_bw()+
  geom_violin(trim=FALSE) +
  scale_y_continuous(name = "MEK1 pS217/pS221") +
  scale_fill_manual(values = c("royalblue", "forestgreen"))+
  ggtitle(sprintf("%.3e", LUSC_RPPA_MEK1_pS217S221_wilcox$p.value))+
  xlab("Subtype")+
  geom_boxplot(width = 0.1)+
  theme(axis.text=element_text(size = 10, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 10, face="bold", hjust= 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor= element_blank(),
        legend.position = "none")

#CHK2
LUSC_RPPA_CHK2_wilcox<-wilcox.test(LUSC_rppa_clust$CHK2~LUSC_rppa_clust$Km_K3, LUSC_rppa_clust) 
LUSC_RPPA_CHK2<-ggplot(LUSC_rppa_clust, aes(x=as.factor(Km_K3), y=CHK2, fill = as.factor(Km_K3)))+
  theme_bw()+
  geom_violin(trim=FALSE) +
  scale_y_continuous(name = "Relative CHK2 Expression") +
  scale_fill_manual(values = c("royalblue", "forestgreen"))+
  ggtitle(sprintf("%.3e", LUSC_RPPA_CHK2_wilcox$p.value))+
  xlab("Subtype")+
  geom_boxplot(width = 0.1)+
  theme(axis.text=element_text(size = 10, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 10, face="bold", hjust= 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor= element_blank(),
        legend.position = "none")

#YB1_pS102
LUSC_RPPA_YB1_pS102_wilcox<-wilcox.test(LUSC_rppa_clust$YB1_pS102~LUSC_rppa_clust$Km_K3, LUSC_rppa_clust) 
LUSC_RPPA_YB1_pS102<-ggplot(LUSC_rppa_clust, aes(x=as.factor(Km_K3), y=YB1_pS102, fill = as.factor(Km_K3)))+
  theme_bw()+
  geom_violin(trim=FALSE) +
  scale_y_continuous(name = "YB1 pS102") +
  scale_fill_manual(values = c("royalblue", "forestgreen"))+
  ggtitle(sprintf("%.3e", LUSC_RPPA_YB1_pS102_wilcox$p.value))+
  xlab("Subtype")+
  geom_boxplot(width = 0.1)+
  theme(axis.text=element_text(size = 10, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 10, face="bold", hjust= 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor= element_blank(),
        legend.position = "none")

#RPPA revisions - additional integrin activation signals
#JNK_pT183Y185
LUSC_RPPA_pJNK_wilcox<-wilcox.test(LUSC_rppa_clust$JNK_pT183Y185~LUSC_rppa_clust$Km_K3, LUSC_rppa_clust) 
LUSC_RPPA_pJNK<-ggplot(LUSC_rppa_clust, aes(x=as.factor(Km_K3), y=JNK_pT183Y185, fill = as.factor(Km_K3)))+
  theme_bw()+
  geom_violin(trim=FALSE) +
  scale_y_continuous(name = "JNK pT183/Y185") +
  scale_fill_manual(values = c("royalblue", "forestgreen"))+
  ggtitle(sprintf("%.3e", LUSC_RPPA_pJNK_wilcox$p.value))+
  xlab("Subtype")+
  geom_boxplot(width = 0.1)+
  theme(axis.text=element_text(size = 10, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 10, face="bold", hjust= 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor= element_blank(),
        legend.position = "none")

#Paxillin
LUSC_RPPA_paxillin_wilcox<-wilcox.test(LUSC_rppa_clust$PAXILLIN~LUSC_rppa_clust$Km_K3, LUSC_rppa_clust) 
LUSC_RPPA_paxillin<-ggplot(LUSC_rppa_clust, aes(x=as.factor(Km_K3), y=PAXILLIN, fill = as.factor(Km_K3)))+
  theme_bw()+
  geom_violin(trim=FALSE) +
  scale_y_continuous(name = "PAXILLIN") +
  scale_fill_manual(values = c("royalblue", "forestgreen"))+
  ggtitle(sprintf("%.3e", LUSC_RPPA_paxillin_wilcox$p.value))+
  xlab("Subtype")+
  geom_boxplot(width = 0.1)+
  theme(axis.text=element_text(size = 10, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 10, face="bold", hjust= 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor= element_blank(),
        legend.position = "none")

#PKCdelta
LUSC_RPPA_PKCdelta_pS664_wilcox<-wilcox.test(LUSC_rppa_clust$PKCDELTA_pS664~LUSC_rppa_clust$Km_K3, LUSC_rppa_clust) 
LUSC_RPPA_PKCdelta_pS664<-ggplot(LUSC_rppa_clust, aes(x=as.factor(Km_K3), y=PKCDELTA_pS664, fill = as.factor(Km_K3)))+
  theme_bw()+
  geom_violin(trim=FALSE) +
  scale_y_continuous(name = "PKCdelta pS664") +
  scale_fill_manual(values = c("royalblue", "forestgreen"))+
  ggtitle(sprintf("%.3e", LUSC_RPPA_PKCdelta_pS664_wilcox$p.value))+
  xlab("Subtype")+
  geom_boxplot(width = 0.1)+
  theme(axis.text=element_text(size = 10, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 10, face="bold", hjust= 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor= element_blank(),
        legend.position = "none")

#EGFR-pY1068
LUSC_RPPA_EGFR_pY1068_wilcox<-wilcox.test(LUSC_rppa_clust$EGFR_pY1068~LUSC_rppa_clust$Km_K3, LUSC_rppa_clust) 
LUSC_RPPA_EGFR_pY1068<-ggplot(LUSC_rppa_clust, aes(x=as.factor(Km_K3), y=EGFR_pY1068, fill = as.factor(Km_K3)))+
  theme_bw()+
  geom_violin(trim=FALSE) +
  scale_y_continuous(name = "EGFR pY1068") +
  scale_fill_manual(values = c("royalblue", "forestgreen"))+
  ggtitle(sprintf("%.3e", LUSC_RPPA_EGFR_pY1068_wilcox$p.value))+
  xlab("Subtype")+
  geom_boxplot(width = 0.1)+
  theme(axis.text=element_text(size = 10, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 10, face="bold", hjust= 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor= element_blank(),
        legend.position = "none")

#EGFR-pY1173
LUSC_RPPA_EGFR_pY1173_wilcox<-wilcox.test(LUSC_rppa_clust$EGFR_pY1173~LUSC_rppa_clust$Km_K3, LUSC_rppa_clust) 
LUSC_RPPA_EGFR_pY1173<-ggplot(LUSC_rppa_clust, aes(x=as.factor(Km_K3), y=EGFR_pY1173, fill = as.factor(Km_K3)))+
  theme_bw()+
  geom_violin(trim=FALSE) +
  scale_y_continuous(name = "EGFR pY1173") +
  scale_fill_manual(values = c("royalblue", "forestgreen"))+
  ggtitle(sprintf("%.3e", LUSC_RPPA_EGFR_pY1173_wilcox$p.value))+
  xlab("Subtype")+
  geom_boxplot(width = 0.1)+
  theme(axis.text=element_text(size = 10, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 10, face="bold", hjust= 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor= element_blank(),
        legend.position = "none")

# Supp Fig 8A integrin adhesome --------------------------------------------

#read in integrin adhesome list of proteins
integrin_adhesome<-read.csv(file= "./Horton_etal_ConsensusAdhesome_41556_2015_BFncb3257_MOESM10_ESM_formatted.csv", header = T)

#reads in with the first row blank due to merged header cells
integrin_adhesome<-integrin_adhesome[-1,]

#reads in with extra rows
integrin_adhesome<-integrin_adhesome[-c(61:64),]

#use GSVA to apply a signature to the expression matrix
set.seed(20131113)
LUSC_gsva_integrin_adhesome = gsva(LUSC_rnaseq_logcpmTNT, list(integrin_adhesome$Gene.name))

LUSC_TvNT_sigmatrisome_ccannot2_gsvaord_integ<-data.frame(LUSC_TvNT_sigmatrisome_ccannot2_gsvaord_integ,
                                                          Integrin_adhesome = t(LUSC_gsva_integrin_adhesome))

LUSC_TvNT_sigmatrisome_ccannot2_gsvaord_integTonly <- subset(LUSC_TvNT_sigmatrisome_ccannot2_gsvaord_integ,Km_K3!=2)
LUSC_gsva_integrin_adhesome_wilcox<-wilcox.test(LUSC_TvNT_sigmatrisome_ccannot2_gsvaord_integTonly$Integrin_adhesome~LUSC_TvNT_sigmatrisome_ccannot2_gsvaord_integTonly$Km_K3,
                                                LUSC_TvNT_sigmatrisome_ccannot2_gsvaord_integTonly) 
LUSC_gsva_integrin_adhesome<-ggplot(LUSC_TvNT_sigmatrisome_ccannot2_gsvaord_integTonly,
                                    aes(x=as.factor(Km_K3), y=Integrin_adhesome, fill = as.factor(Km_K3)))+
  theme_bw()+
  geom_violin(trim=FALSE) +
  scale_y_continuous(name = "Integrin Adhesome Score") +
  scale_fill_manual(values = c("royalblue", "forestgreen"))+
  ggtitle(sprintf("%.3e", LUSC_gsva_integrin_adhesome_wilcox$p.value))+
  xlab("Subtype")+
  geom_boxplot(width = 0.1)+
  theme(axis.text=element_text(size = 10, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 10, face="bold", hjust= 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor= element_blank(),
        legend.position = "none")

# analysis of CNV and MAF files ---------------------------------------------------
###FGFR amplification
#FGFR1 is amplified at chr8:38 387 813 - 38 445 509

#read CNV (Hg18) for FGF annotations
LUSC_cnv<-read.table(file = "./LUSC.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg18__seg.seg.txt", sep = "\t", header = TRUE)

LUSC_cnv_chr8<-LUSC_cnv[which(LUSC_cnv$Chromosome=="8"),]
segment_starts<-as.numeric(levels(as.factor(LUSC_cnv_chr8$Start)))

toohigh<-which(segment_starts>38387813)[1]
toolow<-which(segment_starts<38387813)
toolow<-toolow[length(toolow)]

startposition<-segment_starts[toolow:toohigh][1]

LUSC_cnv_chr8_FGFR1seg<-LUSC_cnv_chr8[which(LUSC_cnv_chr8$Start==startposition),] #diesnt work because segments are different for each patient

#split Ch8 into a list by patients
LUSC_cnv_chr8$Sample<-as.factor(LUSC_cnv_chr8$Sample)
LUSC_cnv_ch8_patient<-split(LUSC_cnv_chr8, f=as.factor(LUSC_cnv_chr8$Sample))

AllFGFRsegments<-data.frame(matrix(NA,nrow = length(LUSC_cnv_ch8_patient), ncol = 6))

sample<-c()
chromosome<-c()
start<-c()
ends<-c()
Num_probes<-c()
segment_mean<-c()

for (i in length(LUSC_cnv_ch8_patient)-1){
  toolow<-which(LUSC_cnv_ch8_patient[[i]]$Start<38387813)
  toolow<-toolow[length(toolow)]
  
  startposition<-LUSC_cnv_ch8_patient[[i]]$Start[toolow][1]
  FGFRsegment<-as.data.frame(LUSC_cnv_ch8_patient[[i]][which(LUSC_cnv_ch8_patient[[i]]$Start==startposition),])
  FGFRsegment<-unname(FGFRsegment)
  sample[i]<-FGFRsegment[,1]
  chromosome[i]<-FGFRsegment[,2]
  start[i]<-FGFRsegment[,3]
  ends[i]<-FGFRsegment[,4]
  Num_probes[i]<-FGFRsegment[,5]
  segment_mean[i]<-FGFRsegment[,6]
  
}
colnames(AllFGFRsegments)<-colnames(LUSC_cnv_chr8)

AllFGFRsegments<-lapply(LUSC_cnv_ch8_patient, function(i){
  toolow<-which(i$Start<38387813)
  toolow<-toolow[length(toolow)]
  
  startposition<-i$Start[toolow][1]
  FGFRsegment<-as.data.frame(i[which(i$Start==startposition),])
})
AllFGFRsegmentsdf<-bind_rows(AllFGFRsegments)

hist(AllFGFRsegmentsdf$Segment_Mean, breaks = 20)

#from teh distribution consider <-0.7 as deletion and >0.7 as gain
AllFGFRsegmentsdf$assessment<-rep("WT", length = nrow(AllFGFRsegmentsdf))
AllFGFRsegmentsdf$assessment[which(AllFGFRsegmentsdf$Segment_Mean<=-0.7)]<-"Loss"
AllFGFRsegmentsdf$assessment[which(AllFGFRsegmentsdf$Segment_Mean>=0.7)]<-"Gain"

#merge FGFR1 CNV with cluster info

#tumor only 
AllFGFRsegmentsdf$Acc<-substr(AllFGFRsegmentsdf$Sample,1,17)
AllFGFRsegmentsdfT<-AllFGFRsegmentsdf[grep("-01[A-Z]-", AllFGFRsegmentsdf$Acc),]

AllFGFRsegmentsdfT$Acc<-substr(AllFGFRsegmentsdfT$Acc,1,12)

AllFGFRsegmentsdfT$Acc<-gsub("-", "\\.",AllFGFRsegmentsdfT$Acc)

AllFGFRsegmentsdfT_clust<-merge(AllFGFRsegmentsdfT, LUSC_OS_cc2, by.x = "Acc", by.y = "row.names")

#create a contingency table for cluster 1 vs 3 and Gain vs WT
AllFGFRsegmentsdfT_clustK1K3<-subset(AllFGFRsegmentsdfT_clust, Km_K3 !=2)
FGFR1_clust_table<-table(AllFGFRsegmentsdfT_clustK1K3$assessment, AllFGFRsegmentsdfT_clustK1K3$Km_K3)

FGFR1_clust_table_fishers<-fisher.test(FGFR1_clust_table) #pval = 0.1935 so not significant enriched for FGFR amplification in K3 vs K1

#WES data
#LUSC maf files have been merged in Terminal
source("http://bioconductor.org/biocLite.R")
BiocManager::install("ComplexHeatmap")
BiocManager::install("VariantAnnotation")
BiocManager::install("Biostrings")

library(maftools)

#Individual MAF files downloaded from GDAC and concatenated into one merged file using terminal
#manually reformated the barcode numbers in excel to match maf file format
clin<-read.csv(file = './LUSC_AllMatrisome_KMeansClusterAnnotations.csv', header = TRUE)
LUSCwes<- read.maf(maf = './TCGALUSC_merged.maf.txt')

barcodes<-as.character(getSampleSummary(LUSCwes)$Tumor_Sample_Barcode)

#manually regenerate the cluster information in Excel to match the DNA barcodes
LUSC_clusterfileD<-read.csv('./LUSC_concensusclusterannotations_210114_DNAbarcodeannotations.csv')
colnames(LUSC_clusterfileD)[1]<-"Tumor_Sample_Barcode" #need the ID columns to match between the maf file and the clinical data

#now need to re-order the clinical dataset so that the Accs are in the same order as the WES files

#manually formatted the cluster information in excel to match the RNAseq and WES barcode details
LUSC_clusterfileDord<-LUSC_clusterfileD[match(barcodes,LUSC_clusterfileD$Tumor_Sample_Barcode),]
LUSC_clusterfileDord$Km_K3<-as.factor(paste0("K", LUSC_clusterfileDord$Km_K3))

#add to this the FGFR status
AllFGFRsegmentsdfT_clustK1K3
LUSC_clusterfileDordv2<-merge(LUSC_clusterfileDord, AllFGFRsegmentsdfT_clustK1K3[,c(1,8,47:49)], by="Acc")
LUSC_clusterfileDordv2<-LUSC_clusterfileDordv2[match(barcodes,LUSC_clusterfileDordv2$Tumor_Sample_Barcode),]

LUSCwes_clustv2<- read.maf(maf = './TCGALUSC_190821_merged.maf.txt',
                           clinicalData = LUSC_clusterfileDordv2) #clinical data with FGFR status from below

#note the flags among the top ten genes: TTN, MUC16, USH2A, SYNE1

LUSCwes_clust #gives the summary information
getSampleSummary(LUSCwes_clust) #shows the sample summary data
getGeneSummary(LUSCwes_clust) #for gene level summary
getClinicalData(LUSCwes_clust) #clinical data associated with samples
getFields(LUSCwes_clust) #shows the data types saved in the maf object

getClinicalData(LUSCwes_clustv2)
#plot maf summary
plotmafSummary(maf = LUSCwes_clust, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

#draw oncoplot
oncoplot(maf = LUSCwes_clust, top =20)
#from Campbell et al the major LUSC genes are
LUSC_mutatedgenes<-c("TP53", "CDKN2A", "NFE2L2", "PTEN", "MLL2", "RB1", "FAT1", "NOTCH1", "RASA1", "NF1", "ARID1A", "KDM6A", "PIK3CA", "CUL3", "HRAS", "IRF6", "FBXW7", "ARHGAP35", "PASK", "NSD")
oncostrip(maf = LUSCwes_clust, genes = c("TP53", "KRAS", "KEAP1", "COL11A1", "NFE2L2", "CUL3"), clinicalFeatures = 'Km_K3')
oncostrip(maf = LUSCwes_clust, genes = LUSC_mutatedgenes, clinicalFeatures = 'Km_K3')

#recolor the K1 and K3 values to green and blue
Km_K3colors<-LUSC_colour2$Cluster_Km
Km_K3colors[3]<-"forestgreen"
FGFRcolors<-c("WT"="black", "Gain" = "red2","Loss" = "royalblue")
names(Km_K3colors)<-paste0("K", names(Km_K3colors))
clustcolors = list(Km_K3 = Km_K3colors, assessment = FGFRcolors)

LUSC_mutatedgenes<-c("TP53", "CDKN2A", "NFE2L2", "PTEN", "MLL2", "RB1", "FAT1", "NOTCH1", "RASA1", "NF1", "ARID1A", "KDM6A", "PIK3CA", "CUL3", "HRAS", "IRF6", "FBXW7", "ARHGAP35", "PASK", "NSD")

print(oncostrip(maf = LUSCwes_clustv2, draw_titv=TRUE,genes = LUSC_mutatedgenes, clinicalFeatures = c('Km_K3','assessment'),annotationColor = clustcolors, sortByAnnotation=TRUE))

#generate a second plot that also contains the most highly differentially association mutations
#re-perform the association just in case the old one below is incorrect
LUSCwes_clustv2KmK3<-clinicalEnrichment(LUSCwes_clustv2, clinicalFeature = 'Km_K3')
LUSCwes_clustv2KmK3_paircomp<-LUSCwes_clustv2KmK3$pairwise_comparision
LUSCwes_clustv2KmK3_grpcomp<-LUSCwes_clustv2KmK3$groupwise_comparision

#select signicant genes by fdr in K1 vs K3
LUSCwes_clustv2KmK3_paircompsig<-LUSCwes_clustv2KmK3_paircomp[which(fdr<0.05),]

#generate a plot of significant genes
print(oncostrip(maf = LUSCwes_clustv2, genes = LUSCwes_clustv2KmK3_paircompsig$Hugo_Symbol,
                clinicalFeatures = c('Km_K3','assessment'),annotationColor = clustcolors, sortByAnnotation=TRUE,
                fontSize = 0.5))

#generate a plot of significant genes with fdr<0.01
print(oncostrip(maf = LUSCwes_clustv2, genes = LUSCwes_clustv2KmK3_paircompsig$Hugo_Symbol[1:16],
                clinicalFeatures = c('Km_K3','assessment'),annotationColor = clustcolors, sortByAnnotation=TRUE))
LUSCwes_clust_ascsv<-subsetMaf(maf = LUSCwes_clust, mafObj = FALSE)

#perform mutational enrichment analysis on K5Absorb cluster
LUSCwes_clustKmK3<-clinicalEnrichment(LUSCwes_clust, clinicalFeature = 'Km_K3')
LUSCwes_clustKmK3_paircomp<-LUSCwes_clustKmK3$pairwise_comparision
LUSCwes_clustKmK3_grpcomp<-LUSCwes_clustKmK3$groupwise_comparision

#compare mutational burden between the two clusters
LUSCwes_clustv2_titv = titv(maf = LUSCwes_clustv2, plot = FALSE, useSyn = TRUE)

#sum the raw counts for mutations
LUSCwes_clustv2_titv_rawcounts<-LUSCwes_clustv2_titv$raw.counts
LUSCwes_clustv2_titv_rawcounts$sum<-apply(LUSCwes_clustv2_titv_rawcounts[,2:7],1,sum)

#merge with the cluster info
LUSCwes_clustv2_titv_rawcounts<-merge(LUSCwes_clustv2_titv_rawcounts, LUSC_clusterfileDordv2, by = 'Tumor_Sample_Barcode')
LUSCwes_clustv2_tmbclust_wilcox<-wilcox.test(LUSCwes_clustv2_titv_rawcounts$sum~ LUSCwes_clustv2_titv_rawcounts$Km_K3) #p=0.01401

#violin plot of TMB for clusters
LUSC_tmbbyclust_vio<-ggplot(LUSCwes_clustv2_titv_rawcounts,
                            aes(x=as.factor(Km_K3), y=sum, fill = as.factor(Km_K3)))+
  theme_bw()+
  geom_violin(trim=FALSE) +
  scale_y_continuous(name = "Total Number of Mutations") +
  scale_x_discrete(name = "", limits = c("K1", "K3"), labels = c("ECM-Low", "ECM-High"))+
  scale_fill_manual(values = c("royalblue", "forestgreen"))+
  ggtitle(sprintf("%.3e", LUSCwes_clustv2_tmbclust_wilcox$p.value))+
  xlab("Matreotype")+
  geom_boxplot(width = 0.1)+
  theme(axis.text=element_text(size = 10, face="bold"),
        axis.title=element_text(size=10, face="bold"),
        plot.title=element_text(size = 10, face="bold", hjust= 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor= element_blank(),
        legend.position = "none")

#plot titv summary
plotTiTv(res = LUSCwes_clustv2_titv)

#check if there is a difference in transversion or transition frequencies between clusters
LUSCwes_clustv2_titv_clust<-merge(LUSCwes_clustv2_titv$TiTv.fractions,LUSC_clusterfileDordv2, by = 'Tumor_Sample_Barcode') 
LUSCwes_clustv2_ti_wilcox<-wilcox.test(LUSCwes_clustv2_titv_clust$Ti~ LUSCwes_clustv2_titv_clust$Km_K3) #p=0.733
LUSCwes_clustv2_tv_wilcox<-wilcox.test(LUSCwes_clustv2_titv_clust$Tv~ LUSCwes_clustv2_titv_clust$Km_K3) #p=0.734
#no significant difference between clusters.

