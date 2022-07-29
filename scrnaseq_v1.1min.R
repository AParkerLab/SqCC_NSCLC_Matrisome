#scRNAseq analysis of SqCC Tumours

#Dataset:
#Lambrechts et al. Phenotypic Molding of stromal cells in the lung tumor microenvironment
#Nature Medicine 24(8):1277-1289

#have downloaded loom files from 
#https://gbiomed.kuleuven.be/english/research/50000622/laboratories/54213024/scRNAseq-NSCLC

# establish environment ----------------------------------------------------
library(knitr)
library(renv)
library(loomR)
library(hdf5r)
library(Seurat)

# load all cells counts ---------------------------------------------------


# load in data ------------------------------------------------------
all.data <- Read10X(
  data.dir = file.path("./lung_scRNAseq/2096-Lungcancer_counts/", "LC_counts/")
)
# Initialize the Seurat object with the raw (non-normalized data).
all <- CreateSeuratObject(
  counts = all.data, 
  project = "allcells", 
  min.cells = 3, 
  min.features = 201
)

# QC metrics --------------------------------------------------------
#percentMT genes
all[["percent.mt"]] <- Seurat::PercentageFeatureSet(all, pattern = "^MT-")
# Show QC metrics for the first 5 cells
head(all@meta.data, 5)

# Visualize QC metrics as a violin plot
plot_list <- VlnPlot(
  all, 
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), 
  ncol = 3
)
plot_list
typeof(plot_list)
class(plot_list)
plot_list[[1]] + ggtitle('Feature RNA count') #gives plot by identity

plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


# Filter out low and high gene expression and MT content ---------------------------------
all <- subset(
  all, 
  subset = nFeature_RNA > 101 & nFeature_RNA < 6000 & percent.mt < 10
)


# normalise data ----------------------------------------------------------
#first need to normalise the data then can scale to regress ot total cellular count and mito read count
all <- NormalizeData(all)

#saveRDS(all, "./allcells_seuratobj_norm_211219.rds")
#all<-readRDS("./allcells_seuratobj_norm_211219.rds")
# scale data --------------------------------------------------------------
all <- ScaleData(all, vars.to.regress =c("nCount_RNA","percent.mt"), features = rownames(all))

# select variably expressed genes -----------------------------------------

all<-FindVariableFeatures(all,
                          selection.method = "vst",
                          mean.cutoff = c(0.125,3))

#identify the top 10 most variable genes
top10 <- head(
  VariableFeatures(all), 
  10
)

#run PCA on these variable genes
all <- RunPCA(
              all, 
              features = VariableFeatures(object = all),
              npcs = 41
              )

#visualse the PCA results
print(all[["pca"]], dims = 1:5, nfeatures = 5)
#visualise the dimensional loadings
VizDimLoadings(all, dims = 1:2, reduction="pca")
#DimPlot
DimPlot(all, reduction = "pca")

#determine the dimensionality of the dataset
ElbowPlot(all) #8 would be appropriate

#perform a second dimensional reduction using the runTSNE function and default settings
all<-RunTSNE(all, dims = 1:8, seed.use = 8482)

DimPlot(all, reduction = 'tsne', pt.size = 1) +ggtitle("tSNE with default Perplexity (30)")
#looks similar to the manuscript but not completely the same

# cluster the cells with 8 dimensions -------------------------------------
all<-FindNeighbors(all, dims = 1:8)
all@graphs$RNA_snn[1:10,1:10]


all <-FindClusters(all, resolution = 0.5) 
DimPlot(all, reduction = 'tsne')


# Automatically Define Cell Types -----------------------------------------
all.markers<-FindAllMarkers(all, only.pos = TRUE, min.pct = 0.15, logfc.threshold = 1.32)
all.markers %>%
  group_by(cluster) %>%
  slice_max(n=2, order_by = avg_log2FC)

#see how it looks using cluster feature genes in Lambrechts et al.

VlnPlot(all, features = c("CLDN18", "FOLR1", "AQP4", "PEBP4"), pt.size = FALSE) #alveolar =cluster 7
VlnPlot(all, features = c("CLDN5", "FLT1", "CDH5", "RAMP2"), pt.size = FALSE) #endothelial =cluster 8
VlnPlot(all, features = c("CAPS", "TMEM190", "PIFO", "SNTN"), pt.size = FALSE) #epithelial =cluster 18
VlnPlot(all, features = c("COL1A1", "DCN", "COL1A2", "C1R"), pt.size = FALSE) #fibroblast = cluster 13
VlnPlot(all, features = c("CD79A", "IGKC", "IGLC3", "IGHG3"), pt.size = FALSE) #Bcell = cluster 12 but a bit messy
VlnPlot(all, features = c("LYZ", "MARCO", "CD68", "FCGR3A"), pt.size = FALSE) #myeloid = clusters 2,6,16 but a bit messy
VlnPlot(all, features = c("CD3D", "TRBC1", "TRBC2", "TRAC"), pt.size = FALSE) #T-cell = clusters 0,1,3,4,14,but a bit messy



# reassign cluster IDs ----------------------------------------------------
all@meta.data$seurat_clusters #contains the existing seurat (numbered clusters)

all[["Original_Seurat_Clusters"]]<-Idents(all)
all<-StashIdent(all, save.name = "Original_Seurat_Clusters")
#now reassign the ids based on my manual annotations
cell.labels<-all@ident

new.cluster.ids<-c("Tcell_1", "Tcell_2", "Myeloid_1","Tcell_3", "Tcell_4","Unknown_1","Tcell_4",
                   "Alveolar", "Endothelial", "Unknown_2", "Unknown_3", "Unknown_4", "Bcell",
                   "Fibroblast", "Tcell_5", "Unknown_5", "Myeloid_2", "Unknown_6", "Epithelial" )
names(new.cluster.ids)<-levels(all)
all<-RenameIdents(all, new.cluster.ids)
#this has now stored this data under 'active.ident' so can be changed out for whatever is necessary

#create a new metadata slot for these identities and saved them just in case
all[["Manual_Clusters_v1"]]<-Idents(all)
all<-StashIdent(all, save.name = "Manual_Clusters_v1")

DimPlot(all, reduction = 'tsne', label = TRUE) +NoLegend()

# Extract expression data for manual merging ------------------------------

all_cellinfo<-all@active.ident
all_cellinfodf<-data.frame(CellID = names(all@active.ident),
                           ManualClusterID = all@active.ident,
                           SeuratClusterID = all@meta.data$Original_Seurat_Clusters)

all_scaledata<-all@assays[["RNA"]]@scale.data

#merge into one dataframe and save
all_scaledatat<-t(all_scaledata)
alltogether<-merge(all_cellinfodf, all_scaledatat, by.x = 'CellID', by.y = 'row.names')

#read in existing metadata annotations for the cells
all_metadata<-read.csv("/srv/scratch/z3403160/lung_scRNAseq/2097-Lungcancer_metadata.csv")
all_metadata$CellType<-gsub("_", " ",all_metadata$CellType)

#now merge with data above
alltogetherv2<-merge(all_metadata, alltogether, by.x = "Cell", by.y = "CellID")

# define score function ---------------------------------------------------
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

# Apply Risk Signature ---------------------------------------------------------
#use the dataframe to generate the scores. 

#read in the risk signature
LUSC_lasso1<-read.csv("./LUSC_TvNT_DEGcoremat_lassomodelalphahalflambda1se_univariatecoef.csv", sep = ",", row.names = 1)


#create a submatrix with just these genes
alltogetherv2_LUSCrisk<-alltogetherv2[,which(colnames(alltogetherv2) %in% rownames(LUSC_lasso1))]
alltogetherv2_LUSCrisk<-alltogetherv2_LUSCrisk[,order(match(colnames(alltogetherv2_LUSCrisk), rownames(LUSC_lasso1[-1,])))]

#separate into positive OR and negative OR and apply the score to each separately
negORgenes<-rownames(LUSC_lasso1[which(LUSC_lasso1$Coef<0),])
posORgenes<-rownames(LUSC_lasso1[which(LUSC_lasso1$Coef>0),])

alltogetherv2_LUSCrisk_posOR<-alltogetherv2_LUSCrisk[,which(colnames(alltogetherv2_LUSCrisk) %in% posORgenes)]
alltogetherv2_LUSCrisk_negOR<-alltogetherv2_LUSCrisk[,which(colnames(alltogetherv2_LUSCrisk) %in% negORgenes)]

LUSC_lasso1_posOR<-LUSC_lasso1[which(rownames(LUSC_lasso1) %in% posORgenes),]
LUSC_lasso1_negOR<-LUSC_lasso1[which(rownames(LUSC_lasso1) %in% negORgenes),]
#use absolute value of coefficient for negOR so that large positive numbers reflect high scores for negative OR genes
LUSC_lasso1_negOR$Coef<--LUSC_lasso1_negOR$Coef

#apply risk score spearately to positive and negative OR genes
#posOR genes
alltogetherv2_LUSCriskscore_posOR<-LRscore(LUSC_lasso1_posOR, alltogetherv2_LUSCrisk_posOR)
alltogetherv2_LUSCriskscoreval_posOR<-alltogetherv2_LUSCriskscore_posOR[[2]]

#negOR genes
alltogetherv2_LUSCriskscore_negOR<-LRscore(LUSC_lasso1_negOR, alltogetherv2_LUSCrisk_negOR)
alltogetherv2_LUSCriskscoreval_negOR<-alltogetherv2_LUSCriskscore_negOR[[2]]

#append the scores to the original dataframe
alltogetherv3<-data.frame(alltogetherv2,
                          LUSCRiskScore_PosOR = alltogetherv2_LUSCriskscoreval_posOR,
                          LUSCRiskScore_NegOR = alltogetherv2_LUSCriskscoreval_negOR)
alltogetherv3$CellType<-as.factor(alltogetherv3$CellType)
levels(alltogetherv3$CellType)<-gsub("_", " ", levels(alltogetherv3$CellType))

levels(alltogetherv3$CellType)<-c("Alv.", "B Cell", "Cancer", "Endo.", "Epi.", "Eryth.", "Fib.", "Mast", "Mye.", "T Cell")


#plot violin plot for each cell type for pos and negative
scrna_LUSCrisk_posORbycelltype_vio<-ggplot(alltogetherv3, aes(x=as.factor(CellType), y=LUSCRiskScore_PosOR, fill = as.factor(CellType)))+
                                            theme_bw()+
                                            geom_violin(trim=TRUE, width = 1.8, size = 0.25) +
                                            xlab("")+
                                            ylab("Positive Expression Score")+
                                            theme(axis.text=element_text(size = 8, face="bold"),
                                                  axis.title=element_text(size=8, face="bold"),
                                                  plot.title=element_text(size = 12, face="bold", hjust= 0.5),
                                                  panel.grid.major = element_blank(),
                                                  panel.grid.minor= element_blank(),
                                                  legend.position= 'none')

scrna_LUSCrisk_negORbycelltype_vio<-ggplot(alltogetherv3, aes(x=as.factor(CellType), y=LUSCRiskScore_NegOR, fill = as.factor(CellType)))+
                                            theme_bw()+
                                            geom_violin(trim=TRUE, width = 1.8, size = 0.25) +
                                            xlab("")+
                                            ylab("Negative Expression Score")+
                                            theme(axis.text=element_text(size = 8, face="bold"),
                                                  axis.title=element_text(size=, face="bold"),
                                                  plot.title=element_text(size = 12, face="bold", hjust= 0.5),
                                                  panel.grid.major = element_blank(),
                                                  panel.grid.minor= element_blank(),
                                                  legend.position= 'none')


library("cowplot")
scrna_LUSCrisk_plots<-ggdraw() +
  draw_plot(scrna_LUSCrisk_posORbycelltype_vio, x = 0, y = 0.5, width = 1, height = 0.5)+
  draw_plot(scrna_LUSCrisk_negORbycelltype_vio, x = 0, y = 0, width = 1, height = 0.5)

scrna_LUSCrisk_posORbycelltype_vlnplot<-VlnPlot(all, features = "LUSCRiskScore_PosOR", pt.size = FALSE, stack = FALSE) +
  theme(legend.position='none')+
  ylab("Positive Score")+
  ggtitle(NULL)+
  theme(
    axis.text.x=element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()
  )

scrna_LUSCrisk_negORbycelltype_vlnplot<-VlnPlot(all, features = "LUSCRiskScore_NegOR", pt.size = FALSE, stack = FALSE) +
  theme(legend.position='none')+
  ylab("Negative Score")+
  ggtitle(NULL)

library("cowplot")
scrna_LUSCrisk_vlnplots<-ggdraw() +
  draw_plot(scrna_LUSCrisk_posORbycelltype_vlnplot, x = 0, y = 0.5, width = 1, height = 0.5)+
  draw_plot(scrna_LUSCrisk_negORbycelltype_vlnplot, x = 0, y = 0, width = 1, height = 0.5)

scrna_LUSCrisk_vlnplotsv2<- scrna_LUSCrisk_posORbycelltype_vlnplot +
  scrna_LUSCrisk_negORbycelltype_vlnplot +
  plot_layout(ncol = 1)

#plot for LUSC only and compare tumour adn non-tumour tissue
scrnaLUSConly_LUSCrisk_posORbycelltype_vlnplot<-VlnPlot(all_LUSC, features = "LUSCRiskScore_PosOR",
                                                        split.by="all_LUSC_TvNT",
                                                        split.plot=TRUE,
                                                        pt.size = FALSE, stack = FALSE) +
  #theme(legend.position='none')+
  ylab("Positive Score")+
  ggtitle(NULL)+
  theme(
    axis.text.x=element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()
  )

scrnaLUSConly_LUSCrisk_negORbycelltype_vlnplot<-VlnPlot(all_LUSC, features = "LUSCRiskScore_NegOR",
                                                        split.by="all_LUSC_TvNT",
                                                        split.plot=TRUE,
                                                        pt.size = FALSE, stack = FALSE) +
  #theme(legend.position='none')+
  ylab("Negative Score")+
  ggtitle(NULL)

scrnaLUSConly_LUSCrisk_vlnplotsTvNT<- scrnaLUSConly_LUSCrisk_posORbycelltype_vlnplot +
                                      scrnaLUSConly_LUSCrisk_negORbycelltype_vlnplot +
                                      plot_layout(ncol = 1)

#plot all cells as tsnes
all_tsne_LUSCriskpos<-FeaturePlot(all,
                                  features="LUSCRiskScore_PosOR", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n=11, name = "RdBu"))) +
  ggtitle("Positive Expression Score")
all_tsne_LUSCriskneg<-FeaturePlot(all,
                                  features="LUSCRiskScore_NegOR", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n=11, name = "RdBu"))) +
  ggtitle("Negative  Expression Score")
scrna_LUSCrisk_tsneplots<-all_tsneplot+
  all_tsne_LUSCriskpos +
  all_tsne_LUSCriskneg+
  plot_layout(ncol=3)

#plot LUSC cells only as tsnes
allLUSC_tsne_LUSCriskpos<-FeaturePlot(all_LUSC,
                                  features="LUSCRiskScore_PosOR", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n=11, name = "RdBu"))) +
  ggtitle("Positive Expression Score")
allLUSC_tsne_LUSCriskneg<-FeaturePlot(all_LUSC,
                                  features="LUSCRiskScore_NegOR", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n=11, name = "RdBu"))) +
  ggtitle("Negative  Expression Score")
#arrange using patchwork insted of cowplot
scrnaLUSConly_LUSCrisk_tsneplots<-allLUSC_tsneplot+
  allLUSC_tsne_LUSCriskpos +
  allLUSC_tsne_LUSCriskneg+
  plot_layout(ncol=3)

#see if you can extract the information of the scores to run a wilcox test for TvNT
all_LUSC_riskdfforstats<-data.frame(all_LUSC_TvNT = all_LUSC$all_LUSC_TvNT,
                                    CellType = all_LUSC$CellType,
                                    LUSCRiskScore_PosOR = all_LUSC$LUSCRiskScore_PosOR,
                                    LUSCRiskScore_NegOR = all_LUSC$LUSCRiskScore_NegOR,
                                    LUSCProgScore_PosHR = all_LUSC$LUSCProgScore_PosHR,
                                    LUSCProgScore_NegHR = all_LUSC$LUSCProgScore_NegHR,
                                    CorrClust1_1 = all_LUSC$CorrClust1_1,
                                    CorrClust2_1 = all_LUSC$CorrClust2_1,
                                    CorrClust3_1 = all_LUSC$CorrClust3_1,
                                    CorrClust4_1 = all_LUSC$CorrClust4_1,
                                    Intergrins1 = all_LUSC$Intergrins1,
                                    Integrins_sigLigRec_1 = all_LUSC$Integrins_sigLigRec_1,
                                    MitogenRecs_1 = all_LUSC$MitogenRecs_1,
                                    MitogenLigs_1 = all_LUSC$MitogenLigs_1)

#wilcox test for different cell types in risk score
#fibroblasts
testingdf<-all_LUSC_riskdfforstats[which(all_LUSC_riskdfforstats$CellType == "Fibroblast"),]
all_LUSC_posriskscore<-wilcox.test(testingdf$LUSCRiskScore_PosOR~testingdf$all_LUSC_TvNT)

all_LUSC_negriskscore<-wilcox.test(testingdf$LUSCRiskScore_NegOR~testingdf$all_LUSC_TvNT)

#endothelial
testingdf<-all_LUSC_riskdfforstats[which(all_LUSC_riskdfforstats$CellType == "EC"),]
all_LUSC_negriskscore<-wilcox.test(testingdf$LUSCRiskScore_NegOR~testingdf$all_LUSC_TvNT)

all_LUSC_posriskscore<-wilcox.test(testingdf$LUSCRiskScore_PosOR~testingdf$all_LUSC_TvNT)

#alveolar
testingdf<-all_LUSC_riskdfforstats[which(all_LUSC_riskdfforstats$CellType == "Alveolar"),]
all_LUSC_negriskscore<-wilcox.test(testingdf$LUSCRiskScore_NegOR~testingdf$all_LUSC_TvNT)

all_LUSC_posriskscore<-wilcox.test(testingdf$LUSCRiskScore_PosOR~testingdf$all_LUSC_TvNT)

#myeloid
testingdf<-all_LUSC_riskdfforstats[which(all_LUSC_riskdfforstats$CellType == "Myeloid"),]
all_LUSC_posriskscore<-wilcox.test(testingdf$LUSCRiskScore_PosOR~testingdf$all_LUSC_TvNT)
all_LUSC_negriskscore<-wilcox.test(testingdf$LUSCRiskScore_NegOR~testingdf$all_LUSC_TvNT)

#Tcell
testingdf<-all_LUSC_riskdfforstats[which(all_LUSC_riskdfforstats$CellType == "T cell"),]
all_LUSC_posriskscore<-wilcox.test(testingdf$LUSCRiskScore_PosOR~testingdf$all_LUSC_TvNT)
all_LUSC_negriskscore<-wilcox.test(testingdf$LUSCRiskScore_NegOR~testingdf$all_LUSC_TvNT)

#mast cell
testingdf<-all_LUSC_riskdfforstats[which(all_LUSC_riskdfforstats$CellType == "Mast cell"),]
all_LUSC_posriskscore<-wilcox.test(testingdf$LUSCRiskScore_PosOR~testingdf$all_LUSC_TvNT)
all_LUSC_negriskscore<-wilcox.test(testingdf$LUSCRiskScore_NegOR~testingdf$all_LUSC_TvNT)

#B cell
testingdf<-all_LUSC_riskdfforstats[which(all_LUSC_riskdfforstats$CellType == "B cell"),]
all_LUSC_posriskscore<-wilcox.test(testingdf$LUSCRiskScore_PosOR~testingdf$all_LUSC_TvNT)
#B cell TvNT for pos risk score:
all_LUSC_negriskscore<-wilcox.test(testingdf$LUSCRiskScore_NegOR~testingdf$all_LUSC_TvNT)

#Cancer
testingdf<-all_LUSC_riskdfforstats[which(all_LUSC_riskdfforstats$CellType == "Cancer"),]
all_LUSC_posriskscore<-wilcox.test(testingdf$LUSCRiskScore_PosOR~testingdf$all_LUSC_TvNT)
all_LUSC_negriskscore<-wilcox.test(testingdf$LUSCRiskScore_NegOR~testingdf$all_LUSC_TvNT)

#wilcox test for different cell types in correlation clusters
#fibroblasts
testingdf<-all_LUSC_riskdfforstats[which(all_LUSC_riskdfforstats$CellType == "Fibroblast"),]
all_LUSC_corrclust1_score<-wilcox.test(testingdf$CorrClust1_1~testingdf$all_LUSC_TvNT) 
all_LUSC_corrclust2_score<-wilcox.test(testingdf$CorrClust2_1~testingdf$all_LUSC_TvNT) 
all_LUSC_corrclust3_score<-wilcox.test(testingdf$CorrClust3_1~testingdf$all_LUSC_TvNT) 
all_LUSC_corrclust4_score<-wilcox.test(testingdf$CorrClust4_1~testingdf$all_LUSC_TvNT) 

#endothelial cells
testingdf<-all_LUSC_riskdfforstats[which(all_LUSC_riskdfforstats$CellType == "EC"),]
all_LUSC_corrclust1_score<-wilcox.test(testingdf$CorrClust1_1~testingdf$all_LUSC_TvNT) 
all_LUSC_corrclust2_score<-wilcox.test(testingdf$CorrClust2_1~testingdf$all_LUSC_TvNT) 
all_LUSC_corrclust3_score<-wilcox.test(testingdf$CorrClust3_1~testingdf$all_LUSC_TvNT) 
all_LUSC_corrclust4_score<-wilcox.test(testingdf$CorrClust4_1~testingdf$all_LUSC_TvNT) 

#epithelial cells are only Non-tumor so can't compare tumor vs non-tumor
testingdf<-all_LUSC_riskdfforstats[which(all_LUSC_riskdfforstats$CellType == "Epithelial"),]

#myeloid cells
testingdf<-all_LUSC_riskdfforstats[which(all_LUSC_riskdfforstats$CellType == "Myeloid"),]
all_LUSC_corrclust1_score<-wilcox.test(testingdf$CorrClust1_1~testingdf$all_LUSC_TvNT) 
all_LUSC_corrclust2_score<-wilcox.test(testingdf$CorrClust2_1~testingdf$all_LUSC_TvNT) 
all_LUSC_corrclust3_score<-wilcox.test(testingdf$CorrClust3_1~testingdf$all_LUSC_TvNT) 
all_LUSC_corrclust4_score<-wilcox.test(testingdf$CorrClust4_1~testingdf$all_LUSC_TvNT)

#Cancer
testingdf<-all_LUSC_riskdfforstats[which(all_LUSC_riskdfforstats$CellType == "Cancer"),]
all_LUSC_corrclust1_score<-wilcox.test(testingdf$CorrClust1_1~testingdf$all_LUSC_TvNT) 
all_LUSC_corrclust2_score<-wilcox.test(testingdf$CorrClust2_1~testingdf$all_LUSC_TvNT) 
all_LUSC_corrclust3_score<-wilcox.test(testingdf$CorrClust3_1~testingdf$all_LUSC_TvNT) 
all_LUSC_corrclust4_score<-wilcox.test(testingdf$CorrClust4_1~testingdf$all_LUSC_TvNT) 

#T-Cell
testingdf<-all_LUSC_riskdfforstats[which(all_LUSC_riskdfforstats$CellType == "T cell"),]
all_LUSC_corrclust1_score<-wilcox.test(testingdf$CorrClust1_1~testingdf$all_LUSC_TvNT) 
all_LUSC_corrclust2_score<-wilcox.test(testingdf$CorrClust2_1~testingdf$all_LUSC_TvNT) 
all_LUSC_corrclust3_score<-wilcox.test(testingdf$CorrClust3_1~testingdf$all_LUSC_TvNT) 
all_LUSC_corrclust4_score<-wilcox.test(testingdf$CorrClust4_1~testingdf$all_LUSC_TvNT) 

#Mast Cell
testingdf<-all_LUSC_riskdfforstats[which(all_LUSC_riskdfforstats$CellType == "Mast cell"),]
all_LUSC_corrclust1_score<-wilcox.test(testingdf$CorrClust1_1~testingdf$all_LUSC_TvNT) 
all_LUSC_corrclust2_score<-wilcox.test(testingdf$CorrClust2_1~testingdf$all_LUSC_TvNT) 
all_LUSC_corrclust3_score<-wilcox.test(testingdf$CorrClust3_1~testingdf$all_LUSC_TvNT) 
all_LUSC_corrclust4_score<-wilcox.test(testingdf$CorrClust4_1~testingdf$all_LUSC_TvNT) 

#B Cell
testingdf<-all_LUSC_riskdfforstats[which(all_LUSC_riskdfforstats$CellType == "B cell"),]
all_LUSC_corrclust1_score<-wilcox.test(testingdf$CorrClust1_1~testingdf$all_LUSC_TvNT) 
all_LUSC_corrclust2_score<-wilcox.test(testingdf$CorrClust2_1~testingdf$all_LUSC_TvNT) 
all_LUSC_corrclust3_score<-wilcox.test(testingdf$CorrClust3_1~testingdf$all_LUSC_TvNT) 
all_LUSC_corrclust4_score<-wilcox.test(testingdf$CorrClust4_1~testingdf$all_LUSC_TvNT) 

#Alveolar
testingdf<-all_LUSC_riskdfforstats[which(all_LUSC_riskdfforstats$CellType == "Alveolar"),]
all_LUSC_corrclust1_score<-wilcox.test(testingdf$CorrClust1_1~testingdf$all_LUSC_TvNT) 
all_LUSC_corrclust2_score<-wilcox.test(testingdf$CorrClust2_1~testingdf$all_LUSC_TvNT) 
all_LUSC_corrclust3_score<-wilcox.test(testingdf$CorrClust3_1~testingdf$all_LUSC_TvNT) 
all_LUSC_corrclust4_score<-wilcox.test(testingdf$CorrClust4_1~testingdf$all_LUSC_TvNT) 

#Erythroblast
#are T only in SqCC it seems so can't run the TvNT comparison 
testingdf<-all_LUSC_riskdfforstats[which(all_LUSC_riskdfforstats$CellType == "Erythroblast"),]

# Matrisome expressed by normal fibroblasts --------------------

#subset the LUSC matrix to just include the fibroblasts 
all_LUSC_fib<-subset(all_LUSC, idents = "Fibroblast")

#reassign fibroblast identities based on tumour or non-tumour derivation
Idents(all_LUSC_fib)<-all_LUSC_fib@meta.data$all_LUSC_TvNT

all_LUSC_fib_TvNTmarkers<-FindMarkers(all_LUSC_fib, ident.1 = "T", ident.2 = "N")

#subset to include only matrisomal genes
matrisome_master<-read.csv(file = './matrisome_hs_masterlist.csv', header = TRUE)
matrisome_master_min<-matrisome_master[1:4]

#filter by matrisome
all_LUSC_fib_TvNTmarkers_mat<-all_LUSC_fib_TvNTmarkers[which(rownames(all_LUSC_fib_TvNTmarkers) %in% matrisome_master_min$Gene.Symbol),]

#filter by core matrisome
all_LUSC_fib_TvNTmarkers_coremat<-all_LUSC_fib_TvNTmarkers[which(rownames(all_LUSC_fib_TvNTmarkers) %in% matrisome_master_min[which(matrisome_master_min$Division=="Core matrisome"),]$Gene.Symbol),]

#top 5 upregulated genes: COL3A1, COL1A1, CTHRC1, COL1A2, POSTN
#top 5 downregulated genes: SFTPC, SPARCL1, PRELP, A2M, CYR61
allLUSC_fib_TvNTup_vio<-VlnPlot(all_LUSC_fib, features = c("COL3A1", "COL1A1", "CTHRC1", "COL1A2", "POSTN"),
                                #split.plot=TRUE,
                                #split.by="all_LUSC_TvNT",
                                pt.size = FALSE) #+ #theme(legend.position='none')+

allLUSC_fib_TvNTdn_vio<-VlnPlot(all_LUSC_fib, features = c("SFTPC", "SPARCL1", "PRELP", "A2M", "CYR61"),
                                  #split.plot=TRUE,
                                  #split.by="all_LUSC_TvNT",
                                  pt.size = FALSE)# + #theme(legend.position='none')+
 
allLUSC_fib_TvNTupanddn_vio<-VlnPlot(all_LUSC_fib, features = c("COL3A1", "COL1A1", "CTHRC1", "COL1A2", "POSTN",
                                                           "SFTPC", "SPARCL1", "PRELP", "A2M", "CYR61"),
                                ncol = 5,
                                pt.size = FALSE)

#top 5 up and top 5 down from core matrisome only
allLUSC_fib_TvNTupanddn_coremat_vio<-VlnPlot(all_LUSC_fib, features = c("COL3A1", "COL1A1", "CTHRC1", "COL1A2", "POSTN",
                                                                "SPARCL1", "PRELP", "CYR61", "CTGF", "TNXB"),
                                     ncol = 5,
                                     pt.size = FALSE)

# Add scores to metadata in seurat object. --------------------------------
alltogetherv3<-alltogetherv2
colnames(alltogetherv3[,25671:25674]) #contains the score values

#check cell order
rownames(all@assays$RNA) #is the gene names
colnames(all@assays$RNA) #is the cell barcodes

colnames(all@assays$RNA)[1:10]
alltogetherv3$Cell[1:10] #these look like they are in the same order but check to make sure

alltogetherv3<-alltogetherv3[order(match(alltogetherv3$Cell, colnames(all@assays$RNA))),]
#this looks like it's in the same order as the data in the object

#add score values as metadata
all@meta.data<-cbind(all@meta.data, alltogetherv3[,25671:25674])

#add cell assignments from the original manuscript
#all_metadata is the info from the paper
all_metadata_min<-all_metadata[which(all_metadata$Cell %in% colnames(all@assays$RNA)),] #because contains details for all cells
all_metadata_min<-all_metadata_min[order(match(all_metadata_min$Cell, colnames(all@assays$RNA))),]

all@meta.data<-cbind(all@meta.data, all_metadata_min)

#assign cell Idents as the original manuscript data
Idents(all)<-'CellType'

# Check cell type expression for correlation clusters ---------------------
LUSC_corrorder_names<-read.csv(file= "./TCGA_LUSC_CoreMatCorrClusters.csv", row.names = 1)

#generate each as a gene list and perform AddModuleScore to calculate enrichment.
LUSC_clust1<-LUSC_corrorder_names[which(LUSC_corrorder_names$correlationcluster==1),]
LUSC_clust2<-LUSC_corrorder_names[which(LUSC_corrorder_names$correlationcluster==2),]
LUSC_clust3<-LUSC_corrorder_names[which(LUSC_corrorder_names$correlationcluster==3),]
LUSC_clust4<-LUSC_corrorder_names[which(LUSC_corrorder_names$correlationcluster==4),]

#generate an enrichment score for 
all<-AddModuleScore(object = all,
                    features = list(LUSC_clust1$genes),
                    name = "CorrClust1_")
all<-AddModuleScore(object = all,
                    features = list(LUSC_clust2$genes),
                    name = "CorrClust2_")
all<-AddModuleScore(object = all,
                    features = list(LUSC_clust3$genes),
                    name = "CorrClust3_")
all<-AddModuleScore(object = all,
                    features = list(LUSC_clust4$genes),
                    name = "CorrClust4_")

library(RColorBrewer)
library(tidyverse)

#Plot TSNE by clusters
all_tsneplot<-DimPlot(all, reduction = 'tsne', label = TRUE) #+NoLegend()

all_tsne_corrclust1<-FeaturePlot(all,
                                features="CorrClust1_1", label = TRUE, repel = TRUE) +
                                scale_colour_gradientn(colours = rev(brewer.pal(n=11, name = "RdBu"))) +
                                ggtitle("Cluster 1")
all_vio_corrclust1<-VlnPlot(all, features = "CorrClust1_1", pt.size = FALSE) +
                            theme(legend.position='none')+
                            ylab("Cluster 1")+
  ggtitle(NULL)+
  theme(
    axis.text.x=element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()
  )

all_tsne_corrclust2<-FeaturePlot(all,
                            features="CorrClust2_1", label = TRUE, repel = TRUE) +
                            scale_colour_gradientn(colours = rev(brewer.pal(n=11, name = "RdBu"))) +
                            ggtitle("Cluster 2")
all_vio_corrclust2<-VlnPlot(all, features = "CorrClust2_1", pt.size = FALSE) +
                            theme(legend.position='none')+
                            ylab("Cluster 2")+
  ggtitle(NULL)+
  theme(
    axis.text.x=element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()
  )

all_tsne_corrclust3<- FeaturePlot(all,
                                  features="CorrClust3_1", label = TRUE, repel = TRUE) +
                                  scale_colour_gradientn(colours = rev(brewer.pal(n=11, name = "RdBu"))) +
                                  ggtitle("Cluster 3")
all_vio_corrclust3<-VlnPlot(all, features = "CorrClust3_1", pt.size = FALSE) +
                            theme(legend.position='none')+
                            ylab("Cluster 3")+
  ggtitle(NULL)+
  theme(
    axis.text.x=element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()
  )

all_tsne_corrclust4<-FeaturePlot(all,
                                features="CorrClust4_1", label = TRUE, repel = TRUE) +
                                scale_colour_gradientn(colours = rev(brewer.pal(n=11, name = "RdBu"))) +
                                ggtitle("Cluster 4")
all_vio_corrclust4<-VlnPlot(all, features = "CorrClust4_1", pt.size = FALSE) +
                    theme(legend.position='none')+
                    ylab("Cluster 4")+
                    ggtitle(NULL)

library(patchwork)

all_tsne_corr_arrange<- all_tsneplot + {
                            all_tsne_corrclust1 +
                              all_tsne_corrclust2 +
                              plot_layout(ncol=1)
                          } + {
                            all_tsne_corrclust3 +
                            all_tsne_corrclust4 +
                            plot_layout(ncol=1)
                          }+
  plot_layout(widths = c(2,1,1))

all_celltypevio_corr_arrange<- all_vio_corrclust1+
                                all_vio_corrclust2+
                                all_vio_corrclust3+
                                all_vio_corrclust4+
                                plot_layout(ncol = 1)

#run plots just for LUSC cells 
#Plot TSNE by clusters
library(viridis)
allLUSC_tsne_corrclust1<-FeaturePlot(all_LUSC,
                                 features="CorrClust1_1", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n=11, name = "RdBu"))) +
  ggtitle("Cluster 1")
allLUSC_vio_corrclust1<-VlnPlot(all_LUSC, features = "CorrClust1_1",
                                split.plot=TRUE,
                                split.by="all_LUSC_TvNT",
                                pt.size = FALSE) + 
  ylab("Cluster 1")+
  ggtitle(NULL)+
  theme(
    axis.text.x=element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()
  )


allLUSC_tsne_corrclust2<-FeaturePlot(all_LUSC,
                                 features="CorrClust2_1", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n=11, name = "RdBu"))) +
  ggtitle("Cluster 2")
allLUSC_vio_corrclust2<-VlnPlot(all_LUSC, features = "CorrClust2_1",
                                split.plot=TRUE,
                                split.by="all_LUSC_TvNT",
                                pt.size = FALSE) +
  ylab("Cluster 2")+
  ggtitle(NULL)+
  theme(
    axis.text.x=element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()
  )

allLUSC_tsne_corrclust3<- FeaturePlot(all_LUSC,
                                  features="CorrClust3_1", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n=11, name = "RdBu"))) +
  ggtitle("Cluster 3")
allLUSC_vio_corrclust3<-VlnPlot(all_LUSC, features = "CorrClust3_1",
                                split.plot = TRUE,
                                split.by="all_LUSC_TvNT",
                                pt.size = FALSE) +
  ylab("Cluster 3")+
  ggtitle(NULL)+
  theme(
    axis.text.x=element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()
  )

allLUSC_tsne_corrclust4<-FeaturePlot(all_LUSC,
                                 features="CorrClust4_1", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n=11, name = "RdBu"))) +
  ggtitle("Cluster 4")
allLUSC_vio_corrclust4<-VlnPlot(all_LUSC, features = "CorrClust4_1",
                                split.plot = TRUE,
                                split.by="all_LUSC_TvNT",
                                pt.size = FALSE) +
  ylab("Cluster 4")+
  ggtitle(NULL)

library(patchwork)

allLUSC_tsne_corr_arrange<- allLUSC_tsneplot + {
  allLUSC_tsne_corrclust1 +
    allLUSC_tsne_corrclust2 +
    plot_layout(ncol=1)
} + {
  allLUSC_tsne_corrclust3 +
    allLUSC_tsne_corrclust4 +
    plot_layout(ncol=1)
}+
  plot_layout(widths = c(2,1,1))

allLUSC_celltypevio_corr_arrange<- allLUSC_vio_corrclust1+
  allLUSC_vio_corrclust2+
  allLUSC_vio_corrclust3+
  allLUSC_vio_corrclust4+
  plot_layout(ncol = 1, guides="collect")