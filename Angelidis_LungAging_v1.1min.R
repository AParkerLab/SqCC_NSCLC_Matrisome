#Lung Aging Analysis

#Angelidis, I. et al. An atlas of the aging lung mapped by single cell transcriptomics and deep tissue proteomics. Nat. Commun. 10, 963 (2019).

# Angelidis load in data --------------------------------------------------
#data provided as normalised counts

Angelidis_rna<-read.table(file = "./GSE124872_raw_counts_whole_lung_bulk.txt", sep = "\t", header = T)
#this is a count matrix for the RNA

# Angelidis read in risk and prognosis score signatures ------------------------------
LUSC_lasso2<-read.csv(file = "./LUSC_TvNT_DEGcoremat_lassomodelalphahalflambda1se_univariatecoef.csv", sep = ",", header = T, row.names = 1)

#define the function to apply the score
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

# convert RNAseq ENSG to HUGO symbols --------------------------------------------
ensg_genes<-Angelidis_rna$Geneid

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")

require(biomaRt)
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

listAttributes(human)
listAttributes(mouse)

ensg_musgenestohshugo = getLDS(attributes = c("ensembl_gene_id"),
                            #attributesL = "hsapiens_homolog_ensembl_gene",
                            filters = "ensembl_gene_id",
                            values = ensg_genes ,
                            mart = mouse,
                            attributesL = c("hgnc_symbol","entrezgene_id", "ensembl_gene_id"),
                            martL = human,
                            uniqueRows=T)#,
                            #useCache = FALSE)


Angelidis_hugo<-merge(ensg_musgenestohshugo, Angelidis_rna, by.x = "Gene.stable.ID", by.y = "Geneid")
Angelidis_hugomin<-Angelidis_hugo[,-c(1,3:9)]
Angelidis_hugomint<-t(Angelidis_hugomin[,-1])
colnames(Angelidis_hugomint)<-Angelidis_hugomin$HGNC.symbol

library("edgeR")
#store the raw counts matrices in an object called DGElist
Angelidis_DGElist<-DGEList(counts = t(Angelidis_hugomint), group = c("old","old", "old", "young", "young", "young"))

#Filter out genes with less than 10 counts i.e. keep genes with more than 10 counts in at least two samples
#library sizes are approximately 60 million so cpm=0.6 corresponds to 10 counts per sample
Angelidis_keepfilterbyexpr<-filterByExpr(Angelidis_DGElist)

Angelidis_DGElistf<-Angelidis_DGElist[Angelidis_keepfilterbyexpr, , keep.lib.sizes = FALSE]

#Normalize using the TMM method which will calculate nromalization factors for each
Angelidis_DGElistfTMM<-calcNormFactors(Angelidis_DGElistf)

#use edgeR to convert the counts to logCPM
library(limma)
Angelidis_logcpm<-cpm(Angelidis_DGElistfTMM, log = TRUE, prior.count = 3) #prior count is used to damp th variances of the logs of low counts


#z transform the expression data
Angelidis_logcpmz<-apply(Angelidis_logcpm,1,scale)
rownames(Angelidis_logcpmz)<-colnames(Angelidis_logcpm)

# Angelidis bulk RNASeq Apply risk score to all samples ----------------
Angelidis_logcpmz_risk<-Angelidis_logcpmz[, which(colnames(Angelidis_logcpmz) %in% rownames(LUSC_lasso2))]

LUSC_lasso2_Angelidis<-LUSC_lasso2[which(rownames(LUSC_lasso2) %in% colnames(Angelidis_logcpmz)),]

Angelidis_logcpmz_risk<-Angelidis_logcpmz_risk[, order(match(colnames(Angelidis_logcpmz_risk), rownames(LUSC_lasso2_Angelidis)))]

Angelidis_riskscore<-LRscore(LUSC_lasso2_Angelidis, Angelidis_logcpmz_risk)

Angelidis_riskscoredf<-data.frame(Accs = rownames(Angelidis_logcpmz_risk),
                                  RiskScore = Angelidis_riskscore[[2]])

