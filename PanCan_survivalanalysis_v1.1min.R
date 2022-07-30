
#survival analysis for PanCan ECM Matreotypes
  
colnames(pancan_rnaseqlTminz_survsig_scoredf3)[36]<-"Score"
pancan_rnaseqlTminz_survsig_scoredf3$Score<-as.factor(pancan_rnaseqlTminz_survsig_scoredf3$Score)

library(survminer)
library(survival)
#now run survival analysis on this pancancer dataset
pancan_rnaseqlTminz_survsig_scoredf4<-pancan_rnaseqlTminz_survsig_scoredf3[-which(pancan_rnaseqlTminz_survsig_scoredf3$type =="LUSC"),]
tumourtypes<-levels(as.factor(pancan_rnaseqlTminz_survsig_scoredf4$type))

tumourrows<-list()
tumouraccs<-list()
for (i in 1:length(tumourtypes)){
  tumourrows[[i]]<-which(pancan_rnaseqlTminz_survsig_scoredf4$type==tumourtypes[i])
  tumouraccs[[i]]<-pancan_rnaseqlTminz_survsig_scoredf4$Accs[which(pancan_rnaseqlTminz_survsig_scoredf4$type==tumourtypes[i])]
}
names(tumourrows)<-tumourtypes
names(tumouraccs)<-tumourtypes


# Molecular Subtype and Survival - Euclidean ------------------------------------------
#overall survival
  p<-list()
  p.val<-c()
  #some tumours have two or three molecular subtypes assigned
  colourpalette<-c("K1" = "royalblue","K2" ="gold","K3" = "forestgreen")
  
  for (i in 1:length(tumourrows)){
    iformula <- as.formula(sprintf("Surv(OS.time, OS) ~ as.factor(Distmethod_cluster)"))
    pancan_progscoremed_fit<-survfit(iformula,
                                     data = pancan_rnaseqlTminz_survsig_scoredf4[tumourrows[[i]],],
                                     na.action = na.omit)
    
    pvalue<-survdiff(iformula,
                     data = pancan_rnaseqlTminz_survsig_scoredf4[tumourrows[[i]],],
                     na.action = na.omit)
    p.val[i] <- format(round((1 - pchisq(pvalue$chisq, length(pvalue$n) - 1)),4))
    
    pancan_progscoremed_fit$call$formula <- eval(pancan_progscoremed_fit$call$formula) 
    
    #define a specific colour palette for this plot
    clusterspresent<-droplevels(pancan_rnaseqlTminz_survsig_scoredf4$Distmethod_cluster[tumourrows[[i]]])
    
    tumourcolourpalette<-colourpalette[names(colourpalette) %in% levels(clusterspresent)]
    
    p[[i]] <- ggsurvplot_list(pancan_progscoremed_fit,
                              data = pancan_rnaseqlTminz_survsig_scoredf4,
                              pval = p.val[i],
                              title = names(tumourrows)[i],
                              risk.table = FALSE,
                              xlab = "Days Since Diagnosis",
                              ylab = "Survival Probability",
                              palette = unname(tumourcolourpalette))
  }
  
  
  #repeat with Disease-Specific survival
 
  p<-list()
  p.val<-c()
  for (i in 1:length(tumourrows)){
    iformula <- as.formula(sprintf("Surv(DSS.time, DSS) ~ as.factor(Distmethod_cluster)"))
    pancan_progscoremed_fit<-survfit(iformula,
                                     data = pancan_rnaseqlTminz_survsig_scoredf4[tumourrows[[i]],],
                                     na.action = na.omit)
    
    pvalue<-survdiff(iformula,
                     data = pancan_rnaseqlTminz_survsig_scoredf4[tumourrows[[i]],],
                     na.action = na.omit)
    p.val[i] <- format(round((1 - pchisq(pvalue$chisq, length(pvalue$n) - 1)),4))
    
    pancan_progscoremed_fit$call$formula <- eval(pancan_progscoremed_fit$call$formula) 
    
    #define a specific colour palette for this plot
    clusterspresent<-droplevels(pancan_rnaseqlTminz_survsig_scoredf4$Distmethod_cluster[tumourrows[[i]]])
    
    tumourcolourpalette<-colourpalette[names(colourpalette) %in% levels(clusterspresent)]
    
    p[[i]] <- ggsurvplot_list(pancan_progscoremed_fit,
                              data = pancan_rnaseqlTminz_survsig_scoredf4,
                              pval = p.val[i],
                              title = names(tumourrows)[i],
                              risk.table = FALSE,
                              xlab = "Days Since Diagnosis",
                              ylab = "Survival Probability",
                              palette = unname(tumourcolourpalette))
  }
  
  

  #Disease-free survival
  #redefine tumourrows for new dataframe
  
  p<-list()
  p.val<-c()
  for (i in 1:length(tumourrowsDFI)){
    iformula <- as.formula(sprintf("Surv(DFI.time, DFI) ~ as.factor(Distmethod_cluster)"))
    pancan_progscoremed_fit<-survfit(iformula,
                                     data = pancan_rnaseqlTminz_survsig_scoredf4DFI[tumourrowsDFI[[i]],],
                                     na.action = na.omit)
    
    pvalue<-survdiff(iformula,
                     data = pancan_rnaseqlTminz_survsig_scoredf4DFI[tumourrowsDFI[[i]],],
                     na.action = na.omit)
    p.val[i] <- format(round((1 - pchisq(pvalue$chisq, length(pvalue$n) - 1)),4))
    
    pancan_progscoremed_fit$call$formula <- eval(pancan_progscoremed_fit$call$formula) 
    
    #define a specific colour palette for this plot
    clusterspresent<-droplevels(pancan_rnaseqlTminz_survsig_scoredf4DFI$Distmethod_cluster[tumourrowsDFI[[i]]])
    
    tumourcolourpalette<-colourpalette[names(colourpalette) %in% levels(clusterspresent)]
    tumourcolourpalette<-unname(tumourcolourpalette)
    
    p[[i]] <- ggsurvplot_list(pancan_progscoremed_fit,
                              data = pancan_rnaseqlTminz_survsig_scoredf4DFI,
                              pval = p.val[i],
                              title = names(tumourrowsDFI)[i],
                              risk.table = FALSE,
                              xlab = "Days Since Diagnosis",
                              ylab = "Survival Probability",#,
                              palette = tumourcolourpalette)#,
  }
  
# Mol subtype HR ----------------------------------------------------------
  #overall survival
  pancan_molsubtype_coxph<-list()
  p.value<-matrix(NA, ncol =1, nrow = length(tumourtypes))
  wald.test<-matrix(NA, ncol =1, nrow = length(tumourtypes))
  beta<-matrix(NA, ncol =1, nrow = length(tumourtypes))
  HR<-matrix(NA, ncol =1, nrow = length(tumourtypes))
  HR.confint.lower<-matrix(NA, ncol =1, nrow = length(tumourtypes))
  HR.confint.upper<-matrix(NA, ncol =1, nrow = length(tumourtypes))
  
  for (i in 1:length(tumourrows)){
    
    pancan_molsubtype_coxph[[i]]<-coxph(Surv(OS.time, OS) ~ Distmethod_cluster , data = pancan_rnaseqlTminz_survsig_scoredf4[tumourrows[[i]],])
    p.value[i]<-summary(pancan_molsubtype_coxph[[i]])$wald["pvalue"]
    wald.test[i]<-summary(pancan_molsubtype_coxph[[i]])$wald["test"]
    beta[i]<-summary(pancan_molsubtype_coxph[[i]])$coef[nrow(summary(pancan_molsubtype_coxph[[i]])$conf.int),1] #to take K3 only#coeficient beta
    HR[i] <-summary(pancan_molsubtype_coxph[[i]])$coef[nrow(summary(pancan_molsubtype_coxph[[i]])$conf.int),2] #to take K3 only #HR
    HR.confint.lower[i] <- summary(pancan_molsubtype_coxph[[i]])$conf.int[nrow(summary(pancan_molsubtype_coxph[[i]])$conf.int),"lower .95"] #to take K3 only
    HR.confint.upper[i] <- summary(pancan_molsubtype_coxph[[i]])$conf.int[nrow(summary(pancan_molsubtype_coxph[[i]])$conf.int),"upper .95"] #to take K3 only
  }
  pancan_molsubtype_OS_coxph<-data.frame(Tumourtype = tumourtypes,
                                     pval = p.value,
                                     waldtest = wald.test,
                                     beta = beta,
                                     HR = HR,
                                     HR.confint.lower = HR.confint.lower,
                                     HR.confint.upper = HR.confint.upper)
  rownames(pancan_molsubtype_OS_coxph) <-tumourtypes
  
  #Plot the HR
  pancan_molsubtype_OS_coxph_HR <- ggplot(pancan_molsubtype_OS_coxph, aes(x = HR, y = Tumourtype)) + 
    geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") + 
    geom_errorbarh(aes(xmax = HR.confint.upper, xmin = HR.confint.lower), size = .5, height = 
                     .2, color = "gray50") +
    geom_point(size = 3.5, color ="red2") +
    theme_bw()+
    theme(panel.grid.minor = element_blank()) +
    ylab("") +
    xlab("Hazard Ratio") +
    ggtitle("K3 Molecular Subtype Overall Survival")
  pancan_molsubtype_OS_coxph_HR
  
  
  #disease-specific survival
  pancan_molsubtype_coxph<-list()
  p.value<-matrix(NA, ncol =1, nrow = length(tumourtypes))
  wald.test<-matrix(NA, ncol =1, nrow = length(tumourtypes))
  beta<-matrix(NA, ncol =1, nrow = length(tumourtypes))
  HR<-matrix(NA, ncol =1, nrow = length(tumourtypes))
  HR.confint.lower<-matrix(NA, ncol =1, nrow = length(tumourtypes))
  HR.confint.upper<-matrix(NA, ncol =1, nrow = length(tumourtypes))
  
  for (i in 1:length(tumourrows)){
    
    pancan_molsubtype_coxph[[i]]<-coxph(Surv(DSS.time, DSS) ~ Distmethod_cluster , data = pancan_rnaseqlTminz_survsig_scoredf4[tumourrows[[i]],])
    p.value[i]<-summary(pancan_molsubtype_coxph[[i]])$wald["pvalue"]
    wald.test[i]<-summary(pancan_molsubtype_coxph[[i]])$wald["test"]
    beta[i]<-summary(pancan_molsubtype_coxph[[i]])$coef[nrow(summary(pancan_molsubtype_coxph[[i]])$conf.int),1] #to take K3 only#coeficient beta
    HR[i] <-summary(pancan_molsubtype_coxph[[i]])$coef[nrow(summary(pancan_molsubtype_coxph[[i]])$conf.int),2] #to take K3 only #HR
    HR.confint.lower[i] <- summary(pancan_molsubtype_coxph[[i]])$conf.int[nrow(summary(pancan_molsubtype_coxph[[i]])$conf.int),"lower .95"] #to take K3 only
    HR.confint.upper[i] <- summary(pancan_molsubtype_coxph[[i]])$conf.int[nrow(summary(pancan_molsubtype_coxph[[i]])$conf.int),"upper .95"] #to take K3 only
  }
  pancan_molsubtype_DSS_coxph<-data.frame(Tumourtype = tumourtypes,
                                         pval = p.value,
                                         waldtest = wald.test,
                                         beta = beta,
                                         HR = HR,
                                         HR.confint.lower = HR.confint.lower,
                                         HR.confint.upper = HR.confint.upper)
  rownames(pancan_molsubtype_DSS_coxph) <-tumourtypes
  write.csv(pancan_molsubtype_DSS_coxph, file = "./PanCan_LUSCMolSubtypeCox_DSS_211213.csv", sep = ",")
  
  #Plot the HR
  pancan_molsubtype_DSS_coxphv2<-pancan_molsubtype_DSS_coxph[order(pancan_molsubtype_DSS_coxph$HR, decreasing=TRUE),]
  
  pancan_molsubtype_DSS_coxphv2<-pancan_molsubtype_DSS_coxphv2[-6,]
  
  #reorder levels of tumour type 
  pancan_molsubtype_DSS_coxphv2$Tumourtype<-reorder(pancan_molsubtype_DSS_coxphv2$Tumourtype, pancan_molsubtype_DSS_coxphv2$HR)
  pancan_molsubtype_DSS_coxphv2$sigcolours<-ifelse(pancan_molsubtype_DSS_coxphv2$pval<0.05, "red2", "royalblue")
  
  pancan_molsubtype_DSS_coxph_HR <- ggplot(pancan_molsubtype_DSS_coxphv2, aes(x = HR, y = Tumourtype)) + 
    geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") + 
    geom_errorbarh(aes(xmax = HR.confint.upper, xmin = HR.confint.lower), size = .5, height = 
                     .2, color = "gray50") +
    geom_point(size = 3.5, colour = pancan_molsubtype_DSS_coxphv2$sigcolours) +
    theme_bw()+
    theme(panel.grid.minor = element_blank()) +
    ylab("") +
    xlab("Hazard Ratio") +
    ggtitle("Disease-Specific Survival")
  pancan_molsubtype_DSS_coxph_HR

  
  #disease-free survival
  pancan_molsubtype_coxph<-list()
  p.value<-matrix(NA, ncol =1, nrow = length(tumourtypesDFI))
  wald.test<-matrix(NA, ncol =1, nrow = length(tumourtypesDFI))
  beta<-matrix(NA, ncol =1, nrow = length(tumourtypesDFI))
  HR<-matrix(NA, ncol =1, nrow = length(tumourtypesDFI))
  HR.confint.lower<-matrix(NA, ncol =1, nrow = length(tumourtypesDFI))
  HR.confint.upper<-matrix(NA, ncol =1, nrow = length(tumourtypesDFI))
  
  for (i in 1:length(tumourrowsDFI)){
    
    pancan_molsubtype_coxph[[i]]<-coxph(Surv(DFI.time, DFI) ~ Distmethod_cluster , data = pancan_rnaseqlTminz_survsig_scoredf4DFI[tumourrowsDFI[[i]],])
    p.value[i]<-summary(pancan_molsubtype_coxph[[i]])$wald["pvalue"]
    wald.test[i]<-summary(pancan_molsubtype_coxph[[i]])$wald["test"]
    beta[i]<-summary(pancan_molsubtype_coxph[[i]])$coef[nrow(summary(pancan_molsubtype_coxph[[i]])$conf.int),1] #to take K3 only#coeficient beta
    HR[i] <-summary(pancan_molsubtype_coxph[[i]])$coef[nrow(summary(pancan_molsubtype_coxph[[i]])$conf.int),2] #to take K3 only #HR
    HR.confint.lower[i] <- summary(pancan_molsubtype_coxph[[i]])$conf.int[nrow(summary(pancan_molsubtype_coxph[[i]])$conf.int),"lower .95"] #to take K3 only
    HR.confint.upper[i] <- summary(pancan_molsubtype_coxph[[i]])$conf.int[nrow(summary(pancan_molsubtype_coxph[[i]])$conf.int),"upper .95"] #to take K3 only
  }
  pancan_molsubtype_DFS_coxph<-data.frame(Tumourtype = tumourtypesDFI,
                                          pval = p.value,
                                          waldtest = wald.test,
                                          beta = beta,
                                          HR = HR,
                                          HR.confint.lower = HR.confint.lower,
                                          HR.confint.upper = HR.confint.upper)
  rownames(pancan_molsubtype_DFS_coxph) <-tumourtypesDFI
  
  #Plot the HR
   pancan_molsubtype_DFS_coxph_HR <- ggplot(reorder(pancan_molsubtype_DFS_coxph, -HR), aes(x = HR, y = Tumourtype)) + 
    geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") + 
    geom_errorbarh(aes(xmax = HR.confint.upper, xmin = HR.confint.lower), size = .5, height = 
                     .2, color = "gray50") +
    geom_point(size = 3.5, color ="red2") +
    theme_bw()+
    theme(panel.grid.minor = element_blank()) +
    ylab("") +
    xlab("Hazard Ratio") +
    ggtitle("K3 Molecular Subtype Disease-Free Survival")
  pancan_molsubtype_DFS_coxph_HR
  
  
  #Progression -free Survival
  pancan_molsubtype_coxph<-list()
  p.value<-matrix(NA, ncol =1, nrow = length(tumourtypesPFI))
  wald.test<-matrix(NA, ncol =1, nrow = length(tumourtypesPFI))
  beta<-matrix(NA, ncol =1, nrow = length(tumourtypesPFI))
  HR<-matrix(NA, ncol =1, nrow = length(tumourtypesPFI))
  HR.confint.lower<-matrix(NA, ncol =1, nrow = length(tumourtypesPFI))
  HR.confint.upper<-matrix(NA, ncol =1, nrow = length(tumourtypesPFI))
  
  for (i in 1:length(tumourrowsPFI)){
    
    pancan_molsubtype_coxph[[i]]<-coxph(Surv(PFI.time, PFI) ~ Distmethod_cluster , data = pancan_rnaseqlTminz_survsig_scoredf4PFI[tumourrowsPFI[[i]],])
    p.value[i]<-summary(pancan_molsubtype_coxph[[i]])$wald["pvalue"]
    wald.test[i]<-summary(pancan_molsubtype_coxph[[i]])$wald["test"]
    beta[i]<-summary(pancan_molsubtype_coxph[[i]])$coef[nrow(summary(pancan_molsubtype_coxph[[i]])$conf.int),1] #to take K3 only#coeficient beta
    HR[i] <-summary(pancan_molsubtype_coxph[[i]])$coef[nrow(summary(pancan_molsubtype_coxph[[i]])$conf.int),2] #to take K3 only #HR
    HR.confint.lower[i] <- summary(pancan_molsubtype_coxph[[i]])$conf.int[nrow(summary(pancan_molsubtype_coxph[[i]])$conf.int),"lower .95"] #to take K3 only
    HR.confint.upper[i] <- summary(pancan_molsubtype_coxph[[i]])$conf.int[nrow(summary(pancan_molsubtype_coxph[[i]])$conf.int),"upper .95"] #to take K3 only
  }
  pancan_molsubtype_PFS_coxph<-data.frame(Tumourtype = tumourtypesPFI,
                                          pval = p.value,
                                          waldtest = wald.test,
                                          beta = beta,
                                          HR = HR,
                                          HR.confint.lower = HR.confint.lower,
                                          HR.confint.upper = HR.confint.upper)
  rownames(pancan_molsubtype_PFS_coxph) <-tumourtypesPFI
  
  #Plot the HR for PFS
  pancan_molsubtype_PFS_coxphv2<-pancan_molsubtype_PFS_coxph[-c(11,26),]
  pancan_molsubtype_PFS_coxphv2$Tumourtype<-reorder(pancan_molsubtype_PFS_coxphv2$Tumourtype, pancan_molsubtype_PFS_coxphv2$HR)
  pancan_molsubtype_PFS_coxphv2$sigcolours<-ifelse(pancan_molsubtype_PFS_coxphv2$pval<=0.055, "red2", "royalblue")
  
  
  pancan_molsubtype_PFS_coxph_HR <- ggplot(pancan_molsubtype_PFS_coxphv2, aes(x = HR, y = Tumourtype)) + 
    geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") + 
    geom_errorbarh(aes(xmax = HR.confint.upper, xmin = HR.confint.lower), size = .5, height = 
                     .2, color = "gray50") +
    geom_point(size = 3.5, colour = pancan_molsubtype_PFS_coxphv2$sigcolours) +
    theme_bw()+
    theme(panel.grid.minor = element_blank()) +
    ylab("") +
    xlab("Hazard Ratio") +
    ggtitle("Progression-Free Survival")
  pancan_molsubtype_PFS_coxph_HR