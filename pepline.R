####################################
####Prepared for Moldesign
#####################################
#####################################
#### Step 1. Load data.
#####################################
library(minfi)
library(ChAMP)
library(ENmix)
library(DMRcate)

#####load data
dataDir<-IDAT_DIR
targets <- read.metharray.sheet(dataDir)
rgSet <- read.metharray.exp(targets = targets,extended = TRUE,force = FALSE)

############################
########Step 2: sex check
###########################

sampleNames(rgSet)=rgSet[[1]]
mset<-preprocessRaw(rgSet)
GMSet<-mapToGenome(mset)
sex<-getSex(GMSet, cutoff = -2)
sex_out<-paste(RESULT_QC_DIR,"/sex_check_results.txt",sep="")
write.table(sex,sex_out,sep="\t")
#####################################
#### Step 3. Remove low quality samples according to detection P value.
#####################################

mset<-preprocessRaw(rgSet)
gmSet<-mapToGenome(mset)

tmp = getBeta(gmSet, "Illumina",offset = 100)

detP <- detectionP(rgSet)
detP<-detP[rownames(tmp),]
detPcut = DETECTION_PVALUE_CUTOFF

tmp[detP >= detPcut] <- NA 
numfail <- matrix(colMeans(is.na(tmp)))
rownames(numfail) <- colnames(detP)
colnames(numfail) <- "Failed CpG Fraction"
outfile = paste(RESULT_QC_DIR, "/Failed_CpG_fraction.txt", sep = "")
write.table(numfail,outfile,sep = "\t")

SampleCutoff= SAMPLE_CUT_OFF

RemainSample <- which(numfail < SAMPLE_CUT_OFF) 
rgSet <- rgSet[,RemainSample]
gmSet <- gmSet[,RemainSample]
detP <- detP[,RemainSample]
tmp <- tmp[,RemainSample]



#####################################
##### Step 4. Background, dye-bias correction and functional quantile normalization
#####################################

## Beta before Functional Normalization

message("\n Step 3: Start Functional Normalization!\n")

GmSet_FunNorm<-preprocessFunnorm(rgSet, nPCs = 2, bgCorr = TRUE, dyeCorr = TRUE, keepCN = TRUE, ratioConvert = FALSE, verbose = TRUE)


###########################
#### Filter_bad_snp_mulitHit_ChrXY_probes.
#####################################
#### Step 5. filter probes with >0.01 detP in more than 5% samples.
#####################################
ProbeCutoff = PROBE_CUT_OFF 
GmSet<-GmSet_FunNorm[rownames(gmSet),]
GmSet_detP = GmSet[rowSums(is.na(tmp)) <= ProbeCutoff*ncol(detP),]
tmp_detP = getBeta(GmSet_detP, "Illumina",offset = 100)

#####################################
###### Step 6. filter out probes with less than 3 beads in at least 5% of samples per probe.
#####################################

### mybeadcount function from ChAMP R package https://bioconductor.org/packages/release/bioc/vignettes/ChAMP/inst/doc/ChAMP.html

mybeadcount <- function(x)
{
  #select out bead count data frame
  getNBeads(x) -> nb
  locusNames <- getManifestInfo(x, "locusNames")
  bc_temp <- matrix(NA_real_, ncol = ncol(x), nrow = length(locusNames),
                    dimnames = list(locusNames, sampleNames(x)))
  
  TypeII.Name <- getProbeInfo(x, type = "II")$Name
  bc_temp[TypeII.Name, ] <- nb[getProbeInfo(x, type = "II")$AddressA,]
  
  TypeI <- getProbeInfo(x, type = "I")
  
  bc_temp->bcB
  bc_temp->bcA    
  
  bcB[TypeI$Name, ] <- nb[TypeI$AddressB,]
  bcA[TypeI$Name, ] <- nb[TypeI$AddressA,]
  
  which(bcB< BEAD_COUNTS) -> bcB3
  which(bcA< BEAD_COUNTS) -> bcA3
  bcA->bcA2
  bcB->bcB2
  bcA2[bcA3]<-NA
  bcA2[bcB3]<-NA
  
  data.frame(bcA2)->bc
  bc
}

bc = mybeadcount(rgSet)
beadCutoff = BEAD_CUTOFF 
bc2 = bc[rowSums(is.na(bc)) < beadCutoff*(ncol(bc)),]

remain_cpg<-row.names(bc2)
GmSet_BC = GmSet_detP[featureNames(GmSet_detP) %in% remain_cpg,]
tmp_BC = getBeta(GmSet_BC, "Illumina",offset = 100)

#####################################
#### Step 7. remove SNP related CpGs
#####################################

data(EPIC.manifest.hg19)
maskname <- rownames(EPIC.manifest.hg19)[which(EPIC.manifest.hg19$MASK_general == TRUE)]
GmSet_SNP = GmSet_BC[!featureNames(GmSet_BC) %in% maskname,]
tmp_SNP = getBeta(GmSet_SNP, "Illumina",offset = 100)

#####################################
##### Step 8. filter multi-hit probes
#####################################
data(multi.hit)
GmSet_multi <-GmSet_SNP[!featureNames(GmSet_SNP) %in% multi.hit$TargetID,]
tmp_multi = getBeta(GmSet_multi, "Illumina",offset = 100)

#####################################
#### Step 9. filter out all probes located in chromosome X and Y.
#####################################
annotation <- getAnnotation(GmSet_multi)
sex_probe <- rownames(annotation)[annotation$chr %in% c("chrX", "chrY")]
GmSet_XY <- GmSet_multi[!featureNames(GmSet_multi) %in% sex_probe, ]
tmp_XY<-getBeta(GmSet_XY, "Illumina",offset = 100)


#####################################
#### EWAS
#####################################
library(minfi)
library(ggExtra)
#Corrective covariable
d = read.csv("immune_pure.csv",header = TRUE,row.names = 1,check.names = F)
sample_clin_feature = read.csv("clinical_messagecsv",header = TRUE, check.names = F)
methy= d[,as.character(sample_clin_feature$ID)]

methy = as.matrix(methy)
lm_model <- lm(t(methy) ~ as.factor(sample_clin_feature$gender)+sample_clin_feature$age+
                 sample_clin_feature$leveledu+sample_clin_feature$firstepi+
                 sample_clin_feature$monthesofdur)
tmp_lm <- t(lm_model$res)+rowMeans(methy)


dmp_pos <- dmpFinder(tmp_lm,pheno = sample_clin_feature$positive_diff,type = "continuous")
write.csv(dmp_pos,"positiveS0.csv")


dmp_neg <- dmpFinder(tmp_lm,pheno = sample_clin_feature$negativeS0,type = "continuous")
write.csv(dmp_neg,"negativeS0.csv")

dmp_gen <- dmpFinder(tmp_lm,pheno = sample_clin_feature$generalS0,type = "continuous")
write.csv(dmp_gen,"generalS0.csv")

dmp_total <- dmpFinder(tmp_lm,pheno = sample_clin_feature$w0panssts,type = "continuous")
write.csv(dmp_panss,"totalS0.csv")
#drug response
m_model2 <- lm(t(methy) ~ as.factor(sample_clin_feature$gender)+sample_clin_feature$age+
                 sample_clin_feature$leveledu+sample_clin_feature$firstepi+
                 sample_clin_feature$monthesofdur+sample_clin_feature$Dose_conversion+
                 sample_clin_feature$w0panssts)
tmp_lm2 <- t(lm_model2$res)+rowMeans(methy)

dmp_diff_Effectiveness <- dmpFinder(tmp_lm2,pheno = sample_clin_feature$Effectiveness,type = "continuous")
write.csv(dmp_diff_Effectiveness,"response.csv")


#####################################
#### Random forest model training
#####################################


library("mlr3")
library("mlr3viz")
library("mlr3verse")
library("ggplot2")
dmp_venn = read.csv("venn_dmp_all_cg.csv")
sample_clin_feature_cls = read.csv("clinical_message.csv",check.names = "F")
methy_veen = methy[as.character(dmp_venn$ID),]
methy_veen = as.data.frame(t(methy_veen))
methy_veen = cbind(sample_clin_feature_cls,methy_veen)
methy_veen = methy_veen[,-1]


CVgroup <- function(k,datasize){
  set.seed(10)
  cvlist <- list()
  n <- rep(1:k,ceiling(datasize/k))[1:datasize]    #将数据分成K份，并生成的完成数据集n
  temp <- sample(n,datasize)   #把n打乱
  x <- 1:k
  dataseq <- 1:datasize
  cvlist <- lapply(x,function(x) dataseq[temp==x])  #dataseq中随机生成k个随机有序数据列
  return(cvlist)
}


a = CVgroup(5,nrow(methy_veen))


methy_veen_sub1 = methy_veen[rownames(methy_veen)[a[[1]]],]
methy_veen_sub2 = methy_veen[rownames(methy_veen)[a[[2]]],]
methy_veen_sub3 = methy_veen[rownames(methy_veen)[a[[3]]],]
methy_veen_sub4 = methy_veen[rownames(methy_veen)[a[[4]]],]
methy_veen_sub5 = methy_veen[rownames(methy_veen)[a[[5]]],]
methy_veen_train1 = rbind(rbind(rbind(methy_veen_sub2,methy_veen_sub3),methy_veen_sub4),methy_veen_sub5)
methy_veen_train2 = rbind(rbind(rbind(methy_veen_sub1,methy_veen_sub3),methy_veen_sub4),methy_veen_sub5)
methy_veen_train3 = rbind(rbind(rbind(methy_veen_sub1,methy_veen_sub2),methy_veen_sub4),methy_veen_sub5)
methy_veen_train4 = rbind(rbind(rbind(methy_veen_sub1,methy_veen_sub2),methy_veen_sub3),methy_veen_sub5)
methy_veen_train5 = rbind(rbind(rbind(methy_veen_sub1,methy_veen_sub2),methy_veen_sub3),methy_veen_sub4)

# fold1
task_data = as_task_classif(methy_veen, target = "Effectiveness", id = "train1")
task_data_train1 = as_task_classif(methy_veen_train1, target = "Effectiveness", id = "train1")
task_data_test1 = as_task_classif(methy_veen_sub1, target = "Effectiveness", id = "train1")
learner_rf1 = lrn("classif.ranger",importance = "permutation",num.trees = 500, mtry=200)
learner_rf1$train(task_data_train1)
measure = msr("classif.acc")
learner_rf1$predict(task_data_test1)$score(measure)
import1 = learner_rf1$importance()

# fold2
task_data_train2 = as_task_classif(methy_veen_train2, target = "Effectiveness", id = "train2")
task_data_test2 = as_task_classif(methy_veen_sub2, target = "Effectiveness", id = "train2")
learner_rf2 = lrn("classif.ranger",importance = "permutation",num.trees = 500, mtry=200)
learner_rf2$train(task_data_train2)
measure = msr("classif.acc")
learner_rf2$predict(task_data_test2)$score(measure)
import2 = learner_rf2$importance()

# fold3
task_data_train3 = as_task_classif(methy_veen_train3, target = "Effectiveness", id = "bbb")
task_data_test3 = as_task_classif(methy_veen_sub3, target = "Effectiveness", id = "bbb")
learner_rf3 = lrn("classif.ranger",importance = "permutation",num.trees = 500, mtry=100)
learner_rf3$train(task_data_train3)
measure = msr("classif.acc")
learner_rf3$predict(task_data_test3)$score(measure)
import3 = learner_rf3$importance()

# fold4
task_data_train4 = as_task_classif(methy_veen_train4, target = "Effectiveness", id = "bbb")
task_data_test4 = as_task_classif(methy_veen_sub4, target = "Effectiveness", id = "bbb")
learner_rf4 = lrn("classif.ranger",importance = "permutation",num.trees = 500, mtry=100)
learner_rf4$train(task_data_train4)
measure = msr("classif.acc")
learner_rf4$predict(task_data_test4)$score(measure)
import4 = learner_rf4$importance()

# fold5
task_data_train5 = as_task_classif(methy_veen_train5, target = "Effectiveness", id = "bbb")
task_data_test5 = as_task_classif(methy_veen_sub5, target = "Effectiveness", id = "bbb")
learner_rf5 = lrn("classif.ranger",importance = "permutation",num.trees = 500, mtry=100)
learner_rf5$train(task_data_train5)
measure = msr("classif.acc")
learner_rf5$predict(task_data_test5)$score(measure)
import5 = learner_rf5$importance()

##############
import = (import1+import2+import3+import4+import5)/5
importance = as.data.frame(import, keep.rownames = TRUE)
colnames(importance) = "Importance"
importance$Feature = rownames(importance)
rownames(importance) = seq(1:nrow(importance))
write.csv(importance,"importance.csv")


#Rf MODEL IMPORTANCE RANK GSEA
gene = read.csv("cgimpor_rank.csv")
colnames(gene)[2] = "SYMBOL"
gene = gene [!duplicated(gene$SYMBOL),]
geneList_df = bitr(gene$SYMBOL,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
geneList_df = merge(geneList_df,gene,by="SYMBOL")
geneList = as.numeric(geneList_df$Importance)
names(geneList) = geneList_df$ENTREZID
geneList=sort(geneList,decreasing = TRUE) 
Go_gseresult <- gseGO(geneList, 'org.Hs.eg.db', keyType = "ENTREZID", ont="all", nPerm = 100, minGSSize = 10, maxGSSize = 100, pvalueCutoff=1)
write.csv(as.data.frame(Go_gseresult),"GO.csv")
ridgeplot(Go_gseresult,50,fill = "pvalue", core_enrichment = TRUE)
KEGG_gseresult <- gseKEGG(geneList, nPerm = 100, minGSSize = 10, maxGSSize = 100, pvalueCutoff=1)
write.csv(as.data.frame(KEGG_gseresult),"KEGG.csv")