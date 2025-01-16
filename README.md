# Multiomics
#Methylation
library(knitr)
library(limma)
library(minfi)

#if(!require("BiocManager",quietly = TRUE))
#  install.packages("BiocManager")

###("IlluminaHumanMethylation450kanno.ilmn12.hg19")

###BiocManager::install("IlluminaHumanMethylationEPICv2anno.20a1.hg38",force = TRUE)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

###BiocManager::install("IlluminaHumanMethylationEPICv2manifest",force = TRUE)
library(IlluminaHumanMethylationEPICv2manifest)
library(RColorBrewer)

BiocManager::install("missMethy")

library(missMethy) ###not done
library(minfiData)

###BiocManager::install("minfiData")

library(Gviz)
###BiocManager::install("Gviz")

library(DMRcate)
###BiocManager::install("DMRcate")
library(stringr)


###dataDirectory <- system.file("", package = "methylationArrayAnalysis")
list.files("E:/Longitudinal_multiomics_space/data_yufei/idat", recursive =  TRUE)

###列出文件夹下所有文件
###list.files("dataDirectory", recursive = TRUE)
########################################################################################################################
###925K注释信息
help(package="IlluminaHumanMethylationEPICv2anno.20a1.hg38")
ann925k <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
head(ann925k)

###读取rawData
targets <- read.metharray.sheet("D:/2025/Longitudinal_multiomics_space/data_yufei", pattern = "sample_data.csv")
head(targets)

dataDirectory <- system.file("D:/2025/Longitudinal_multiomics_space/data_yufei/idat_yufei", package = "ChAMPdata")
list.files("E:/Longitudinal_multiomics_space/data_yufei/idat_yufei", recursive =  TRUE)


targets$ID <- paste(targets$Sample_Group,targets$Sample_Name,sep=".")

head(targets)
###myload <- champ.load("E:/Longitudinal_multiomics_space/data_yufei/idat_yufei", arraytype="EPICv2")
###计算探针的P值  
rgSet <- read.metharray.exp(targets = targets)
rgSet
head(rgSet)

###给予样本细节名字
targets$ID <- paste(targets$Sample_Group,targets$Sample_Name,sep=".")
sampleNames(rgSet) <- targets$Slide
rgSet
##############################################################################################################

###Quality control

###计算P-value  
detP <- detectionP(rgSet)
head(detP)
#write.csv(detP,file = "detP_out.csv")



# examine mean detection p-values across all samples to identify any failed samples
pal <- brewer.pal(8,"Dark2")
par(omi =c(0.1,0,0,0),#调整图形外边距
    mar =c(5,5,4,3),#调整图形内边距
    mfrow=c(1,2))
barplot(colMeans(detP), col=pal[factor(targets$Sample_Group)], las=1, 
        cex.names=0.8, ylab="")
###text(x=-1, y=0.003, srt = 90, xpd=TRUE, labels="Mean detection p-values") 



abline(h=0.05,col="red")
legend("topright", legend=levels(factor(targets$Sample_Group)), fill=pal,
       bg="white",cex = 0.8)

barplot(colMeans(detP), col=pal[factor(targets$Sample_Group)], las=1, 
        cex.names=0.8, ylim=c(0,0.002), ylab="")
abline(h=0.05,col="red")
legend("topright", legend=levels(factor(targets$Sample_Group)), fill=pal, 
       bg="white",cex = 0.5)

qcReport(rgSet, sampNames=targets$ID, sampGroups=targets$Sample_Group, 
         pdf="qcReport.pdf")

pd <- pData(rgSet)
keep = apply(detP,1,function(x){all(x < 0.01)})
head(keep)



###去除低质量的样品
keep <- colMeans(detP) < 0.05
rgSet <- rgSet[,keep]
rgSet
###去除目标数据中的低质量样品
targets <- targets[keep,]
targets[,1:5]
###去除检测到的p值的低质量样品
detP <- detP[,keep]
dim(detP)
###
#> dim(detP)
#[1] 936990    156


#######Normalisation
# normalize the data; this results in a GenomicRatioSet object

mSetSq <- preprocessQuantile(rgSet) 
mSetSq1 <- preprocessFunnorm(rgSet) 

# create a MethylSet object from the raw data for plotting
mSetRaw <- preprocessRaw(rgSet)

# visualise what the data looks like before and after normalisation
par(mfrow=c(1,2))
densityPlot(rgSet, sampGroups=targets$Sample_Group,main="Raw", legend=FALSE)
legend("top", legend = levels(factor(targets$Sample_Group)), 
       text.col=brewer.pal(8,"Dark2"),cex=0.8)
densityPlot(getBeta(mSetSq), sampGroups=targets$Sample_Group,
            main="Normalized", legend=FALSE)
legend("top", legend = levels(factor(targets$Sample_Group)), 
       text.col=brewer.pal(8,"Dark2"),cex=0.8)
############################################################# Data exploration

# MDS plots to look at largest sources of variation
library(limma)

par(mfrow=c(1,2))
plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)])
###legend("topright", legend=levels(factor(targets$Sample_Group)), text.col=pal,
###      cex=0.6, bty = "n")

plotMDS(getM(mSetSq), top=1000, gene.selection="common",  
        col=pal[factor(targets$Sample_Name)])
###legend("topright", legend=levels(factor(targets$Sample_Name)), text.col=pal,
###        cex=0.6,bty = "n")

# Examine higher dimensions to look at other sources of variation
par(mfrow=c(1,3))
plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)], dim=c(1,3))
###legend("topleft", legend=levels(factor(targets$Sample_Group)), text.col=pal, 
###       cex=0.7, bg="white")

plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)], dim=c(2,3))
###legend("topleft", legend=levels(factor(targets$Sample_Group)), text.col=pal,
###      cex=0.7,bg="white")

plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)], dim=c(3,4))
###legend("topleft", legend=levels(factor(targets$Sample_Group)), text.col=pal,
###       cex=0.7, bg="white")



################################################################################
###Filtering
# ensure probes are in the same order in the mSetSq and detP objects
detP <- detP[match(featureNames(mSetSq),rownames(detP)),] 

# remove any probes that have failed in one or more samples
keep <- rowSums(detP < 0.01) == ncol(mSetSq) 
table(keep)

mSetSqFlt <- mSetSq[keep,]
mSetSqFlt
# if your data includes males and females, remove probes on the sex chromosomes
keep <- !(featureNames(mSetSqFlt) %in% ann925k$Name[ann925k$chr %in% 
                                                      c("chrX","chrY")])
table(keep)
mSetSqFlt <- mSetSqFlt[keep,]

# remove probes with SNPs at CpG site
mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
mSetSqFlt

# exclude cross reactive probes 
###xReactiveProbes <- read.csv(file=paste(dataDirectory,
###                                       "48639-non-specific-probes-Illumina450k.csv",
###                                       sep="/"), stringsAsFactors=FALSE)
###keep <- !(featureNames(mSetSqFlt) %in% xReactiveProbes$TargetID)
###table(keep)
mSetSqFlt <- mSetSqFlt[keep,] 
mSetSqFlt

par(mfrow=c(1,2))
plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Group)], cex=0.8)
###legend("right", legend=levels(factor(targets$Sample_Group)), text.col=pal,cex=0.65, bg="white")

plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Name)])
###legend("right", legend=levels(factor(targets$Sample_Name)), text.col=pal,cex=0.7, bg="white")
par(mfrow=c(1,3))
# Examine higher dimensions to look at other sources of variation
plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Name)], dim=c(1,3))

plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Name)], dim=c(2,3))


plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common", 
        col=pal[factor(targets$Sample_Name)], dim=c(3,4))

##################################################################
# calculate M-values for statistical analysis
mVals <- getM(mSetSqFlt)
head(mVals[,1:5])
bVals <- getBeta(mSetSqFlt)
head(bVals[,1:5])
par(mfrow=c(1,2))
densityPlot(bVals, sampGroups=targets$Sample_Group, main="Beta values", 
            legend=FALSE, xlab="Beta values", cex.main=0.6, cex.sub=0.6, cex.lab=0.6, cex.axis=0.6)
densityPlot(mVals, sampGroups=targets$Sample_Group, main="M-values", 
            legend=FALSE, xlab="M values", cex.main=0.6, cex.sub=0.6, cex.lab=0.6, cex.axis=0.6)


######################################################################
###Probe-wise differential methylation analysis
# this is the factor of interest
time <- factor(targets$Sample_Group)
# this is the individual effect that we need to account for
individual <- factor(targets$Sample_ID) 

# use the above to create a design matrix
design <- model.matrix(~0+time+individual, data=targets)

#write.csv(design,file = "design.csv", row.names = FALSE)

colnames(design) <- c(levels(time),levels(individual)[-1])

# fit the linear model 
fit <- lmFit(mVals, design)
# create a contrast matrix for specific comparisons
contMatrix <- makeContrasts(time2-time4,
                            time4-time6,
                            time6-time7,
                            time7-time8,
                            time8-time11,
                            levels=design)
contMatrix

# fit the contrasts
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)

# look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit2))
##############################################################################
# get the table of results for the first contrast (time2 - time4)
ann925kSub <- ann925k[match(rownames(mVals),ann925k$Name),
                      c(1:4,12:19,24:ncol(ann925k))]
DMPs <- topTable(fit2, num=Inf, coef=1, genelist=ann925kSub)

write.csv(DMPs, file = "DMPs.csv", row.names = FALSE)
head(DMPs)

##########################################################################
# plot the top 4 most significantly differentially methylated CpGs 
par(mfrow=c(2,2))
sapply(rownames(DMPs)[1:4], function(cpg){
  plotCpg(bVals, cpg=cpg, pheno=targets$Sample_Group, ylab = "Beta values")
})


###################################################################




